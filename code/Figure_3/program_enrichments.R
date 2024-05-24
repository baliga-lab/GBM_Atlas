# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2024
# Email: sturkarslan@systemsbiology.org
#
# Date: 2024-04-22
#
# Script Name: program_enrichment
#
# Script Description: Given the program information from MINER, build program enrichment sumamry plots and output
#
# Script version: v1
# --------------------------------------------

library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ReactomePA)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(jsonlite)
library(ggsurvfit)
library(survival)
library(survRM2)
library(gghighlight)
library(cowplot)
library("grid")
library("ggplotify")
library(ComplexHeatmap)
library(visNetwork)
library(tidyverse)
# Get todays date for date stamp
today <- Sys.Date()
date_stamp <- format(today, format="%b_%d_%Y")

source("code/program_functional_enrichment.R")
source("code/load_program_survival.R")
source("code/plot_program_survival.R")
source("code/plot_program_cmflows.R")
source("code/plot_mutation_cmflows.R")

### Get expression pattern for selected Programs
# Load model files
overExpressedMembersMatrix <- read_csv("../GBM-Model-052022/overExpressedMembers.csv") |>
  column_to_rownames("...1")

underExpressedMembersMatrix <- read_csv("../GBM-Model-052022/underExpressedMembers.csv") |>
  column_to_rownames("...1")
# get the matrix file
dfr = (overExpressedMembersMatrix-underExpressedMembersMatrix) |>
  rownames_to_column("Regulon") |>
  mutate(Regulon = paste0("R-",Regulon))

## Process Programs for annotating heatmap rows
## We only show disease relevant programs in the heatmap
disease_programs <- read_csv("../GBM-Model-052022/CoxProportionalHazardsPrograms.csv") |>
  filter(`p-value` <= 0.05) |>
  dplyr::rename("Program_ID" = "...1")

## Load disease relevant regulons
disease_regulons <- read_csv("../GBM-Model-052022/diseaseRelevantRegulonsSummary.csv")

# Load programs ans states from model
programs_json = read_json(simplifyVector = T,"../GBM-Model-052022/transcriptional_programs.json")
states_json = read_json(simplifyVector = T, "../GBM-Model-052022/transcriptional_states.json")

# get ordering of the programs
df2 <- dfr |>
  arrange(factor(Regulon, levels=paste0("R-",as.vector(unlist(programs_json)))))

# get ordering of states
df3 <- df2 |>
  dplyr::select(as.vector(unlist(states_json)))

# unlist states dictionary
state_df <- data.frame()
for(i in 1:length(states_json)){
  tmp1 <- data.frame(states_json[i]) |>
    dplyr::rename("regulon" = 1) |>
    mutate(state = names(states_json[i]))
  state_df <- bind_rows(state_df, tmp1)
}

### State annotations
state_df <- state_df |>
  column_to_rownames("regulon")

# Final formatting fr annotations
state_df_annot <- state_df |>
  rownames_to_column("regulon") |>
  rownames_to_column("row_id") |>
  dplyr::select(state)

state_list = list()
for(i in unique(state_df_annot$state)){
  state_list[i] = list(which(state_df_annot$state == i))
}

## Load regulons disease programs
regulons_programs_genes <- read_csv(file="../GBM-Model-052022/regulonsProgramsGenesMapping.csv")

## Load drug database
drugs_mapped_to_network <- read_csv(file="../GBM-Model-052022/Drugs_Mapped_to_Network_Mar_03_2023.csv")

# Load CMFlows
cmflows_disease_genes <- read_csv("../GBM-Model-052022/diseaseRelevantCMFlowsGenes.csv")
cmflows_disease_pathways <- read_csv("../GBM-Model-052022/diseaseRelevantCMFlowsPathways.csv")
### Load miRNA cmflows
cmflows_disease_mirna <- read_csv("../GBM-Model-052022/filteredCausalResults_miRNA.csv")

## Load mutation data to filter mutations for 5 % of TCGA patients
## TCGA Gene Mutations
mutated_genes_df <- read_csv("/Volumes/omics4tb2/kkannan/mutations.GBM.TCGA.2019.05.01.csv") |>
  dplyr::rename("Mutation" = "...1") |>
  #dplyr::select(-Gene, -GeneID ) |>
  pivot_longer(!c(Mutation, Gene, GeneID), names_to = "Patient_ID", values_to = "Status")

## Mutation percents
mutation_percents_genes <- mutated_genes_df |>
  group_by(Mutation, Status) |>
  summarize(count = n()) |>
  filter(Status == 1) |>
  ungroup()|>
  mutate(percent = round((count*100)/247)) |>
  filter(percent >= 5)

cmflows_disease_genes_filtered <- cmflows_disease_genes |>
  filter(Mutation %in% mutation_percents_genes$Mutation)


cmflows_all <- bind_rows(cmflows_disease_genes, cmflows_disease_pathways)
# Load Cancer hallmarks
program_hallmarks <- read_csv("../GBM-Model-052022/cancerHallmarks_enrichments_in_programs.csv")
identifiers <- read_delim("../GbmMINER/data/identifier_mappings.txt", delim = "\t")

######## FUNCTIONS  ###########################

plot_program_summary <- function(in_program){

  ## Get disease relevant program regulons
  program_regulons <- regulons_programs_genes |>
    dplyr::filter(Programs == in_program ) |>
    dplyr::filter(Regulon_ID %in% disease_regulons$Regulon_ID) |>
    dplyr::pull(Regulon_ID) |>
    unique()

  program_genes <- regulons_programs_genes |>
    dplyr::filter(Regulon_ID %in% program_regulons) |>
    dplyr::pull(Gene) |>
    unique()

  ## Get expression matrix ordered by programs and states
  plot_df1 <- dfr |>
    arrange(factor(Regulon, levels=paste0("R-",as.vector(unlist(programs_json))))) |>
    dplyr::filter(Regulon %in% paste0("R-",program_regulons))|>
    dplyr::filter(Regulon %in% paste0("R-",disease_regulons$Regulon_ID))|>
    column_to_rownames("Regulon") |>
    dplyr::select(as.vector(unlist(states_json)))

  # funtion for drawing state annotations
  panel_fun = function(index, names) {
    grid.rect(gp = gpar(fill = "orange", col = "white"))
    grid.text(paste0(names), 0.5, 0.5, gp=gpar(fontsize = 9))
  }

  state_annotations <- HeatmapAnnotation(name = "States",
                                         State = anno_block(
                                           align_to = state_list,
                                           panel_fun = panel_fun,
                                           gp = gpar(fill = "orange", fontsize=10),
                                           which = "column",
                                         ),
                                         gap = unit(2, "mm"),
                                         which="column")

## Heatmap of program activity
 pp1 <-  Heatmap(plot_df1,#left_annotation=program_annotations,
          use_raster = T,
          heatmap_legend_param = list(
            title = "Activity"),
          #column_labels = ticks_df$label,
          cluster_row_slices = F,
          row_title = "Regulons",
          column_title = "Patients",
          row_title_side = "left",
          column_title_side = "bottom",
          row_dend_reorder = F,
          show_row_dend = F,
          raster_resize_mat = T,
          raster_magick_filter = F,
          show_column_dend = F,
          show_column_names = F,
          show_row_names = T,
          cluster_columns = F,
          cluster_rows = F,
          #column_split = state_df_annot,
          #split = prog_df_annot,
          #split=prog_df_disease_final,cluster_row_slices = F,
          col=colorRampPalette(c("blue", "white", "red"))(100),
          #left_annotation = program_annotations,
          top_annotation = state_annotations#, #columnAnnotation(df=state_df_annot,  gp = gpar(color = "yellow"))
          #bottom_annotation = patient_annotations
  )



 # pp1 <- pheatmap(plot_df1, scale = "none", cluster_rows = F, cluster_cols = F, annotation_col = col_annot,show_colnames = F,annotation_names_col = T, main = in_program, legend = F)

  # cmf_sum <- cmflows_all |>
  #   dplyr::filter(Regulon %in% program_regulons)

  ## Load Program survival info
  res1 <- plot_program_survival(in_program)

  ## Program Enrichments
  res2 <- program_functional_enrichment(in_program)

  ## CMFlows plots
  net_fig <- plot_program_cmflows(in_program)

  #### Reformat plots if they are missing
  if(length(res2$reactome_plot$data$ID) != 0){
    p1 = res2$reactome_plot
  } else {
    p1 = NULL
  }

  if(length(res2$hallmark_plot$data$ID) != 0){
    p2 = res2$hallmark_plot
  } else {
    p2 = NULL
  }

  if(length(res2$curated_plot$data$ID) != 0){
    p3 = res2$curated_plot
  } else {
    p3 = NULL
  }

  middle_row <- plot_grid(p1,p2, ncol=2, nrow=1)
  bottom_row <- plot_grid(p3, as.grob(res1), ncol=2, nrow=1)

  final_plot <- cowplot::plot_grid(as.grob(pp1), middle_row, bottom_row, ncol=1, nrow=3)

  cancer_hallmarks_filtered <- res2$cancer_hallmarks |>
    dplyr::filter(wangScore > 0.8) |>
    pull(hallmark)


  cat("Program: ", in_program,
      #" Mutations: ", length(unique(cmf_sum$MutationSymbol)),
      #" Regulators: ", length(unique(cmf_sum$RegulatorSymbol)),
     # " Regulons: ", length(unique(cmf_sum$Regulon)),
      " Genes: ", length(unique(program_genes)),
      " Cancer Hallmarks: ", length(cancer_hallmarks_filtered),
      "\n"
      )

  ## Collate all enrichment results
  all_enrichments <- bind_rows(res2$reactome_enrichment, res2$hallmark_enrichment, res2$curated_enrichment)

  #setwd("~/Documents/GitHub/GBM-Model-052022/")
  htmlwidgets::saveWidget(net_fig, file = paste0("../GBM-Model-052022/Program_Enrichment_Summaries/PR-",in_program, "_CNFlow_summary.html"),selfcontained = T,libdir = "~/Documents/GitHub/GBM-Model-052022/Program_Enrichment_Summaries/html_tmp",)

  ggsave(final_plot, file = paste0("../GBM-Model-052022/Program_Enrichment_Summaries/PR-",in_program, "_enrichment_summary.pdf"), device = "pdf", width = 11, height=11, units = "in")

  write_csv(all_enrichments, file = paste0("../GBM-Model-052022/Program_Enrichment_Summaries/PR-",in_program, "_enrichment_summary.csv"))

  return(final_plot)
}

