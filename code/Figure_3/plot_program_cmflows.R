# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2023
# Email: sturkarslan@systemsbiology.org
#
# Date: 2023-03-08
#
# Script Name: plot_program_cmflows
#
# Script Description: given the program plots summary of cmflows associated
#
# Script version: v1.0.0
# --------------------------------------------
library(tidyverse)
library(visNetwork)

## Data files
## Load regulons disease programs
regulons_programs_genes <- read_csv(file="../GBM-Model-052022/regulonsProgramsGenesMapping.csv",show_col_types = FALSE)

## Load disease relevant regulons
disease_regulons <- read_csv("../GBM-Model-052022/diseaseRelevantRegulonsSummary.csv",show_col_types = FALSE)

# Load CMFlows genes
cmflows_disease_genes <- read_csv("../GBM-Model-052022/diseaseRelevantCMFlowsGenes.csv",show_col_types = FALSE)
# Load CMFlows pathways
cmflows_disease_pathways <- read_csv("../GBM-Model-052022/diseaseRelevantCMFlowsPathways.csv",show_col_types = FALSE)
### Load miRNA cmflows
cmflows_disease_mirna <- read_csv("../GBM-Model-052022/filteredCausalResults_miRNA.csv",show_col_types = FALSE)

## Load mutation data to filter mutations for 5 % of TCGA patients
## TCGA Gene Mutations
mutated_genes_df <- read_csv("/Volumes/omics4tb2/kkannan/mutations.GBM.TCGA.2019.05.01.csv",show_col_types = FALSE) |>
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

## Load mutation pathway data to filter mutations for 5 % of TCGA patients
## TCGA Pathway Mutations
mutated_pathways_df <- read_csv("/Volumes/omics4tb2/kkannan/NonSilentPathwayMutationsGbmRNAAndMicroarray07242020.csv",show_col_types = FALSE) |>
  dplyr::rename("Mutation" = "...1") |>
  dplyr::select(where(~ any(. != 0))) |>
  #dplyr::select(-Gene, -GeneID ) |>
  pivot_longer(!c(Mutation), names_to = "Patient_ID", values_to = "Status")

## Calculate percentage of the  mutated pathways across TCGA cohort
mutation_percents_pathways <- mutated_pathways_df |>
  group_by(Mutation, Status) |>
  summarize(count = n()) |>
  filter(Status == 1) |>
  ungroup()|>
  mutate(percent = round((count*100)/174)) |>
  filter(percent >= 5)

# filter the cmflow pathways for 5%
cmflows_disease_pathways_filtered <- cmflows_disease_pathways |>
  filter(Mutation %in% mutation_percents_pathways$Mutation) |>
  #pull(Mutation) |>
  unique()

# Combine filtered cmflows
cmflows_all <- bind_rows(cmflows_disease_genes_filtered, cmflows_disease_pathways_filtered)

## Load drug database
drugs_mapped_to_network <- read_csv(file="../GBM-Model-052022/Drugs_Mapped_to_Network_Mar_03_2023.csv",show_col_types = FALSE)

# Load Cancer hallmarks
program_hallmarks <- read_csv("../GBM-Model-052022/cancerHallmarks_enrichments_in_programs.csv",show_col_types = FALSE)

## Load identfier mappings
identifiers <- read_delim("../GbmMINER/data/identifier_mappings.txt", delim = "\t",show_col_types = FALSE)


###############################################
### Given the program name plots CMFlows
plot_program_cmflows <- function(in_program){

  ## Get disease relevant program regulons
  program_regulons <- regulons_programs_genes |>
    dplyr::filter(Programs == in_program ) |>
    dplyr::filter(Regulon_ID %in% disease_regulons$Regulon_ID) |>
    dplyr::pull(Regulon_ID) |>
    unique()

  # Get Mutations targeting TFs from filtered CMflows regulons
  mutation_genes_count <- cmflows_disease_genes_filtered |>
    filter(Regulon %in% program_regulons) |>
    pull(MutationSymbol) |>
    unique() |>
    length()

  # Get Mutation Pathways targeting TFs
  mutation_pathways_count <- cmflows_disease_pathways_filtered |>
    filter(Regulon %in% program_regulons) |>
    pull(Pathway) |>
    unique() |>
    length()

  # Get Mutations targeting miRNAs
  mutation_genes_mirna_count <- cmflows_disease_genes_filtered |>
    filter(grepl("mir", RegulatorSymbol)) |>
    filter(Regulon %in% program_regulons) |>
    pull(MutationSymbol) |>
    unique() |>
    length()

  # Get Mutation Pathways targeting miRNAs
  mutation_pathways_mirna_count <- cmflows_disease_pathways_filtered |>
    filter(grepl("mir", RegulatorSymbol)) |>
    filter(Regulon %in% program_regulons) |>
    pull(Pathway) |>
    unique() |>
    length()

  # get CMFlow regulator count
  regulator_count <- cmflows_all |>
    filter(Regulon %in% program_regulons) |>
    pull(RegulatorSymbol) |>
    unique() |>
    length()

  # get miRNA Count
  mirna_count <- cmflows_all |>
    filter(grepl("mir", RegulatorSymbol)) |>
    filter(Regulon %in% program_regulons) |>
    pull(RegulatorSymbol) |>
    unique() |>
    length()

  regulon_list <- cmflows_all |>
    filter(Regulon %in% program_regulons) |>
    pull(Regulon) |>
    unique()

  gene_count <- regulons_programs_genes |>
    filter(Regulon_ID %in% regulon_list) |>
    pull(Gene) |>
    unique() |>
    length()

  program_count <- regulons_programs_genes |>
    filter(Regulon_ID %in% regulon_list) |>
    pull(Programs) |>
    unique() |>
    length()

  # Get drug flows
  drug_filter <- drugs_mapped_to_network |>
    filter(Programs == in_program) |>
    dplyr::filter(Regulon_ID %in% disease_regulons$Regulon_ID) |>
    filter(!is.na(maxPhaseCancer))

  drug_count <- drug_filter |>
    pull(Drug) |>
    unique() |>
    length()

  hallmark_count <- program_hallmarks |>
    filter(program == in_program) |>
    filter(wangScore > 0.8) |>
    pull(hallmark) |>
    unique() |>
    length()

  # Create nodes dataframe
  nodes <- data.frame(
    id = c(paste0("MG-",mutation_genes_count),
           paste0("MP-",mutation_pathways_count),
           paste0("T-",regulator_count),
           paste0("mir-",mirna_count),
           paste0("R-",length(regulon_list)),
           paste0("P-",program_count),
           paste0("G-",gene_count),
           paste0("H-",hallmark_count),
           paste0("D-",drug_count)
    ),
  label = c(mutation_genes_count,
            mutation_pathways_count,
            regulator_count,
            mirna_count,
            length(regulon_list),
            program_count,
            gene_count,
            hallmark_count,
            drug_count
  ),
  group = c(paste0("MG-",mutation_genes_count),
            paste0("MP-",mutation_pathways_count),
            paste0("T-",regulator_count),
            paste0("mir-",mirna_count),
            paste0("R-",length(regulon_list)),
            paste0("P-",program_count),
            paste0("G-",gene_count),
            paste0("H-",hallmark_count),
            paste0("D-",drug_count)
  ),
  level = c(1,1,2,2,3,4,5,6,7),
  #shape = c("triangleDown","triangle","diamond","square", "circle","dot","star","hexagon"),
  color = c("green","green","red","red","blue","orange","gray","violet","pink" )
  )


  # Create edges dataframe
  edges <- data.frame(
    from = c(paste0("MG-",mutation_genes_count),
             paste0("MP-",mutation_pathways_count),
             paste0("MG-",mutation_genes_count),
             paste0("MP-",mutation_pathways_count),
             paste0("T-",regulator_count),
             paste0("mir-",mirna_count),
             paste0("R-",length(regulon_list)),
             paste0("P-",program_count),
             paste0("G-",gene_count),
             paste0("H-",hallmark_count),
             paste0("D-",drug_count)
    ),
    to = c(paste0("T-",regulator_count),
           paste0("T-",regulator_count),
           paste0("mir-",mirna_count),
           paste0("mir-",mirna_count),
           paste0("R-",length(regulon_list)),
           paste0("R-",length(regulon_list)),
           paste0("P-",program_count),
           paste0("G-",gene_count),
           paste0("H-",hallmark_count),
           paste0("D-",drug_count),
           "end"
    )
  )

  # Convert node labels to character
  nodes$label <- as.character(nodes$label)

  net_fig <- visNetwork(nodes, edges, width = "400px",  height = "200px", main = paste0("CM Flows for Program ", in_program)) |>
    visGroups(groupname = paste0("MG-",mutation_genes_count),
              shape = "icon",
              icon = list(code = "f107",color="green")) |>
    visGroups(groupname = paste0("MP-",mutation_pathways_count),
              shape = "icon",
              icon = list(code = "f103", color="green")) |>
    visGroups(groupname = paste0("T-",regulator_count),
              shape = "icon",
              icon = list(code = "f0d7", color="red")) |>
    visGroups(groupname = paste0("mir-",mirna_count),
              shape = "icon",
              label = "miRNA",
              icon = list(code = "f0d8", color="red4")) |>
    visGroups(groupname = paste0("R-",length(regulon_list)),
              shape = "icon",
              icon = list(code = "f0c8", color="blue")) |>
    visGroups(groupname = paste0("P-",program_count),
              shape = "icon",
              icon = list(code = "f0c9", color="orange")) |>
    visGroups(groupname = paste0("G-",gene_count),
              shape = "dot",
              color="lightgray") |>
    #icon = list(code = "f110", color="lightgray")) |>
    visGroups(groupname = paste0("H-",hallmark_count),
              shape = "icon",
              icon = list(code = "f02e", color="purple")) |>
    visGroups(groupname = paste0("D-",drug_count),
              shape = "icon",
              icon = list(code = "f0a3", color="pink"))  %>%
    addFontAwesome() |>
    visNodes(font="20px",fixed = T) |>
    visEdges(arrows = list(to = list(enabled = TRUE, width=2))) %>%
    visHierarchicalLayout(direction = "LR",levelSeparation = 100)

  return(net_fig)
}
