# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2024
# Email: sturkarslan@systemsbiology.org
#
# Date: 2024-05-21
#
# Script Name: analyze_depMap_PRISM.R
#
# Script Description: Analysis of Cancer Dependency map GBM cell lines and their drug response in PRISM screen.
#
# Script version:
# --------------------------------------------

library(tidyverse)
library(ggplot2)
library(edgeR)
library(reshape2)
library(gridExtra)
library(fs)

# ## Load Gene mapping identifiers from MINER
identifiers <- read.delim("~/Documents/GitHub/GbmMINER/data/identifier_mappings.txt", sep="\t") |>
  filter(Source == "Entrez Gene ID" | Source == "Gene Name")

identifiers.f <- identifiers |>
  filter(Source == "Entrez Gene ID" | Source == "Gene Name")

# load ide conversion for omics profiles to cell lines
profile_ids <- read_csv("data/OmicsProfiles.csv") |>
  select(ProfileID,ModelID) |>
  unique()

# Load depmap expression data
depmap_exp <- read_csv("data/OmicsExpressionGenesExpectedCountProfile.csv") |>
  dplyr::rename("omics_profile_id" = "...1") |>
  left_join(profile_ids, by=c("omics_profile_id" = "ProfileID")) |>
  relocate(ModelID) |>
  select(-omics_profile_id) |>
  dplyr::rename("depmapId" = "ModelID")
  

# Load GBM Cell line data
depmap_lines <- read_csv("data/cell-line-selector.csv") |>
  filter(lineage3 == "Glioblastoma") |>
  filter(!depmapId %in% c("ACH-000328", "ACH-000455"))

# filter expression for gbm cell lines
depmap_gbm <- depmap_exp |>
  filter(depmapId %in% depmap_lines$depmapId) |>
  column_to_rownames("depmapId") |>
  t() |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  separate(gene, into = c("gene_name", "ensembl_id"), sep = " ") |>
  mutate(ensembl_id = sub("\\(","",ensembl_id)) |>
  mutate(ensembl_id = sub("\\)","",ensembl_id))

# remove NAs
depmap_gbm2 <- depmap_gbm |>
  #left_join(identifiers, by=c("gene_name" = "Name")) |>
  #relocate("Preferred_Name") |>
  #dplyr::rename("ensembl_id" = "Preferred_Name") |>
  dplyr::select( -gene_name) |>
  #group_by(ensembl_id) |>
  #filter(n()>1)
  filter(!is.na(ensembl_id))# |>
  #column_to_rownames("ensembl_id")

## Load TCGA Data for normalization
# load TCGA id mappings from broad firehose data download
tcga_geneid_mapping <- read_delim("/Users/serdarturkaslan/Downloads/gdac.broadinstitute.org_GBM.mRNAseq_Preprocess.Level_3.2016012800.0.0/GBM.uncv2.mRNAseq_geneid_transcriptid_mapping.txt", delim="\t") |>
  separate(1, into=c("symbol","entrez_id"), sep = "[|]",fill = "left", remove = F)

# combine tcga ids and identifers to match across gene_id, entrez_id and symbol
ids.combined <- left_join(tcga_geneid_mapping, identifiers.f, by=c("entrez_id" = "Name")) |>
  dplyr::rename("gene_id" = "HYBRIDIZATION R") |>
  dplyr::select(-`TCGA-02-0047-01`)

# Load TCGA GBM RNASeq data raw counts for the cohort
tcga.raw.RNASeq <- read_delim("/Users/serdarturkaslan/Downloads/gdac.broadinstitute.org_GBM.mRNAseq_Preprocess.Level_3.2016012800.0.0/GBM.uncv2.mRNAseq_raw_counts.txt", delim="\t") |>
  dplyr::rename("gene_id" = 1) |>
  left_join(ids.combined, by="gene_id") |>
  relocate(Preferred_Name) |>
  select(-gene_id, -symbol, -Source, -entrez_id) |>
  dplyr::rename("ensembl_id" = "Preferred_Name") |>
  filter(!is.na(ensembl_id))

# for each patient analyzed, loop and load the RSEM data
all_cases_data <- left_join(tcga.raw.RNASeq, depmap_gbm2, by="ensembl_id") |>
  filter(!ensembl_id %in% c("ENSG00000004866", "ENSG00000011454", "ENSG00000086205", "ENSG00000104064", "ENSG00000110675", "ENSG00000111850", "ENSG00000116721", "ENSG00000127603", "ENSG00000137871", "ENSG00000143226", "ENSG00000157326", "ENSG00000169621", "ENSG00000188649", "ENSG00000189064", "ENSG00000196233", "ENSG00000205777", "ENSG00000223802", "ENSG00000224659", "ENSG00000276725")) |>
  column_to_rownames("ensembl_id") |>
  drop_na()

# perform TMM normalization of patient data
comb_TMM <- DGEList(counts = all_cases_data)
comb_TMM <- calcNormFactors(comb_TMM, method = "TMM")
comb_TMM <- cpm(comb_TMM,log = T)

# convert to dataframe for plotting
comb_TMM_df <- as.data.frame(comb_TMM) |>
  rownames_to_column("gene_id")

# format data for plotting
comb_plot_data <- melt(comb_TMM_df, id.vars = "gene_id") |>
  mutate(type = if_else(grepl("TCGA", variable), "TCGA","Patient"))

# draw a boxplot for TMM and write to pdf
pcomb <- ggplot(comb_plot_data, aes(x=variable, y=value, group=variable, fill=type, color=type))
pcomb <- pcomb + geom_boxplot()
pcomb <- pcomb + labs(title=paste0("TMM normalization with cohort RNASeq"),
                      x = "Patients", y="TMM Normalized Counts (Counts per million)")
pcomb <- pcomb + theme(axis.text.x = element_blank())
pcomb

write_csv(comb_TMM_df, file="data/depmap_expression_TMM_with_cohort.csv")





## Load DCRA for each cell line by reading analyze patient output
cell_ids <- dir_ls("output", recurse = FALSE, type = "directory") |>
as_tibble() |>
separate(value, into=c("output","cell_id"), sep="/") |>
pull(cell_id)



out_all <- tibble()
count = 0
for(cell_id in cell_ids){
  count = count + 1
  path_res <- paste0("output/",cell_id, "/",cell_id,"/",cell_id,"_analyze_patient_output.csv")
  
  # read file
  out1 <- read_csv(path_res) |>
    select(Drug, DrugConstrainedRegulonActivity) |>
    dplyr::rename({{cell_id}} := "DrugConstrainedRegulonActivity") |>
    unique()
  
  if(count < 2){
    out_all <- out1
  } else {
    out_all <- left_join(out_all, out1, by="Drug")
  }
 
}

# Convert it to longer format and change responder assignment
depmap_dcra <- out_all |>
  pivot_longer(!starts_with("Drug"), names_to = "cell_id", values_to = "dcra") |>
  mutate(response = if_else(dcra >= 0.25, "responder", "non-responder"))

# Drug screen Log Fold Change data
drug_screen <- read_csv("data/Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv") |>
  pivot_longer(cols = !starts_with("...1"), names_to = "cell_id", values_to = "lfc")

# Compound metadata
comp_meta <- read_csv("data/Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv") |>
  select(Drug.Name, IDs) |>
  unique()

# Drug final combined list with drug names
drug_final <- drug_screen |>
  inner_join(comp_meta, by=c("...1" = "IDs")) |>
  relocate(Drug.Name) |>
  dplyr::rename("drug_name" = "Drug.Name", "drug_id" = "...1") |>
  filter(cell_id %in% depmap_lines$depmapId) |>
  filter(!is.na(lfc)) |>
  mutate(drug_response = if_else(lfc < 0, "responder", "non-responder")) |>
  mutate(drug_id_merged = paste0(drug_name, "_", drug_id)) #|>
 # mutate(pdgsc = if_else(drug_name %in% overlap_drugs, "pdgsc", "not_pdgsc"))

# Join with DCRA values from network
final_df <- drug_final |>
  inner_join(depmap_dcra, by=c("cell_id" = "cell_id", "drug_name" = "Drug")) |>
  mutate(concurrency = if_else(drug_response == response, "coherent", "noncoherent"))

### If we only use drugs present in PDGSC screen
pdgsc_data <- read_csv("../GBM_Manuscript_2023/data/GBM_drug_sensitivity_pruned.csv")
overlap_drugs <- intersect(pdgsc_data$Compound, final_df$drug_name)

### xCures drugs
xcures_drugs <- toupper(read_csv("../XCures-SNO/data/xCures_time_on_treatment.csv") |> filter(time_on_treatment !=0) |> pull(therapy) |> unique())
xcures_overlap_drugs <- intersect(xcures_drugs, final_df$drug_name)

### Drugs annotated as Cancer drugs
cancer_drugs <- read_delim("data/cancer_drugs_chembl.csv", delim = ";")
cancer_overlap_drugs <- intersect(cancer_drugs$`Parent Molecule Name`, final_df$drug_name)

### Drugs annotated as Cancer Phase 4 drugs
cancer_phase4_drugs <- read_delim("data/cancer_drugs_chembl.csv", delim = ";") |>
  filter(`Max Phase for Indication` == 4)
cancer_phase4_overlap_drugs <- intersect(cancer_phase4_drugs$`Parent Molecule Name`, final_df$drug_name)

### Distribution of LFC values
p <- ggplot(drug_final, aes(x=drug_id_merged, y=lfc, group=cell_id))
p <- p + geom_boxplot()
p


## PDGSC dataframe
final_df2 <- final_df |>
  filter(drug_name %in% overlap_drugs)

## Cancer drugs dataframe
final_df3 <- final_df |>
  filter(drug_name %in% cancer_overlap_drugs)

## Cancer Phase 4 drugs dataframe
final_df4 <- final_df |>
  filter(drug_name %in% cancer_phase4_overlap_drugs)

## xCures drugs dataframe
final_df5 <- final_df |>
  filter(drug_name %in% xcures_overlap_drugs)



res_summary <- function(df,  drugs_label="All"){
  
  coherence_count <- df |>
    group_by(drug_response) |>
    dplyr::summarise(count=n())
  
  ##### All together
  p1 <- ggplot(df, aes(x=response, y=lfc)) + 
    geom_boxplot(aes(color=response, fill=response), alpha=0.5) +
    geom_text(data=coherence_count, aes(x=drug_response, y=-1, label=count)) +
    scale_color_manual(values=c("#4575b4", "#d73027")) +
    scale_fill_manual(values=c("#4575b4", "#d73027")) +
    geom_jitter(width=0.1, height=0.0000001, aes(color=response)) + 
    geom_hline(yintercept=0, linetype="dotted") +
    xlab("gbmMINER predicted response") + ylab("LFC") +
    ggpubr::stat_compare_means() +
    ggtitle(paste0("PRISM Repurposing screens 23Q2 //  Drugs: ",drugs_label,": ",length(unique(df$drug_name)), " GBM Cell Lines: ", length(unique(df$cell_id)))) +
    theme(axis.text = element_text(size=12)) 
  #p1
  
  p2 <- ggplot(df, aes(x=response, y=lfc)) + 
    geom_violin(aes(color=response, fill=response), alpha=0.5) +
    scale_color_manual(values=c("#4575b4", "#d73027")) +
    scale_fill_manual(values=c("#4575b4", "#d73027")) +
    geom_hline(yintercept=0, linetype="dotted") +
    xlab("gbmMINER predicted response") + ylab("LFC") +
    ggpubr::stat_compare_means() +
    ggtitle(paste0("PRISM Repurposing screens 23Q2 //  Drugs: ",drugs_label,": ",length(unique(df$drug_name)), " GBM Cell Lines: ", length(unique(df$cell_id)))) +
    theme(axis.text = element_text(size=12)) 
  #p2
  
  responder_df <- df |> filter(response == "responder")
  non_responder_df <- df |> filter(response == "non-responder")
  
  
  p3 <- ggplot() + 
    ggdist::stat_halfeye(data=responder_df, aes(x=1.5, y=lfc, color=response, fill=response, side="right"), alpha=0.5, justification=-0.15) +
    ggdist::stat_halfeye(data=non_responder_df, aes(x=1, y=lfc, color=response, fill=response, side="left"), alpha=0.5, justification=1.15) +
    ggpubr::stat_compare_means() +
    geom_violin(data=(df |> filter(response == "responder")), aes(x=1.5, y=lfc, color=response, fill=response), alpha=0.5,width=0.2) +
    geom_violin(data=(df |> filter(response == "non-responder")), aes(x=1.0, y=lfc, color=response, fill=response), alpha=0.5,width=0.2) +
    scale_color_manual(values=c("#4575b4", "#d73027")) +
    scale_fill_manual(values=c("#4575b4", "#d73027")) +
    ggpubr::stat_compare_means() +
    ggtitle(paste0("PRISM Repurposing screens 23Q2 //  Drugs: ",drugs_label,": ",length(unique(df$drug_name)), " GBM Cell Lines: ", length(unique(df$cell_id)))) +
    scale_x_discrete(breaks=c(0,1,2), labels= c("","non-responder","responder")) +
    theme(axis.text = element_text(size=12)) +
    xlab("gbmMINER predicted response") + ylab("LFC")
  #p3
  
  
  coherent_drugs <- df |> 
    group_by(drug_name, concurrency) |>
    dplyr::summarise(count=dplyr::n()) |>
    dplyr::filter(concurrency == "coherent") |>
    dplyr::arrange(desc(count)) |>
    mutate(pdgsc = if_else(drug_name %in% overlap_drugs, "pdgsc_drug", "other_drug"))
  
    
  
  p4 <- ggplot(coherent_drugs, aes(x=factor(drug_name, levels=coherent_drugs$drug_name), y=count, fill=pdgsc))
  p4 <- p4 + geom_bar(stat="identity") +
    scale_fill_manual(breaks=c("pdgsc_drug","other_drug"),values=c("#DD9977","#dddddd"), labels=c("pdgsc_drug","other_drug")) +
    theme(axis.text = element_text(size=12), axis.text.x = element_blank()) +
    labs(x="Drugs", y="Concurrency Count") +
    ggtitle(paste0("PRISM Repurposing screens 23Q2 //  Drugs: ",drugs_label,": ",length(unique(df$drug_name)), " GBM Cell Lines: ", length(unique(df$cell_id)))) 
    
    
  #p4
  grid.arrange(p1,p2,p3,p4, ncol=2, nrow=2)
  
}



plot_dist <- function(df, drugs_label){
  
  responder_df <- df |> filter(response == "responder")
  non_responder_df <- df |> filter(response == "non-responder")
  
  p3 <- ggplot() + 
    ggdist::stat_halfeye(data=responder_df, aes(x=1.5, y=lfc, color=response, fill=response, side="right"), alpha=0.5, justification=-0.15) +
    ggdist::stat_halfeye(data=non_responder_df, aes(x=1, y=lfc, color=response, fill=response, side="left"), alpha=0.5, justification=1.15) +
    ggpubr::stat_compare_means() +
    geom_violin(data=(df |> filter(response == "responder")), aes(x=1.5, y=lfc, color=response, fill=response), alpha=0.5,width=0.2) +
    geom_violin(data=(df |> filter(response == "non-responder")), aes(x=1.0, y=lfc, color=response, fill=response), alpha=0.5,width=0.2) +
    scale_color_manual(values=c("#4575b4", "#d73027")) +
    scale_fill_manual(values=c("#4575b4", "#d73027")) +
    ggpubr::stat_compare_means(data=df, aes(x=response, y=lfc), label.x = 0, label.y = 1) +
    ggtitle(paste0("Drugs: ",drugs_label,": ",length(unique(df$drug_name)), " // GBM Cell Lines: ", length(unique(df$cell_id)))) +
    scale_x_discrete(breaks=c(0,1,2), labels= c("","non-responder","responder")) +
    theme(axis.text = element_text(size=12)) +
    xlab("gbmMINER predicted response") + ylab("LFC")
  p3
}

z1 <- plot_dist(df = final_df, drugs_label = "ALL")
z2 <- plot_dist(df = final_df2, drugs_label = "PDGSC")
z3 <- plot_dist(df = final_df3, drugs_label = "CANCER")
z4 <- plot_dist(df = final_df4, drugs_label = "CANCER/PHASE4")
z5 <- plot_dist(df = final_df5, drugs_label = "xCures")

grid.arrange(z1,z2,z3,z4,z5,z6)




### remoive viloin plot
plot_dist2 <- function(df, drugs_label){
  
  responder_df <- df |> filter(response == "responder")
  non_responder_df <- df |> filter(response == "non-responder")
  
  p3 <- ggplot() + 
    ggdist::stat_halfeye(data=responder_df, aes(x=1.05, y=lfc, color=response, fill=response, side="right"), alpha=0.5) +
    ggdist::stat_halfeye(data=non_responder_df, aes(x=0.95, y=lfc, color=response, fill=response, side="left"), alpha=0.5) +
    #geom_violin(data=(df |> filter(response == "responder")), aes(x=1.5, y=lfc, color=response, fill=response), alpha=0.5,width=0.2) +
    #geom_violin(data=(df |> filter(response == "non-responder")), aes(x=1.0, y=lfc, color=response, fill=response), alpha=0.5,width=0.2) +
    scale_color_manual(values=c("#4575b4", "#d73027")) +
    scale_fill_manual(values=c("#4575b4", "#d73027")) +
    ggpubr::stat_compare_means(data=df, aes(x=response, y=lfc, label = after_stat(p.signif)), label.x = 1, label.y = 1) +
    ggtitle(paste0(drugs_label,": ",length(unique(df$drug_name)))) +
    scale_x_discrete(breaks=c(1,1), labels= c("responder","non-responder")) +
    theme(axis.text = element_text(size=12)) +
    xlab("gbmMINER predicted response") + ylab("LFC")
  p3
}

z1 <- plot_dist2(df = final_df, drugs_label = "ALL")
z2 <- plot_dist2(df = final_df2, drugs_label = "PDGSC")
z3 <- plot_dist2(df = final_df3, drugs_label = "CANCER")
z4 <- plot_dist2(df = final_df4, drugs_label = "CANCER/PHASE4")
z5 <- plot_dist2(df = final_df5, drugs_label = "xCures")

grid.arrange(z6,z3,z4,z2,z5,z7, ncol=3, nrow=2)




## Load PDGSC data
pdgsc <- read_csv("../GBM_Manuscript_2023/data/PDGSC_conccurency_final_data.csv")

responder_df_pdgsc <- pdgsc |> filter(net_class == "responder")
non_responder_df_pdgsc <- pdgsc |> filter(net_class == "nonresponder")

z6 <- ggplot() + 
  scale_y_log10() +
  ggdist::stat_halfeye( data=responder_df_pdgsc, aes(x=1.05, y=hts_act, color=net_class, fill=net_class, side="right"), alpha=0.5) +
  ggdist::stat_halfeye( data=non_responder_df_pdgsc, aes(x=0.95, y=hts_act, color=net_class, fill=net_class, side="left"), alpha=0.5) +
   scale_color_manual(values=c("#4575b4", "#d73027")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  ggpubr::stat_compare_means(data=pdgsc, aes(x=net_class, y=hts_act, label = after_stat(p.signif)), label.x = 1) +
  ggtitle(paste0("PDGSC: ",length(unique(pdgsc$Compound)))) +
  scale_x_discrete(breaks=c(1,1), labels= c("responder","nonresponder")) +
  theme(axis.text = element_text(size=12)) +
  xlab("gbmMINER predicted response") + ylab("log10 IC50")
z6


## Load xCures data
xcures <- read_csv("../GBM_Manuscript_2023/data/xCures_concurrency_final_data.csv") |>
  filter(Drug != "BEVACIZUMAB")

responder_df_xcures <- xcures |> filter(response == "responder")
non_responder_df_xcures <- xcures |> filter(response == "non-responder")

z7 <- ggplot() + 
  #geom_violin(data=xcures, aes(x=response, y=time_on_treatment, fill=response))
  #scale_y_log10() +
  ggdist::stat_halfeye(data=responder_df_xcures, aes(x=1.05, y=time_on_treatment, color=response, fill=response, side="right"), alpha=0.5) +
  ggdist::stat_halfeye(data=non_responder_df_xcures, aes(x=0.95, y=time_on_treatment, color=response, fill=response, side="left"), alpha=0.5) +
  scale_color_manual(values=c("#4575b4", "#d73027")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  ggpubr::stat_compare_means(data=xcures, aes(x=response, y=time_on_treatment, label = after_stat(p.signif)), label.x = 1) +
  ggtitle(paste0("xCures: ",length(unique(xcures$Drug)))) +
  scale_x_discrete(breaks=c(1,1), labels= c("responder","non-responder")) +
  theme(axis.text = element_text(size=12)) +
  xlab("gbmMINER predicted response") + ylab("Time on Treatment (days)")
z7





#Upset plot
library(ComplexHeatmap)
library(viridis)

list_prism <- unique(final_df$drug_name)
list_cancer <- unique(cancer_drugs$`Parent Molecule Name`)
list_cancer_phase4 <- unique(cancer_phase4_drugs$`Parent Molecule Name`)
list_pdgsc <- unique(pdgsc_data$Compound)
list_xcures <- unique(xcures_drugs)

list_all <- list(prism=list_prism, cancer=list_cancer, cancer4=list_cancer_phase4,pdgsc=list_pdgsc, xcures=list_xcures)

upset_list <- make_comb_mat(list_all)
col_size = comb_size(upset_list)
row_size = set_size(upset_list)



ht = UpSet(
  upset_list,
  top_annotation = upset_top_annotation(upset_list, ylim = c(0, max(col_size) * 1.1),show_annotation_name = T),
  right_annotation = upset_right_annotation(upset_list, ylim = c(0, max(row_size) * 1.1)),
  pt_size = unit(5, "mm"), lwd = 3,
  comb_col = c(viridis(7))[comb_degree(upset_list)]
  )

ht = draw(ht)

col_od = column_order(ht)
row_od = row_order(ht)

decorate_annotation("intersection_size", {
  grid.text(col_size[col_od], 
            seq_len(length(col_size)), 
            unit(col_size[col_od], "native") + unit(2, "mm"), 
            default.units = "native", just = "bottom",
            gp = gpar(fontsize = 8))
})

decorate_annotation("set_size", {
  grid.text(row_size[row_od], 
            unit(row_size[row_od], "native") + unit(2, "mm"), 
            rev(seq_len(length(row_size))), 
            default.units = "native", just = "bottom", rot = -90,
            gp = gpar(fontsize = 8))
})










p2 <- ggplot(mat_test, aes(x=net_class, y=hts_act)) + 
  geom_boxplot(aes(color=net_class, fill=net_class), alpha=0.5) +
  scale_y_log10() + scale_color_manual(values=c("#4575b4", "#d73027")) +
  scale_fill_manual(values=c("#4575b4", "#d73027")) +
  geom_jitter(width=0.1, height=0.0000001, aes(color=net_class)) + 
  geom_hline(yintercept=hts_cut1, linetype="dotted") +
  xlab("predicted response") + ylab("log10 IC50") +
  ggtitle(paste(drug, "Wilcox test: ", wil_pval, "| Significant")) +
  theme(axis.text = element_text(size=12)) +
  geom_label(inherit.aes = F, data=mat_test |> group_by(net_class) |> summarise(count=n()), aes(x=net_class, y=1, label=count)) 

p2

write_csv(mat_test, file="data/PDGSC_conccurency_final_data.csv")

## Load PDGSC DCRA
pdgsc_activity <- read_csv(file="data/binary_cutoff_HTSvsSYGNAL_v2.csv")

# heatmap
library(ComplexHeatmap)

pdgsc_mx <- pdgsc_activity |>
  filter(Compound %in% toupper(unique(mat_test$Compound))) |>
  #filter(Compound != "BEVACIZUMAB") |>
  pivot_wider(id_cols = "Compound", names_from = "pid", values_from = "net_act") |>
  column_to_rownames("Compound")


Heatmap(pdgsc_mx,
        col=c("#2c7bb6","yellow", "#ca0020"),
        show_row_dend = F,
        show_column_dend = F,
        row_names_gp = gpar(fontsize=6),
        column_names_gp = gpar(fontsize=6),
        heatmap_legend_param = list(
          title = "Activity")
)




#xCures Boxplot
drug_res_df_final_nobev <- drug_res_df |>
  filter(!is.na(response)) |>
  filter(Drug != "BEVACIZUMAB") |>
  unique()

write_csv(drug_res_df, "data/xCures_concurrency_final_data.csv")

p77 <- ggplot(drug_res_df_final_nobev, aes(x=response, y=time_on_treatment, fill=response))
p77 + geom_boxplot() +
  #geom_point() +
  scale_fill_manual(values=c("#4575b4", "#d73027"),labels=c("nonresponder","responder")) +
  #geom_jitter(width=0.1, height=0.0000001, aes(color=response)) +
  geom_label(inherit.aes = F, data=drug_res_df_final_nobev |> group_by(response) |> summarise(count=n()), aes(x=response, y=0, label=count)) +
  labs(x="Response Prediction", y="Time on Treatment", title=paste0(" No BEV // All Patients // Time on Treatment")) +
  scale_color_manual(values=c("#4575b4", "#d73027")) +
  stat_compare_means(label.y = 475, label.x=1.5,method.args = list(alternative="less")) +
  stat_compare_means(label.y = 450, label.x=1.5, method = "t.test") +
  stat_compare_means(label.y = 425, label.x=1.5, method = "kruskal.test") +
  scale_y_break(c(500,1000)) +
  scale_y_continuous(breaks =c(0,100,200,300,400,500))




### XCure heatmap
xcures_mx <- xcures_activity |>
  filter(Drug %in% toupper(my_drugs)) |>
  filter(Drug != "BEVACIZUMAB") |>
  dplyr::select(Drug,patient,DrugConstrainedRegulonActivity) |>
  unique() |>
  pivot_wider(id_cols = "Drug", names_from = "patient", values_from = "DrugConstrainedRegulonActivity") |>
  column_to_rownames("Drug")


Heatmap(xcures_mx,
        col=c("#2c7bb6","yellow", "#ca0020"),
        show_row_dend = F,
        show_column_dend = F,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=6),
        heatmap_legend_param = list(
          title = "Activity")
)

fs::dir_ls("/Volumes/omics4tb2/SWEDISH_PDGSCs/processed_data/run1_2/",recurse = T,regexp = ".genes.results")

1831+1831+1831+1831+1831+1831+
  
  filing jointly







for(i in 1:21){
  assign("Run_",{{i}}) = i +1
  
}










p2 <- ggplot(final_df, aes(x=response, y=lfc, group=response))
p2 <- p2 +  geom_boxplot() +
  ggpubr::stat_compare_means()
p2




