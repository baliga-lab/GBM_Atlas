# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2024
# Email: sturkarslan@systemsbiology.org
#
# Date: 2024-04-22
#
# Script Name: filter_CM_flows
#
# Script Description: Given the CM flows and disease relevant regulons, filter disease relevant CM Flows
#
# Script version: v1
# --------------------------------------------

### Calculate Disease Relevant CM Flows
library(tidyverse)

# Load identifier mappings
id_mappings <- read_delim("../GbmMINER/data/identifier_mappings.txt",delim = "\t") |>
  filter(Source == "Gene Name") |>
  dplyr::select(-Source)

# disease relevant regulons
disease_regulons <- read_csv("../GBM-Model-052022/diseaseRelevantRegulonsSummary.csv")


# 1. Load CM Flows
####################### CM Flows  #######################
# Load causal flows and format
cm_flows_all <- read_csv("../GBM-Model-052022/filteredCausalResults.csv") |>
  dplyr::select(-`...1`) |>
  separate(Mutation, into = c("MutationSymbol", "code", "potential", "somatic"), sep = "_", remove = F) |>
  mutate(Regulon = gsub("R-", "", Regulon)) |>
  #filter(Regulon %in% disease_regulons$Regulon_ID) |>
  left_join(id_mappings, by=c("MutationSymbol" = "Name")) |>
  dplyr::rename("MutationEnsembl" = "Preferred_Name") |>
  left_join(id_mappings, by=c("Regulator" = "Preferred_Name")) |>
  dplyr::rename("RegulatorSymbol" = "Name") |>
  dplyr::rename("RegulatorEnsembl" = "Regulator") |>
  dplyr::select(-code,-potential,-somatic) |>
  dplyr::relocate(Mutation,MutationSymbol,RegulatorSymbol,Regulon,MutationEnsembl,RegulatorEnsembl) |>
  mutate(CMF_ID = paste0("CMF-", str_pad(row_number(), 5, pad = "0")))

# Write to file
write_csv(cm_flows_all, file="../GBM-Model-052022/AllCMFlowsAnnotated.csv")

# Filter for disease relevant regulons
cm_flows_disease <- cm_flows_all |>
  filter(Regulon %in% disease_regulons$Regulon_ID)

# Get CM flows genes
cm_flow_genes_disease <- cm_flows_disease|>
  filter(!grepl("_PID_", Mutation)) |>
  filter(!grepl("SIGNALING|DISSASSEMBLY", Mutation)) |>
  mutate(CMFlow = paste0(MutationSymbol,
                         "_",
                         if_else(MutationRegulatorEdge == 1, "up", "dn"),
                         "_",
                         RegulatorSymbol,
                         "_",
                         if_else(RegulatorRegulon_Spearman_R > 0, "up", "dn"),
                         "_",
                         Regulon
  )
  ) |>
  mutate(CMFlowType = "Mutated Gene") |>
  relocate(CMF_ID, CMFlow)

# Write to file
write_csv(cm_flow_genes_disease, file="diseaseRelevantCMFlowsGenes.csv")


# Get CM flows genes
cm_flow_pathways_disease <- cm_flows_disease |>
  filter(grepl("_PID_", Mutation)) |>
  filter(grepl("SIGNALING|DISSASSEMBLY", Mutation)) |>
  mutate(MutationSymbol = NA) |>
  separate(Mutation, into =c("Pathway"), sep= "_PID_NCI_NATURE_", remove = F) |> # Format mutation to get pathway
  mutate(CMFlow = paste0(Pathway, # create a concatenated CM flow summary
                         "_",
                         if_else(MutationRegulatorEdge == 1, "up", "dn"),
                         "_",
                         RegulatorSymbol,
                         "_",
                         if_else(RegulatorRegulon_Spearman_R > 0, "up", "dn"),
                         "_",
                         Regulon
  )
  ) |>
  mutate(CMFlowType = "Mutated Pathway") |>
  relocate(CMF_ID, CMFlow)

# Write to file
write_csv(cm_flow_pathways_disease, file="diseaseRelevantCMFlowsPathways.csv")

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

## Filter disease relevant CM flows for 5 percent
cm_flow_genes_disease |>
  filter(Mutation %in% mutation_percents_genes$Mutation) |>
  pull(MutationSymbol) |>
  unique()
## # of smatically mutated disease relevant genes: 30


## TCGA Pathway Mutations
mutated_pathways_df <- read_csv("/Volumes/omics4tb2/kkannan/NonSilentPathwayMutationsGbmRNAAndMicroarray07242020.csv") |>
  dplyr::rename("Mutation" = "...1") |>
  dplyr::select(where(~ any(. != 0))) |>
  #dplyr::select(-Gene, -GeneID ) |>
  pivot_longer(!c(Mutation), names_to = "Patient_ID", values_to = "Status")

## Mutation percents pathways
mutation_percents_pathways <- mutated_pathways_df |>
  group_by(Mutation, Status) |>
  summarize(count = n()) |>
  filter(Status == 1) |>
  ungroup()|>
  mutate(percent = round((count*100)/174)) |>
  filter(percent >= 5)

cm_flow_pathways_disease |>
  filter(Mutation %in% mutation_percents_pathways$Mutation) |>
  pull(Mutation) |>
  unique()

## of somatically mutated pathways 49