# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2023
# Email: sturkarslan@systemsbiology.org
#
# Date: 2023-03-08
#
# Script Name: program_functional_enrichment
#
# Script Description: Given the program name calculates functional
# enrichment across several msigDB categories
#
# Script version: v1.0.0
# --------------------------------------------
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ReactomePA)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(GOSemSim)

# Build Annotation data
hsGO <- godata('org.Hs.eg.db', ont="BP", keytype = "ENSEMBL")

##### Function to find functional enrichment
program_functional_enrichment <- function(in_program, disease.relevant=FALSE){

  #read regulons to programs to genes mapping
  regulons_programs_genes <- read_csv("../GBM-Model-052022/regulonsProgramsGenesMapping.csv",show_col_types = FALSE)

  # Read cancer hallmarks GO mapping
  hallmarks_list <- read_rds("../GBM-Model-052022/hallmarks_to_GOTerms.rds",)

  ## Load disease relevant regulons
  disease_regulons <- read_csv("../GBM-Model-052022/diseaseRelevantRegulonsSummary.csv",show_col_types = FALSE)

  # Load MSIGDB Hallmarks set
  h_df = msigdbr(category = "H")
  h_t2g = h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()

  # Load MSIGDB Curated dataset
  c2_df = msigdbr(category = "C2")
  c2_t2g = c2_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()


  if(disease.relevant == TRUE){
    # Get genes for the programs
    my_genes <- regulons_programs_genes |>
      filter(Programs == in_program) |>
      filter(Regulon_ID %in% disease_regulons$Regulon_ID) |>
      pull(Gene) |>
      unique()
  } else {
    # Get genes for the programs
    my_genes <- regulons_programs_genes |>
      filter(Programs == in_program) |>
      pull(Gene) |>
      unique()
  }


  # Convert to entrezID
  my_entrez <-
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = my_genes,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  # Convert to SYMBOL
  my_symbols <-
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = my_genes,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

  # Reactome Enrichment
  cat("Running Reactome Enrichment...\n")
  reactome_enrich <-
    enrichPathway(
      gene = my_entrez$ENTREZID,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      readable = T
    )

  # Hallmarks enrichment
  cat("Running MsigDB Hallmarks Enrichment...\n")
  hallmark_enrich <-
    enricher(
      gene = my_symbols$SYMBOL,
      TERM2GENE = h_t2g,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )

  # Curated Geneset enrichment
  cat("Running MsigDB Curated Geneset Enrichment...\n")
  curated_enrich <-
    enricher(
      gene = my_symbols$SYMBOL,
      TERM2GENE = c2_t2g,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )

  # Plot enrichments
  reactome_plot <- dotplot(reactome_enrich, font.size=8, title="Reactome Enrichment")
  hallmark_plot <- dotplot(hallmark_enrich, font.size=8, title="Hallmarks Enrichment")
  curated_plot <- dotplot(curated_enrich, font.size=8, title="Curated Enrichment")

  ## Cancer Hallmarks
  # run GO BP enrichment
  cat("Running GO Enrichment...\n")
  ego2 <- enrichGO(gene         = my_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

  # get only GO terms
  program_go_terms <- unique(ego2@result$ID)

  # Loop through each hallmark
  cat("Running Cancer Hallmarks Enrichment...\n")
  program_sem_res <- tibble()
  for(myhallmark in names(hallmarks_list)){
    cat("\t\t ", myhallmark, "\n")

    # Run semantic similarity
    gosem_res1 = mgoSim(program_go_terms, hallmarks_list[[myhallmark]], semData = hsGO, measure = "Lin", combine = "max")
    gosem_res2 = mgoSim(program_go_terms, hallmarks_list[[myhallmark]], semData = hsGO, measure = "Wang", combine = "max")
    gosem_res3 = mgoSim(program_go_terms, hallmarks_list[[myhallmark]], semData = hsGO, measure = "Jiang", combine = "max")
    # Build results tibble
    tmp1 <- tibble(program = in_program, hallmark = myhallmark, linScore = gosem_res1, wangScore = gosem_res2, jiangScore = gosem_res3)
    # Add to master table
    program_sem_res <- bind_rows(program_sem_res, tmp1)
  }

  # Returm results
  return(
    list(
      reactome_enrichment = reactome_enrich@result,
      reactome_plot = reactome_plot,
      hallmark_enrichment = hallmark_enrich@result,
      hallmark_plot = hallmark_plot,
      curated_enrichment = curated_enrich@result,
      curated_plot = curated_plot,
      cancer_hallmarks = program_sem_res
    )
  )


}
