library(tidyverse)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

d1<- read_delim("~/Downloads/mc3_gene_level_GBM_mc3_gene_level.txt", delim = "\t")

#### Enrichment of prognostic markers in States
# Read state hazard ratios
states_hr <- read_csv("../GBM-Model-052022/State_HR_plot_data.csv") |>
  dplyr::rename("Patient_ID" = "...1") |>
  group_by(label) |>
  mutate(medianGuanScore = median(GuanScore)) # |>
 #mutate(label = label - 1) ## States start with number zero so we correct for that


# boxplot
p <- ggplot(states_hr, aes(factor(label, levels=unique(states_hr$label)), GuanScore, group = label, fill=medianGuanScore))
p <- p + geom_boxplot()
p <- p + scale_fill_gradient2(low = "blue",mid = "white",high = "red", midpoint = 0.5)
p <- p + theme()
p <- p + labs(x = "State", y = "Guan Score")
p

ggsave(p, filename="../GBM-Model-052022/State_Hazard_Boxplot.pdf", device = "pdf",width = 6.61, height= 1.8, units = "in")

## TCGA Gene Mutations
mutated_genes_df <- read_csv("/Volumes/omics4tb2/kkannan/mutations.GBM.TCGA.2019.05.01.csv") |>
  dplyr::rename("Mutation" = "...1") |>
  pivot_longer(!c(Mutation, Gene, GeneID), names_to = "Patient_ID", values_to = "Status")

## Mutation percents (we filter for mutations that are only in 5% of the TCGA cohort patients)
mutated_genes_df_percent <- mutated_genes_df |>
  group_by(Mutation, Status) |>
  summarize(count = n()) |>
  filter(Status == 1) |>
  ungroup()|>
  mutate(percent = round((count*100)/247)) |>
  filter(percent >= 5)

# Prognostic markers selected
gbm_selected <- c("PTEN","EGFR", "TP53", "NF1", "PIK3R1", "PIK3CA","RB1", "ATRX", "IDH1", "PDGFRA", "TERT","H3F3A", "CDKN2A","MUC16","FLG","MYB","MYBL1","MN1")


## Filter the mutations for 5%
mutated_genes_df_filtered <- mutated_genes_df |>
  #filter(Mutation %in% mutated_genes_df_percent$Mutation)
  filter(Gene %in% gbm_selected)

# Combine mutations and states to get states to patient assignments
mutated_genes_states <- left_join(mutated_genes_df_filtered, states_hr, by = "Patient_ID")

# Sumamrise each state
states_vs_mutations <-  mutated_genes_states |>
  group_by(label, Gene) |>
  summarise(Sum = sum(Status)) |>
  filter(!is.na(label))

### Test enrichment of each TCGA mutation in States
res_df <- data_frame()
for(my_mutation in unique(mutated_genes_df_filtered$Mutation)){
  cat("Mutation: ", my_mutation,"\n")
  # get the patients that have the gene mutated
  m0 <- mutated_genes_df |>
    filter(Mutation == my_mutation) |>
    group_by(Mutation) |>
    filter(Status == 1) |>
    pull(Patient_ID)

  m1 <- m0  |>
    unique() |>
    length()

  n = 516 - m1

  for(my_state in unique(states_hr$label)){
    cat("\tState: ", my_state,"\n")
    k0 <- states_hr |>
      filter(label == my_state) |>
      pull(Patient_ID)
    k1 = k0 |>
      unique() |>
      length()

    x = length(intersect(m0,k0))

    res <- phyper(x,m1,n,k1,lower.tail = F)

    tmp1 <- data_frame(State = my_state, Mutation = my_mutation, pvalue = res)
    res_df <- bind_rows(res_df, tmp1)

  }

}

### Filter res_df for Genes of interest
res_df_selected <- res_df |>
  separate(Mutation, into = c("Gene"), sep = "_") |>
  group_by(Gene) |>
  mutate(padjust = p.adjust(pvalue)) |>
  filter(Gene %in% gbm_selected) |>
  filter(padjust <= 0.05) |>
  pivot_wider(id_cols = Gene, names_from = State, values_from = padjust) |>
  column_to_rownames("Gene")

psel <- pheatmap(scale = "none",res_df_selected, cluster_rows = F, cluster_cols = F, fontsize = 8)

psel <- pheatmap(scale = "none",res_df_selected, cluster_rows = F, cluster_cols = F, fontsize = 8,color = colorRampPalette(c("red","blue","gray"))(18) )

ggsave(psel, filename="../GBM-Model-052022/State_vs_GBM_Prognostic_Mutations_Enrichment.pdf", device = "pdf", width=6.61, height = 1.71, units = "in")



## Load Cancer Driver Genes from IntoGen and filter for GBM
gbm_drivers <- read_delim("../GBM-Model-052022/IntoGen_Compendium_Cancer_Genes_02102023.tsv", delim = "\t") |>
  filter(`CANCER_TYPE` == "GBM") |>
  pull(SYMBOL) |>
  unique()

# filter States for states with driver mutations
# order states based on the states boxplot
drivers_states_vs_mutations <- states_vs_mutations |>
  filter(Gene %in% gbm_drivers) |>
  arrange(desc(Sum), .by_group = TRUE) |>
  pivot_wider(id_cols = Gene, names_from = label, values_from = Sum) |>
  column_to_rownames("Gene") |>
  dplyr::select(unique(p$data$label))


pp <- pheatmap(drivers_states_vs_mutations, cluster_rows = F, cluster_cols = F, display_numbers = T)

ggsave(pp, filename="../GBM-Model-052022/State_vs_GBM_Progrnostic_Mutations.pdf", device = "pdf")

# Prognostic markers selected
gbm_selected <- c("PTEN","EGFR", "TP53", "NF1", "PIK3R1", "PIK3CA","RB1", "ATRX", "IDH1", "PDGFRA", "TERT","H3F3A", "CDKN2A","MUC16","FLG","MYB","MYBL1","MN1")

# heatmap for selected drivers
drivers_states_vs_mutations_sel <- states_vs_mutations |>
  filter(Gene %in% gbm_selected) |>
  arrange(desc(Sum), .by_group = TRUE) |>
  pivot_wider(id_cols = Gene, names_from = label, values_from = Sum) |>
  column_to_rownames("Gene") |>
  dplyr::select(unique(p$data$label))


pp_selected <- pheatmap(drivers_states_vs_mutations_sel, cluster_rows = F, cluster_cols = F, display_numbers = T, number_format = "%.0f",border_color = "#ffffff",angle_col = 0,color = brewer.pal(21, "OrRd"))

ggsave(pp_selected, filename="../GBM-Model-052022/State_vs_GBM_Prognostic_Mutations_Selected.pdf", device = "pdf", width=6.61, height = 1.71, units = "in")



# Load identifier mappings
identifiers <- read_delim("../GbmMINER/data/identifier_mappings.txt",delim = "\t") |>
  filter(Source == "Gene Name") |>
  dplyr::select(-Source)

# Read state hazard ratios
states_hr <- read_csv("../GBM-Model-052022/State_HR_plot_data.csv") |>
  dplyr::rename("Patient_ID" = "...1") |>
  group_by(label) |>
  mutate(medianGuanScore = median(GuanScore)) #|>

# Load model files
overExpressedMembersMatrix <- read_csv("../GBM-Model-052022/overExpressedMembers.csv") |>
  column_to_rownames("...1")

underExpressedMembersMatrix <- read_csv("../GBM-Model-052022/underExpressedMembers.csv") |>
  column_to_rownames("...1")
# get the matrix file
dfr = (overExpressedMembersMatrix-underExpressedMembersMatrix)

#read regulons to programs to genes mapping
regulons_programs_genes <- read_csv("../GBM-Model-052022/regulonsProgramsGenesMapping.csv",show_col_types = FALSE)

## Plot gene state activity
gene_state_activity_boxplot <- function(in.gene, in.state){
  gene_ensembl <- identifiers |>
    filter(Name %in% in.gene) |>
    pull(Preferred_Name)
  
  gene_regulons <- regulons_programs_genes |>
    filter(Gene %in% gene_ensembl) |>
    pull(Regulon_ID) |>
    unique()
  
  regulator_regulons <- regulons_programs_genes |>
    filter(Regulator %in% gene_ensembl) |>
    pull(Regulon_ID) |>
    unique()
  
  # Select State Patients
  state_patients <- states_hr |>
    filter(label == in.state) |>
    pull(Patient_ID) |>
    unique()
  
  # Collect expression data from TCGA for selected patients
  state_df <- dfr |>
    #column_to_rownames("Name") |>
    dplyr::select(any_of(state_patients))
  
  # Select gene regulons
  selected_gene_regulons <- state_df[gene_regulons,] |>
    rownames_to_column("Regulon") |>
    pivot_longer(!c(Regulon), names_to = "Patient_ID", values_to = "Activity") |>
    mutate(Type ="Gene Regulon")
  
  # Select gene regulons
  selected_regulator_regulons <- state_df[regulator_regulons,] |>
    rownames_to_column("Regulon") |>
    pivot_longer(!c(Regulon), names_to = "Patient_ID", values_to = "Activity") |>
    mutate(Type ="Regulator Regulon")
  
  selected_regulons <- bind_rows(selected_gene_regulons, selected_regulator_regulons) |>
    mutate(Gene = in.gene, State = in.state)
  
  return(selected_regulons)
}



in.genes <- gbm_selected
in.state <- "0"

gene_activity_vs_states <- tibble()
for(in.gene in in.genes){
  for(in.state in unique(states_hr$label)){
    cat("Analyzing Gene: ", in.gene, "\n", "\tState: ", in.state, "\n")
    res1 <- gene_state_activity_boxplot(in.gene = in.gene, in.state = in.state)
    gene_activity_vs_states <- bind_rows(gene_activity_vs_states, res1)
  }
}

gene_activity_vs_states_sum <- gene_activity_vs_states |>
  group_by(State, Gene) |>
  mutate(meanActivity = mean(Activity)) |>
  ungroup() |>
  dplyr::select(State, Gene, meanActivity) |>
  unique() |>
  pivot_wider(id_cols = Gene, names_from = State, values_from = meanActivity) |>
  column_to_rownames("Gene")

pheatmap(gene_activity_vs_states_sum, cluster_cols = F)


ggplot(gene_activity_vs_states |> filter(Gene == "MUC16" & Type == 'Gene Regulon'), aes(x=factor(State, levels = unique(states_hr$label)), y=Activity, group=State)) +
  geom_violin()




####### Subtype enrichment
tcga_survival <- read_csv("/Volumes/omics4tb2/kkannan/TCGA_Survival_Gbm.csv")


states_survival <- states_hr |>
  left_join(tcga_survival, by="Patient_ID")

subtypes <- states_survival|>
  group_by(Subtype) |>
  summarise(subtype_total=n())


### Test enrichment of each TCGA mutation in States
sub_df <- data_frame()
for(my_subtype in unique(states_survival$Subtype)){
  cat("Subtype: ", my_subtype,"\n")
  # get the patients that have the gene mutated
  m0 <- states_survival |>
    filter(Subtype == my_subtype) |>
    pull(Patient_ID) |>
    unique()

  m1 <- length(m0)

  n = 516 - m1

  for(my_state in unique(states_hr$label)){
    cat("\tState: ", my_state,"\n")
    k0 <- states_hr |>
      filter(label == my_state) |>
      pull(Patient_ID)
    k1 = k0 |>
      unique() |>
      length()

    x = length(intersect(m0,k0))

    res <- phyper(x,m1,n,k1,lower.tail = F)

    tmp1 <- data_frame(State = my_state, Subtype = my_subtype, pvalue = res)
    sub_df <- bind_rows(sub_df, tmp1)

  }

}

## Sbtype enrichment
sub_plot_df <- sub_df |>
  mutate(padjust = p.adjust(pvalue)) |>
  mutate(minusLogpvalue = -log10(padjust))

p <- ggplot(sub_plot_df, aes(x = factor(State, levels = unique(states_hr$label)), y=minusLogpvalue, group=State, fill=factor(Subtype, levels = rev(c("proneural","neural","g_cimp","classical","mesenchymal"))), alpha=if_else(minusLogpvalue > 1.30103, 1, 0.5))) +
  geom_bar(stat="identity") +
  labs(x="States", y="Count",) +
  theme(legend.text = element_text(size=7)) +
  geom_hline(yintercept = 1.30103, lty="dashed") +
  #scale_fill_brewer(type = "seq",palette = "YlGnBu")
  scale_fill_manual(name="Subtype",
                    values = c(
                      "proneural" = "#ffffbf",
                      "neural" = "#abdda4",
                      "g_cimp" = "#2b83ba",
                      "classical" = "#fdae61",
                      "mesenchymal" = "#d7191c"
                    )
  )


ggsave(p, filename="../GBM-Model-052022/State_vs_Subtypes.pdf", device = "pdf", width=6.61, height = 1.71, units = "in")



subtype_count <- states_survival |>
  select(label,Subtype) |>
  group_by(label,Subtype) |>
  summarise(n=n()) |>
  ungroup() |>
  left_join(subtypes, by="Subtype") |>
  mutate(subtype_norm = n/subtype_total*100) |>
  arrange(match(Subtype, factor(Subtype, levels = rev(c("proneural","neural","g_cimp","classical","mesenchymal")))))


p <- ggplot(subtype_count, aes(x = factor(label, levels = unique(states_hr$label)), y=n, group=label, fill=factor(Subtype, levels = rev(c("proneural","neural","g_cimp","classical","mesenchymal"))))) +
  geom_bar(stat="identity") +
  labs(x="States", y="Count",) +
  theme(legend.text = element_text(size=7)) +
  #scale_fill_brewer(type = "seq",palette = "YlGnBu")
  scale_fill_manual(name="Subtype",
                    values = c(
    "proneural" = "#ffffbf",
    "neural" = "#abdda4",
    "g_cimp" = "#2b83ba",
    "classical" = "#fdae61",
    "mesenchymal" = "#d7191c"
      )
    )



# ggplot(subtype_count |> filter(Subtype %in% c("classical","mesenchymal","proneural")), aes(x = factor(label, levels = unique(states_hr$label)), y=n, group=label, fill=Subtype)) +
#   geom_bar(stat="identity")

ggplot(states_survival, aes(x = factor(label, levels = unique(states_hr$label)), y=GuanScore.x, group=label)) +
  geom_boxplot()
