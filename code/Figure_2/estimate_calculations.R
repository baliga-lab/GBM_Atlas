# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2024
# Email: sturkarslan@systemsbiology.org
#
# Date: 2024-05-17
#
# Script Name: estimate_calculations.R
#
# Script Description: Uses ESTIMATE algorithm to calculate ImmuneScore for States and Programs
#
# Script version:
# --------------------------------------------

library(tidyverse)
library(utils)
library(estimate)
library(CePa)
library(xCell)
library(Rtsne)

# Load identifier mappings
identifiers <- read_delim("../GbmMINER/data/identifier_mappings.txt",delim = "\t")# |>

# Read state hazard ratios
states_hr <- read_csv("../GBM-Model-052022/State_HR_plot_data.csv") |>
  dplyr::rename("Patient_ID" = "...1") |>
  group_by(label) |>
  mutate(medianGuanScore = median(GuanScore))


# Load cohort expression data
exp_data <- read_csv("/Volumes/omics4tb2/SYGNAL/GBM-Serdar/data/GbmMicroRNAMergedWithIDsZScored.csv") |>
  dplyr::rename("Gene" = "...1") |>
  left_join(identifiers, by=c("Gene" = "Preferred_Name")) |>
  filter(Source == "Gene Name") |>
  relocate(Name) |>
  dplyr::select(-Gene, -Source)

exp_data2 <- exp_data |>
  column_to_rownames("Name")

my_colnames <- gsub(".1", "", colnames(exp_data)[-1], fixed = T)
colnames(exp_data2) <- my_colnames
  
write.csv(exp_data2, file="data/GBM_exp_zscored_formatted.csv", row.names = T)

## Combine with common genes if neded
# merged_df <- left_join(common_genes, exp_data, by=c("GeneSymbol" = "Name")) |>
#   column_to_rownames("GeneSymbol") |>
#   dplyr::select(-EntrezID, -Synonyms,-GeneName, -Chromosome)

# Convert to matrix
merged_df <- exp_data |>
  column_to_rownames("Name")

###### Functiion to calculate ImmuneScore, StromalScore and EstimateScore by Using ESTIMATE package
# https://bioinformatics.mdanderson.org/public-software/estimate/
calculate_state_estimate <- function(in.state){
  # Select State Patients
  state_patients <- states_hr |>
    filter(label == in.state) |>
    pull(Patient_ID) |>
    unique()

  # Collect expression data from TCGA for selected patients
  state_df <- merged_df |>
    dplyr::select(any_of(state_patients))

  # write output to gct file
  outputGCT(state_df, "output/state.exp.gct")

  # Run estimte analysis
  res <- estimateScore("output/state.exp.gct",output.ds = "output/state_out")

  # Read output
  tmp_out <- read.gct("output/state_out") |>
    as.data.frame() |>
    rownames_to_column("Type") |>
    pivot_longer(cols = starts_with("TCGA")) |>
    mutate(State = in.state)

  # return results
  return(tmp_out)

}

## Process all the states with ESTIMATE calculations
state_estimate_df <- tibble()
for(in.state in unique(states_hr$label)){
  cat("Analyzing State: ", in.state, "\n" )
  # Run analysis
  est_out <- calculate_state_estimate(in.state)
  # Collect results for all
  state_estimate_df <- bind_rows(state_estimate_df, est_out )
}


## Statistical correlation of meanGuanScores and meanImmuneScore for states
# modify risk scores df
guan_df <- states_hr |>
  select(Patient_ID, GuanScore, label) |>
  dplyr::rename("State" = label) |>
  unique()

# modify immunescores df
estimate_df <- state_estimate_df |>
  filter(Type == "ImmuneScore") |>
  dplyr::rename("Patient_ID" = name, "ImmuneScore" = value) |>
  select(Patient_ID, ImmuneScore, State) |>
  unique() |>
  mutate(Patient_ID = gsub(".", "-", Patient_ID, fixed = T))

# join risk scores and immunescores df
immune_df <- guan_df |>
  full_join(estimate_df, by=c("Patient_ID", "State")) |>
  group_by(State) |>
  mutate(meanGuanScore = mean(GuanScore)) |>
  mutate(meanImmuneScore = mean(ImmuneScore)) |>
  select(State, meanGuanScore, meanImmuneScore) |>
  unique() |>
  mutate(GuanScore_Category = ifelse(meanGuanScore >= 0.5, "high risk", "low risk"),
         ImmuneScore_Category = ifelse(meanImmuneScore > 0, "high risk", "low risk"))


# Assuming you have already loaded the required libraries and have the data in a tibble called "your_tibble".
# If not, please adjust the data frame name accordingly.

# Step 1: Recode GuanScore and ImmuneScore to "high risk" categories

# Step 2: Cross-tabulate the two categorical variables
cross_tab <- table(immune_df$GuanScore_Category, immune_df$ImmuneScore_Category)

# Step 3: Perform the Chi-squared test
chi_squared_test <- chisq.test(cross_tab)

# Step 4: Print the results
print(chi_squared_test)


#Pearson's Chi-squared test with Yates' continuity correction

#data:  cross_tab
#X-squared = 5.2356, df = 1, p-value = 0.02213


# states_risks <- states_hr |>
#   select(label, medianGuanScore) |>
#   dplyr::rename("State" = label) |>
#   unique()
# 
# states_estimates_filt <- state_estimate_df |>
#   filter(Type == "ImmuneScore") |>
#   select(State, value) |>
#   unique() |>
#   group_by(State) |>
#   mutate(medianValue = median(value)) |>
#   select(State, medianValue) |>
#   unique()

# # Boxplot of Estimate scores for ImmuneScore
# # Create plot dataframe
# immune_df <- states_estimates_filt |>
#   left_join(states_risks, by = c("State" = "label")) |>
#   mutate(color = if_else(medianValue > 0, "red","gray"))


# Assuming your data frame is named "df"
# Replace "df" with the actual name of your data frame if it's different.

# Step 1: Load required libraries (if not already installed)
# install.packages("dplyr")   # Uncomment and run this line if you haven't installed the "dplyr" package
# library(dplyr)              # Uncomment and run this line to load the library

# Step 2: Calculate the correlation between ImmuneScore and RiskScore by State
cor_by_state <- immune_df %>%
  group_by(State) %>%
  summarize(correlation = cor(medianGuanScore, medianImmuneScore))

# Step 3: Print the correlation results
print(cor_by_state)


# Assuming you have already calculated the correlation by state and stored it in "cor_by_state".
# If not, please run the previous code to calculate the correlation first.

# Step 1: Extract the correlation values from the "cor_by_state" data frame
cor_values <- cor_by_state$correlation

# Step 2: Use Fisher's Z-transformation to calculate the overall correlation and its standard error
overall_cor <- sum(cor_values) / length(cor_values)
z_transformation <- sqrt(length(cor_values) - 3) * overall_cor / sqrt(1 - overall_cor^2)

# Step 3: Calculate the p-value using the t-distribution
p_value <- 2 * pt(abs(z_transformation), df = length(cor_values) - 2)

# Step 4: Print the overall correlation and p-value
print(paste("Overall Correlation:", overall_cor))
print(paste("P-value:", p_value))




# Boxplot of Estimate scores for ImmuneScore
# Create plot dataframe
immune_df <- state_estimate_df |>
  left_join(states_hr, by = c("State" = "label")) |>
  filter(Type == "ImmuneScore") |>
  group_by(State) |>
  mutate(medianValue = median(value)) |>
  mutate(color = if_else(medianValue > 0, "red","gray"))

# Plot the boxplot
p <- ggplot(immune_df, aes(x = factor(State, levels = unique(states_hr$label)), y=value, group = State, fill=color)) +
  geom_boxplot() +
  scale_fill_manual(values = c("gray","red")) +
  geom_hline(yintercept = 0,lty="dashed" ) +
  labs(x="States", y="ImmuneScore") +
  theme(legend.position = "none")
  #scale_fill_gradient2(low = "blue",mid = "white",high = "red", midpoint = 0) +
# Write to a file
ggsave(p, filename="../GBM-Model-052022/State_ImmuneScore_Boxplot.pdf", device = "pdf",width = 6.61, height= 1.8, units = "in")



###### Functiion to calculate Cell Types by usog XCell ackage
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1#Sec1
calculate_state_celltypes <- function(in.state){
    # Select State Patients
    state_patients <- states_hr |>
      filter(label == in.state) |>
      pull(Patient_ID) |>
      unique()
    # Collect expression data from TCGA for selected patients
    state_df <- exp_data |>
      column_to_rownames("Name") |>
      dplyr::select(any_of(state_patients))
    # Run analysis
    res_xcell <- xCellAnalysis(state_df)
    # collect results
    res_xcell <- as.data.frame(res_xcell) |>
      rownames_to_column("Cell") |>
      pivot_longer(cols = starts_with("TCGA")) |>
      mutate(State = in.state)

    return(res_xcell)
}

## Run cell type analysis for all states
state_celltypes_df <- tibble()
for(in.state in unique(states_hr$label)){
  cat("Analyzing State: ", in.state, "\n" )
  est_out <- calculate_state_celltypes(in.state)
  state_celltypes_df <- bind_rows(state_celltypes_df, est_out )
}

# You can run the analysis for selected cell types
selected_cells <- c("Smooth muscle", "Basophils", "iDC", "MSC","HSC","Th1 cells","NKT","cDC","Th2 cells","Endothelial cells","MicroenvironmentScore","ImmuneScore","StromaScore")

# Boxplot of the results
ggplot(state_celltypes_df |> filter(Cell == "Th2 cells"), aes(x = factor(State, levels = unique(states_hr$label)), y=value, group=State, fill=as.character(State))) +
  geom_boxplot() +
  geom_point()

## Dataframe for creating a heatmap
state_celltypes_summary <- state_celltypes_df |>
  group_by(State, Cell) |>
  dplyr::summarise(mean = mean(value)) |>
  arrange(match(State, states_hr$label)) |>
  mutate(mean = as.numeric(mean)) |>
  filter(Cell %in% selected_cells) |>
  #filter(mean >=0.2) |>
  pivot_wider(id_cols = Cell, names_from = State, values_from = mean) |>
  column_to_rownames("Cell")

# Plot state cell types enrichment
ht1 <- Heatmap(state_celltypes_summary,
               cluster_rows = T,
               cluster_columns = F,
               show_row_dend = F,
               row_names_gp = gpar(fontsize=6)

                )
ht1

pdf(file="../GBM-Model-052022/State_CellType_Enrichment.pdf", width=6.61, height = 1)
ht1
dev.off()


program_celltypes_summary <- program_celltypes_df |>
  group_by(Program, Cell) |>
  dplyr::summarise(mean = mean(value)) |>
  arrange(match(Program, disease_programs$Program_ID)) |>
  mutate(mean = as.numeric(mean)) |>
  #filter(mean >=0.5) |>
  pivot_wider(id_cols = Cell, names_from = Program, values_from = mean) |>
  column_to_rownames("Cell")


#### Calculate enrichmen of cell types in patients where the program is active
## Load Program activity for cohort
program_activity <- read_csv("../GBM-Model-052022/cohortProgramActivity.csv") |>
  dplyr::rename("Program" = "...1")

patients_overactive_programs <- program_activity |>
  pivot_longer(cols=starts_with("TCGA"),names_to = "Patients") |>
  dplyr::filter(value == 1)

disease_programs <- read_csv("../GBM-Model-052022/CoxProportionalHazardsPrograms.csv") |>
  filter(`p-value` <= 0.05) |>
  dplyr::rename("Program_ID" = "...1") |>
  mutate(Program_ID = as.character(Program_ID)) |>
  arrange(desc(HR))


calculate_program_celltypes <- function(in.program){

  # Select State Patients
  program_patients <- patients_overactive_programs |>
    filter(Program == in.program) |>
    pull(Patients) |>
    unique()

  # Collect expression data from TCGA for selected patients
  programs_df <- exp_data |>
    column_to_rownames("Name") |>
    dplyr::select(any_of(program_patients))

  res_xcell <- xCellAnalysis(programs_df)

  res_xcell <- as.data.frame(res_xcell) |>
    rownames_to_column("Cell") |>
    pivot_longer(cols = starts_with("TCGA")) |>
    mutate(Program = in.program)

  return(res_xcell)

}



program_celltypes_df <- tibble()
for(in.program in unique(disease_programs$Program_ID)){
  cat("Analyzing Program: ", in.state, "\n" )

  est_out <- calculate_program_celltypes(in.program)

  program_celltypes_df <- bind_rows(program_celltypes_df, est_out )
}


for(in.cell in unique(program_celltypes_df$Cell)){
  cat("Making plot for ", in.cell,"\n" )
  # Boxplot
  # Create plot dataframe
  program_cells <- program_celltypes_df |>
    left_join(disease_programs, by = c("Program" = "Program_ID")) |>
    filter(Cell == in.cell) |>
    group_by(Program) |>
    mutate(medianValue = median(value)) |>
    mutate(color = if_else(medianValue > 0.5, "red","gray")) |>
    mutate(risk = if_else(HR > 0, "high", "low"))

  # Plot the boxplot
  p <- ggplot(program_cells, aes(x = factor(Program, levels = unique(disease_programs$Program_ID)), y=value, group = Program, fill=color)) +
    geom_boxplot() +
    scale_fill_manual(values = c("gray","red")) +
    geom_hline(yintercept = 0,lty="dashed" ) +
    labs(title = paste0("Cell Type: ", in.cell), x="Programs", y=paste0(in.cell)) +
    theme(legend.position = "none") +
    facet_grid(.~risk,scales = "free_x", space = "free") +
    theme_classic() +
    theme(legend.position = "none", strip.text = element_blank(), panel.spacing.y = unit(0, "pt"), axis.text.x = element_text(angle=90))

  ggsave(p, filename=paste0("../GBM-Model-052022/Cell_Type_Enrichments/Cell_", in.cell, "_enrichment.pdf"), device = "pdf",width = 6.61, height= 3, units = "in",dpi = "retina")


}



### Example Program Immune Score caluclation
in.program = "172"
#### Calculate Program Immune scores
calculate_program_estimates <- function(in.program){
  # Select State Patients
  program_patients <- patients_overactive_programs |>
    filter(Program == in.program) |>
    pull(Patients) |>
    unique()
  # Collect expression data from TCGA for selected patients
  programs_df <- exp_data |>
    column_to_rownames("Name") |>
    dplyr::select(any_of(program_patients))
  # write output to gct file
  outputGCT(programs_df, "output/program.exp.gct")
  # Run estimte analysis
  res <- estimateScore("output/program.exp.gct",output.ds = "output/program_out")
  # Read output
  tmp_out <- read.gct("output/program_out") |>
    as.data.frame() |>
    rownames_to_column("Type") |>
    pivot_longer(cols = starts_with("TCGA")) |>
    mutate(Program = in.program)
  # return results
  return(tmp_out)
}

program_estimates_df <- tibble()
for(in.program in unique(disease_programs$Program_ID)){
  cat("Analyzing Program: ", in.state, "\n" )

  est_out <- calculate_program_estimates(in.program)

  program_estimates_df <- bind_rows(program_estimates_df, est_out )
}

# Boxplot of Estimate scores for ImmuneScore
# Create plot dataframe
immune_df_program <- program_estimates_df |>
  left_join(disease_programs, by = c("Program" = "Program_ID")) |>
  filter(Type == "ImmuneScore") |>
  group_by(Program) |>
  mutate(medianValue = median(value)) |>
  mutate(color = if_else(medianValue > 0, "red","gray")) |>
  mutate(risk = if_else(HR > 0, "high", "low"))

# Plot the boxplot
p <- ggplot(immune_df_program, aes(x = factor(Program, levels = unique(disease_programs$Program_ID)), y=value, group = Program, fill=color)) +
  geom_boxplot() +
  scale_fill_manual(values = c("gray","red")) +
  geom_hline(yintercept = 0,lty="dashed" ) +
  labs(x="Programs", y="ImmuneScore") +
  theme(legend.position = "none") +
  facet_grid(.~risk,scales = "free_x", space = "free") +
  theme_classic() +
  theme(legend.position = "none", strip.text = element_blank(), panel.spacing.y = unit(0, "pt"))
#scale_fill_gradient2(low = "blue",mid = "white",high = "red", midpoint = 0) +
# Write to a file
ggsave(p, filename="../GBM-Model-052022/Program_ImmuneScore_Boxplot.pdf", device = "pdf",width = 6.61, height= 1.2, units = "in",dpi = "retina")