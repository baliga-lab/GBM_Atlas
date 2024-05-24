# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2023
# Email: sturkarslan@systemsbiology.org
#
# Date: 2023-03-08
#
# Script Name: plot_program_survival
#
# Script Description: Plots program survival KM Curves
#
# Script version: v1.0.0
# --------------------------------------------
library(tidyverse)
library(ggsurvfit)
library(survival)
library(survRM2)
library(gghighlight)
library(cowplot)
# load custom code
source("code/load_program_survival.R")

## Load survival info for all programs
survival_df_programs <- load_program_survival()

######## DATA FILES #########
## Load survival data
tcga_survival <- read_csv("../GBM-Model-052022/guanSurvivalDf_TCGA_GBM.csv",show_col_types = FALSE)
## Load Program activity for cohort
program_activity <- read_csv("../GBM-Model-052022/cohortProgramActivity.csv",show_col_types = FALSE) |>
  dplyr::rename("Program" = "...1")

###############################################
##### Function to build survival plot #####
plot_program_survival <- function(in_program){
  # Select relevant program activity
  selected_program_activity <- program_activity |>
    dplyr::filter(Program == in_program) |>
    column_to_rownames("Program")

  ## Get Over active patients for the given program
  program_over_patients <- colnames(selected_program_activity[,selected_program_activity == 1])

  ## Get Neutral active patients for the given program
  program_neutral_patients <- colnames(selected_program_activity[,selected_program_activity == 0])

  ## Get Under active patients for the given program
  program_under_patients <- colnames(selected_program_activity[,selected_program_activity == -1])

  overactive_survival <- tcga_survival |>
    dplyr::filter(Patient_ID %in% program_over_patients ) |>
    mutate(Activity = "Over")

  underactive_survival <- tcga_survival |>
    dplyr::filter(Patient_ID %in% program_under_patients ) |>
    mutate(Activity = "Under")

  neural_survival <- tcga_survival |>
    dplyr::filter(Patient_ID %in% program_neutral_patients ) |>
    mutate(Activity = "Neutral")

  # Combine all activities
  combined_survival <- bind_rows(overactive_survival, underactive_survival, neural_survival)

  # Survival plot
  program_plot <- survfit2(Surv(duration, observed) ~ as.character(Activity), data = combined_survival) |>
    ggsurvfit(theme = theme_light()) +
    annotate("text", x = 2500, y = 0.85, label=paste0("Program: ",in_program)) +
    labs(title = "Survival") +
    scale_color_manual(
      name = "Activity",
      breaks = c("Under", "Neutral", "Over"),
      values = c("Under" = "blue", "Neutral" = "gray", "Over" = "red"))

  print(program_plot)
  return(program_plot)
}

###############################################
# Function to plot survival plot for a program
program_survival <- function(in_program){

  ## Filter for min and max program KM data
  survival_df2_programs <- survival_df_programs |>
    filter(program %in% in_program)

  ## Create tge KM plot for all programs faded except min and max
  program_plot <- survfit2(Surv(duration, observed) ~ as.character(program), data = survival_df2_programs) |>
    ggsurvfit(theme = theme_light()) +
    labs(x="Time (days)") +
    annotate("text", x = 2500, y = 0.85, label=paste0("Program: ",in_program)) +
    scale_color_brewer(type="div", palette = "RdGy") +
    theme(legend.title=element_text(size=10),
          legend.text=element_text(size=10),
          legend.position = c(.50, .65))

  print(program_plot)
}
