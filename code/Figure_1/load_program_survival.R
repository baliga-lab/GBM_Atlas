# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2023
# Email: sturkarslan@systemsbiology.org
#
# Date: 2023-03-08
#
# Script Name: load_program_survival
#
# Script Description: Loads all the survival info for the programs
#
# Script version: v1.0.0
# --------------------------------------------
library(tidyverse)

##### Load Program survival data #####
load_program_survival <- function(){
  # load survival data
  survival_files_programs <- data.frame(path=fs::dir_ls("../GBM-Model-052022/KM_csv_Files/",recurse = T, type = "file",glob = "*program*.csv"),row.names = NULL)

  # Read survival KM plot files
  survival_read_programs <- tibble(file_names = survival_files_programs$path) %>%
    mutate(data = map(survival_files_programs$path, read_csv, show_col_types = FALSE)) %>%
    unnest()

  # create and modify tibble
  survival_df_programs <- survival_read_programs |>
    separate(file_names, into=c("a","b","c","filename"), sep="/", remove = F) |>
    separate(filename, into = c("program","d","e","type"), sep = "_") |>
    dplyr::select(-file_names, -a,-b,-c,-d,-e) |>
    mutate(color = case_when(type == "min" ~ "#DD9900",
                             type == "max" ~ "#33DDEE",
                             TRUE ~ "#DDDDDD"))
}
