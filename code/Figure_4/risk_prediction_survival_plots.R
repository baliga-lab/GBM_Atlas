# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2024
# Email: sturkarslan@systemsbiology.org
#
# Date: 2024-05-21
#
# Script Name: risk_prediction_survival_plots.R
#
# Script Description: Plots survival plots for various datasets
#
# Script version: 
# --------------------------------------------


library(tidyverse)
library(fs)
library(ggsurvfit)
library(survival)
library(survRM2)
library(gghighlight)
library(survminer)
library(ggbreak) 
library(gridExtra)
library(ggpmisc)

## Load xCures CARIS CSV files
data_files <- fs::dir_ls("../XCures-SNO/data/Caris_CSVs_for_ISB_analysis",recurse = T,regexp = ".csv")
#rnaseq_files <- fs::dir_ls(root_folder,recurse = 1)
# reformat files into a table
complete_data_files <- data.frame(path= data_files,row.names = NULL) |>
  separate(path, into=c("d1","d2","d3","d4","patient_id"), sep = "/", remove = F) |>
  separate(patient_id, into=c("patient_id","case_id","date"), sep = "_", remove = F)

# read and collate program activities for LGG cohort
xcures_caris <- map_dfr(complete_data_files$path, read_csv) 

# plot IDH mutation status
xcures_idh <- xcures_caris |>
  #filter(Technology == "Exome") |>
  dplyr::select(`Accession Number`, Biomarker, `Test Result`, Technology) |>
  filter(Biomarker %in% c("IDH1")) |>
  filter(Technology %in% c("NGS Q3","Exome"))|>
  mutate(TestResult = if_else(`Test Result` %in% c(grep("^(W|w)ild",`Test Result`, value=TRUE, perl=TRUE)), "Wild Type", `Test Result`)) |>
  unique()

xcures_mgmt <- xcures_caris |>
  dplyr::select(`Accession Number`, Biomarker, Test,`Test Result`, Technology) |>
  filter(Biomarker %in% c("MGMT")) |>
  filter(Test == "MGMT-Me") |>
  unique() |>
  select("Accession Number", `Test Result`) |>
  dplyr::rename("TestResult_MGMT" = "Test Result")

# Load risk prediction results for xCures new analysis
gene_model_x <- read_csv("data/risk_prediction/gene_model_km_raw_data_xcure_10_2_2023.csv") |>
  mutate(type = "gene_model")

regulon_model_x <- read_csv("data/risk_prediction/regulon_model_km_raw_data_xcure_10_2_2023.csv") |>
  mutate(type = "regulon_model")

program_model_x <- read_csv("data/risk_prediction/program_model_km_raw_data_xcure_10_2_2023.csv") |>
  mutate(type = "program_model")

gene_panel_x <- read_csv("data/risk_prediction/gene_panel_km_raw_data_xcure_10_2_2023.csv") |>
  mutate(type = "gene_panel")|>
  dplyr::rename("Patient_ID" = "...1")


# Load risk prediction results
gene_model <- read_csv("data/risk_prediction/gene_model_km_raw_data.csv") |>
  mutate(type = "gene_model") |>
  filter(dataset != "xCures") |>
  bind_rows(gene_model_x) |> 
  left_join(xcures_idh, by=c("Patient_ID" = "Accession Number")) |>
  left_join(xcures_mgmt, by=c("Patient_ID" = "Accession Number"))
  

regulon_model <- read_csv("data/risk_prediction/regulon_model_km_raw_data.csv") |>
  mutate(type = "regulon_model")|>
  filter(dataset != "xCures") |>
  bind_rows(regulon_model_x) |>
  left_join(xcures_idh, by=c("Patient_ID" = "Accession Number")) |>
  left_join(xcures_mgmt, by=c("Patient_ID" = "Accession Number"))

program_model <- read_csv("data/risk_prediction/program_model_km_raw_data.csv") |>
  mutate(type = "program_model")|>
  filter(dataset != "xCures") |>
  bind_rows(program_model_x) |>
  left_join(xcures_idh, by=c("Patient_ID" = "Accession Number")) |>
  left_join(xcures_mgmt, by=c("Patient_ID" = "Accession Number"))

gene_panel <- read_csv("data/risk_prediction/gene_panel_km_raw_data.csv") |>
  mutate(type = "gene_panel")|>
  filter(dataset != "xCures") |>
  bind_rows(gene_panel_x) |>
  left_join(xcures_idh, by=c("Patient_ID" = "Accession Number")) |>
  left_join(xcures_mgmt, by=c("Patient_ID" = "Accession Number"))



# Combine all predictions into one file for supplementry table
m1 <- bind_rows(gene_panel,gene_model, regulon_model,program_model) |>
  mutate(type = factor(type, levels=c("gene_panel","gene_model","regulon_model","program_model")))
write_csv(m1, file="output/Risk_Prediction_All-Models.csv")

km_plot_function <- function(df, in.medians, type, dataset){
  n_high <- df |>
    filter(riskLabel == "High") |>
    pull(Patient_ID) |>
    unique() |>
    length()
  
  n_low <- df |>
    filter(riskLabel == "Low") |>
    pull(Patient_ID) |>
    unique() |>
    length()
  
  p01 <- survfit2(Surv(duration, observed) ~ riskLabel, data = df) |>
    ggsurvfit(linewidth = 2) +
    add_quantile(y_value = 0.5, color = c("#d73027","#4575b4","#333333"), linewidth = 0.75) +
    scale_ggsurvfit() +
    add_censor_mark(size = 4, alpha = 0.9) +
    scale_color_manual(breaks=c("Low","High"), values=c("#4575b4","#d73027")) +
    add_pvalue(location = "annotation", x = 1250,y=0.95, hjust=0, size=6) +
    geom_text(data=in.medians, aes(x=median, y=0, label=paste0("mOS: ",round(median, digits = 0))), size=5, color=c("#d73027","#4575b4"), fontface="bold", hjust=c(1.05,-0.1)) +
    geom_text(label=paste0("n= ", n_high), x = 1250,y=0.85, hjust=0, size=5, color="#d73027") +
    geom_text(label=paste0("n= ", n_low), x = 1250,y=0.75, hjust=0, size=5,color="#4575b4") +
    labs(title=paste0(type, " | ", dataset), x="Time (days)") +
    #scale_x_continuous(limits = c(0,2000), breaks = c(0,1000,2000)) +
    theme(axis.text = element_text(size=16), axis.title = element_text(size=16), legend.position = "none") 
  p01
}

#### Function to create KM curves
xcures_km_plots <- function(in.type, in.dataset){
  # get the appropariate dataframe for KM curve
  df1 <- m1 |>
    filter(dataset == in.dataset & type == in.type)
  # Build survival model
  tmp1 <- survfit2(Surv(duration, observed) ~ riskLabel, data = df1 )
  # Get median values for plotting
  new_medians <- surv_median(tmp1)
  
  # Gene Panel Plots
  pp1 <- km_plot_function(df = df1, in.medians = new_medians, type = in.type, dataset = in.dataset)
  pp1
}


grid.arrange(
  xcures_km_plots(in.type = "gene_panel", in.dataset = "Gravendeel"),
  xcures_km_plots(in.type = "gene_panel", in.dataset = "TCGA"),
  xcures_km_plots(in.type = "gene_panel", in.dataset = "xCures"),
  xcures_km_plots(in.type = "program_model", in.dataset = "Gravendeel"),
  xcures_km_plots(in.type = "program_model", in.dataset = "TCGA"),
  xcures_km_plots(in.type = "program_model", in.dataset = "xCures"),
  ncol=3,
  nrow = 2
  
)







p01 <- survfit2(Surv(duration, observed) ~ riskLabel, data = m1 |> dplyr::filter(dataset=="xCures" & type=="program_model")) |>
  ggsurvfit(linewidth = 1) +
  add_quantile(y_value = 0.5, color = c("#d73027","#4575b4","#333333"), linewidth = 0.75) +
  scale_ggsurvfit() +
  add_censor_mark(size = 4, alpha = 0.9) +
  scale_color_manual(breaks=c("Low","High"), values=c("#4575b4","#d73027")) +
  add_pvalue(location = "annotation", x = 1,y=0.75, hjust=0, size=5) +
  geom_text(data=medians, aes(x=median, y=0, label=median),color=c("#d73027","#4575b4"), hjust=c(5,1)) +
  geom_text(label=paste0("n= ", length(xcures_all$Patient_ID)), x = 1,y=0.65, hjust=0, size=5) +
  labs(title="xCures // IDH Status", x="Time (days)") 


p <- ggsurvplot(d1, data = as.data.frame(m1),
           combine = FALSE,
           facet.by = c("dataset","type"),
           censor = FALSE,
           legend.title = "Risk",
           legend.labs = c("High","Low"),
           pval = TRUE,
           palette = c("#d73027","#4575b4"),
           xlim = c(0,2000),
           ggtheme = theme_bw()
           )
p+add_quantile(y_value = 0.5, color = c("#d73027","#333333"), linewidth = 0.75) 

risk_summary <- d1 |> tidy_survfit()
write_csv(risk_summary, file="output/Survival_Analysis_Output.csv")


################# Only Program Model #######################
# Get counts of IDH distribution
pm_counts <- program_model |>
  filter(dataset=="xCures") |>
  select(Patient_ID, riskLabel,TestResult) |>
  group_by(riskLabel,TestResult) |>
  mutate(TestResult = case_when(TestResult == "Pathogenic Variant" ~ "IDHmut",
                                TestResult == "Wild Type" ~ "IDHwt",
                                TestResult == NA ~ "NA")) |>
  summarise(count=n()) |>
  ungroup() |>
  group_by(riskLabel) |>
  mutate(total = sum(count)) |>
  mutate(percent = round((count/total)*100, digits = 2))

highrisk_label <- paste0("IDHwt: %",
  pm_counts |> 
    filter(riskLabel=="High" & TestResult == "IDHwt") |> 
    pull(percent),
  "\n IDHmut: %",
  pm_counts |> 
    filter(riskLabel=="High" & TestResult == "IDHmut") |> 
    pull(percent),
  "\n Unknown: %",
  pm_counts |> 
    filter(riskLabel=="High" & is.na(TestResult)) |> 
    pull(percent)
  )

lowrisk_label <- paste0("IDHwt: %",
                         pm_counts |> 
                           filter(riskLabel=="Low" & TestResult == "IDHwt") |> 
                           pull(percent),
                         "\n IDHmut: %",
                         pm_counts |> 
                           filter(riskLabel=="Low" & TestResult == "IDHmut") |> 
                           pull(percent),
                         "\n Unknown: %",
                         pm_counts |> 
                           filter(riskLabel=="Low" & is.na(TestResult)) |> 
                           pull(percent)
)

d3 <- survfit2(Surv(duration, observed) ~ riskLabel, data = program_model |> filter(dataset=="xCures"))
p <- ggsurvplot(d3, data = as.data.frame(program_model),
                combine = FALSE,
                #facet.by = c("dataset"),
                censor = FALSE,
                legend.title = "Risk",
                legend.labs = c("High","Low"),
                pval = TRUE,
                palette = c("#d73027","#4575b4"),risk.table = TRUE,
                xlim = c(0,2000),
                ggtheme = theme_bw()
)
p$plot <- p$plot + 
  annotate(geom="text", label=list(highrisk_label), x=100, y=0.5, color="#d73027") +
  annotate(geom="text", label=list(lowrisk_label), x=1200, y=0.5, color="#4575b4") +
  labs(title="Program model: xCures")


##### IDH mutation status based
d4 <- survfit2(Surv(duration, observed) ~ TestResult,  data = program_model |> filter(dataset=="xCures" & !is.na(TestResult)))
p <- ggsurvplot(d4, data = as.data.frame(program_model |> filter(dataset=="xCures" & !is.na(TestResult))),
                combine = FALSE,
                #facet.by = c("dataset"),
                censor = FALSE,
                legend.title = "Risk",
                #legend.labs = c("High","Low"),
                pval = TRUE,
                palette = c("#d73027","#4575b4"),risk.table = TRUE,
                xlim = c(0,2000),
                ggtheme = theme_bw()
)
p$plot <- p$plot + 
  annotate(geom="text", label=list(highrisk_label), x=100, y=0.5, color="#d73027") +
  annotate(geom="text", label=list(lowrisk_label), x=1200, y=0.5, color="#4575b4") +
  labs(title="Program model: xCures")


##### MGMT Methylation status based
# Survival based on just MGMT metylation
d5 <- survfit2(Surv(duration, observed) ~ TestResult_MGMT, data = program_model |> filter(dataset=="xCures" & TestResult_MGMT != "Equivocal"))
p <- ggsurvplot(d5, data = as.data.frame(program_model |> filter(dataset=="xCures" & TestResult_MGMT !="Equivocal")),
                combine = FALSE,
                #facet.by = c("dataset"),
                censor = FALSE,
                legend.title = "Risk",
                legend.labs = c("Hypermethylated","Not_Metylated"),
                pval = TRUE,
                palette = c("#4575b4","#d73027"),
                xlim = c(0,2000),risk.table = TRUE,
                ggtheme = theme_bw()
)
p$plot <- p$plot + 
 # annotate(geom="text", label=list(highrisk_label), x=100, y=0.5, color="#d73027") +
  #annotate(geom="text", label=list(lowrisk_label), x=1200, y=0.5, color="#4575b4") +
  labs(title="xCures | MGMT Methylation")


# ################# Survival based on just Program model vs MGMT metylation

# Get counts of MGMT Metylation distribution
pm_counts_mgmt <- program_model |>
  filter(dataset=="xCures") |>
  select(Patient_ID, riskLabel,TestResult_MGMT) |>
  group_by(riskLabel,TestResult_MGMT) |>
  summarise(count=n()) |>
  ungroup() |>
  group_by(riskLabel) |>
  mutate(total = sum(count)) |>
  mutate(percent = round((count/total)*100, digits = 2))

highrisk_label_mgmt <- paste0("Not_Methylated: %",
                         pm_counts_mgmt |> 
                           filter(riskLabel=="High" & TestResult_MGMT == "Not Methylated") |> 
                           pull(percent),
                         "\n Hypermethylated: %",
                         pm_counts_mgmt |> 
                           filter(riskLabel=="High" & TestResult_MGMT == "Hypermethylated") |> 
                           pull(percent),
                         "\n Equivocal %",
                         pm_counts_mgmt |> 
                           filter(riskLabel=="High" & TestResult_MGMT == "Equivocal") |> 
                           pull(percent),
                         "\n Unknown: %",
                         pm_counts_mgmt |> 
                           filter(riskLabel=="High" & is.na(TestResult_MGMT)) |> 
                           pull(percent)
)

lowrisk_label_mgmt <- paste0("Not_Methylated: %",
                              pm_counts_mgmt |> 
                                filter(riskLabel=="Low" & TestResult_MGMT == "Not Methylated") |> 
                                pull(percent),
                              "\n Hypermethylated: %",
                              pm_counts_mgmt |> 
                                filter(riskLabel=="Low" & TestResult_MGMT == "Hypermethylated") |> 
                                pull(percent),
                              "\n Equivocal %",
                              pm_counts_mgmt |> 
                                filter(riskLabel=="Low" & TestResult_MGMT == "Equivocal") |> 
                                pull(percent),
                              "\n Unknown: %",
                              pm_counts_mgmt |> 
                                filter(riskLabel=="Low" & is.na(TestResult_MGMT)) |> 
                                pull(percent)
)




d6 <- survfit2(Surv(duration, observed) ~ riskLabel, data = program_model |> filter(dataset=="xCures"))
p <- ggsurvplot(d6, data = as.data.frame(program_model |> filter(dataset=="xCures")),
                combine = FALSE,
                #facet.by = c("dataset"),
                censor = FALSE,
                legend.title = "Risk",
                legend.labs = c("High","Low"),
                pval = TRUE,
                palette = c("#d73027","#4575b4"),
                risk.table = TRUE,
                xlim = c(0,2000),
                ggtheme = theme_bw()
)
p$plot <- p$plot + 
  annotate(geom="text", label=list(highrisk_label_mgmt), x=-100, y=0.5, color="#d73027", hjust=0) +
  annotate(geom="text", label=list(lowrisk_label_mgmt), x=1400, y=0.5, color="#4575b4", hjust=0) +
  labs(title="Program model: xCures | MGMT Methylation")

  

  ###################### TCGA Survival ###################
  # Load TCGA new subtypes
  tcga_subtypes <- read_csv("data/Zakharova_etal_Glioma_reclassification.csv") |>
    mutate(dataset = "TCGA")
  
  # combine new annotations with program model survival
  program_model_tcga <- program_model |>
    filter(dataset == "TCGA") |>
    left_join(tcga_subtypes, by=c("Patient_ID" = "Case_ID")) |>
    mutate(WHO_CNS5_diagnosis = gsub(",", "\n", WHO_CNS5_diagnosis)) |>
    mutate(WHO_CNS5_diagnosis = if_else(WHO_CNS5_diagnosis == "NA. Grade NA.", NA, WHO_CNS5_diagnosis))

  # combine new annotations with program model survival
  regulon_model_tcga <- regulon_model |>
    filter(dataset == "TCGA") |>
    left_join(tcga_subtypes, by=c("Patient_ID" = "Case_ID")) |>
    mutate(WHO_CNS5_diagnosis = gsub(",", "\n", WHO_CNS5_diagnosis)) |>
    mutate(WHO_CNS5_diagnosis = if_else(WHO_CNS5_diagnosis == "NA. Grade NA.", NA, WHO_CNS5_diagnosis))
  
  ## Function to draw 
  survival_counts <- function(input_df, type){
    # Get counts of IDH distribution
    if(type == "IDH"){
      pm_counts_tcga <- input_df |>
        select(Patient_ID, riskLabel, `IDH status`) |>
        group_by(riskLabel,`IDH status`) |>
        mutate(`IDH status` = if_else(`IDH status` == "WT (maf)", "WT", `IDH status`)) |>
        summarise(count=n()) |>
        ungroup() |>
        group_by(riskLabel) |>
        mutate(total = sum(count)) |>
        mutate(percent = round((count/total)*100, digits = 2))
      
      highrisk_label_tcga <- paste0("IDHwt: %",
                                    pm_counts_tcga |> 
                                      filter(riskLabel=="High" & `IDH status` == "WT") |> 
                                      pull(percent),
                                    "\n IDHmut: %",
                                    pm_counts_tcga |> 
                                      filter(riskLabel=="High" & `IDH status` == "Mutant") |> 
                                      pull(percent),
                                    "\n Unknown: %",
                                    pm_counts_tcga |> 
                                      filter(riskLabel=="High" & `IDH status` == "unknown") |> 
                                      pull(percent)
      )
      
      lowrisk_label_tcga <- paste0("IDHwt: %",
                                   pm_counts_tcga |> 
                                     filter(riskLabel=="Low" & `IDH status` == "WT") |> 
                                     pull(percent),
                                   "\n IDHmut: %",
                                   pm_counts_tcga |> 
                                     filter(riskLabel=="Low" & `IDH status` == "Mutant") |> 
                                     pull(percent),
                                   "\n Unknown: %",
                                   pm_counts_tcga |> 
                                     filter(riskLabel=="Low" & `IDH status` == "unknown") |> 
                                     pull(percent)
      )
    }
    
    if(type == "MGMT"){
      pm_counts_tcga <- input_df |>
        select(Patient_ID, riskLabel, `MGMT promoter status`) |>
        group_by(riskLabel,`MGMT promoter status`) |>
       # mutate(`MGMT promoter status` = if_else(`IDH status` == "WT (maf)", "WT", `IDH status`)) |>
        summarise(count=n()) |>
        ungroup() |>
        group_by(riskLabel) |>
        mutate(total = sum(count)) |>
        mutate(percent = round((count/total)*100, digits = 2))
      
      highrisk_label_tcga <- paste0("Methylated: %",
                                    pm_counts_tcga |> 
                                      filter(riskLabel=="High" & `MGMT promoter status` == "Methylated") |> 
                                      pull(percent),
                                    "\n Unmethylated: %",
                                    pm_counts_tcga |> 
                                      filter(riskLabel=="High" & `MGMT promoter status` == "Unmethylated") |> 
                                      pull(percent),
                                    "\n Unknown: %",
                                    pm_counts_tcga |> 
                                      filter(riskLabel=="High" & `MGMT promoter status` == "unknown") |> 
                                      pull(percent)
      )
      
      lowrisk_label_tcga <- paste0("Methylated: %",
                                   pm_counts_tcga |> 
                                     filter(riskLabel=="Low" & `MGMT promoter status` == "Methylated") |> 
                                     pull(percent),
                                   "\n Unmethylated: %",
                                   pm_counts_tcga |> 
                                     filter(riskLabel=="Low" & `MGMT promoter status` == "Unmethylated") |> 
                                     pull(percent),
                                   "\n Unknown: %",
                                   pm_counts_tcga |> 
                                     filter(riskLabel=="Low" & `MGMT promoter status` == "unknown") |> 
                                     pull(percent)
      )
    }
    
    return(list(high_risk = highrisk_label_tcga, low_risk = lowrisk_label_tcga, table=pm_counts_tcga))
    
  }
  
  # Create survival labels
  survival_labels <- survival_counts(input_df = program_model_tcga, type = "IDH")

# survival analysis
  d_tcga_idh <- survfit2(Surv(duration, observed) ~ riskLabel, data = program_model_tcga)
  p_tcga_idh <- ggsurvplot(d_tcga_idh, data = as.data.frame(program_model_tcga),
                  combine = FALSE,
                  #facet.by = c("dataset"),
                  censor = FALSE,
                  legend.title = "Risk",
                  legend.labs = c("High","Low"),
                  pval = TRUE,
                  palette = c("#d73027","#4575b4"),
                  risk.table = TRUE,
                  xlim = c(0,2000),
                  ggtheme = theme_bw()
  )
  p_tcga_idh$plot <- p_tcga_idh$plot + 
    annotate(geom="text", label=list(survival_labels$high_risk), x=-100, y=0.5, color="#d73027", hjust=0) +
    annotate(geom="text", label=list(survival_labels$low_risk), x=1400, y=0.5, color="#4575b4", hjust=0) +
    labs(title="Program model: TCGA | IDH Mutation") 
  
 p_tcga_idh
  
  # Create survival labels
  survival_labels <- survival_counts(input_df = program_model_tcga, type = "MGMT")
  
  p_tcga_mgmt <- ggsurvplot(d_tcga_idh, data = as.data.frame(program_model_tcga),
                           combine = FALSE,
                           #facet.by = c("dataset"),
                           censor = FALSE,
                           legend.title = "Risk",
                           legend.labs = c("High","Low"),
                           pval = TRUE,
                           palette = c("#d73027","#4575b4"),
                           risk.table = TRUE,
                           xlim = c(0,2000),
                           ggtheme = theme_bw()
  )
  p_tcga_mgmt$plot <- p_tcga_mgmt$plot + 
    annotate(geom="text", label=list(survival_labels$high_risk), x=700, y=0.5, color="#d73027", hjust=0) +
    annotate(geom="text", label=list(survival_labels$low_risk), x=1400, y=0.5, color="#4575b4", hjust=0) +
    labs(title="Program model: TCGA | MGMT Status")
  
  
  

  
  # TCGA survival IDH wt vs IDH mutant
  tcga_surv <- read_csv("../GBM-Model-052022/guanSurvivalDf_TCGA_GBM.csv")
  
  tcga_all <- tcga_surv |>
    left_join(tcga_subtypes, by=c("Patient_ID" = "Case_ID")) |>
    mutate(WHO_CNS5_diagnosis = gsub(",", "\n", WHO_CNS5_diagnosis)) |>
    mutate(WHO_CNS5_diagnosis = if_else(WHO_CNS5_diagnosis == "NA. Grade NA.", NA, WHO_CNS5_diagnosis)) |>
    mutate(`IDH status` = if_else(`IDH status` == "WT (maf)", "WT", `IDH status`)) |>
    dplyr::rename("IDHStatus" = "IDH status") |>
    dplyr::filter(IDHStatus != "unknown")
  
  
  tcga_raw_counts <- tcga_all |>
    select(Patient_ID, `IDHStatus`) |>
    group_by(`IDHStatus`) |>
    mutate(`IDHStatus` = if_else(`IDHStatus` == "WT (maf)", "WT", `IDHStatus`)) |>
    summarise(count=n()) |>
    ungroup() |>
    #group_by(`IDHStatus`) |>
    mutate(total = sum(count)) |>
    mutate(percent = round((count/total)*100, digits = 2))
  
  # Create survival labels
  survival_labels <- survival_counts(input_df = tcga_all, type = "IDH")
  
  tcga_raw_idh <- survfit2(Surv(duration, observed) ~ IDHStatus, data = tcga_all)
  p_tcga_raw_idh <- ggsurvplot(tcga_raw_idh, data = as.data.frame(tcga_all),
                           combine = FALSE,
                           #facet.by = c("dataset"),
                           censor = FALSE,
                           legend.title = "IDH Status",
                           legend.labs = c("IDHmut","IDHwt"),
                           pval = TRUE,
                           palette = c("#4575b4","#d73027"),
                           risk.table = TRUE,
                           xlim = c(0,4000),
                           ggtheme = theme_bw()
  )
  p_tcga_raw_idh$plot <- p_tcga_raw_idh$plot + 
    annotate(geom="text", label=list(survival_labels$high_risk), x=-100, y=0.5, color="#d73027", hjust=0) +
    annotate(geom="text", label=list(survival_labels$low_risk), x=1400, y=0.5, color="#4575b4", hjust=0) +
    labs(title="Regulon model: TCGA | IDH Mutation")
  
  
  
  ### Analysis specific to IDH status based separation
  # create Ff fpr IDH wild type
  program_model_tcga_idhwt <- program_model_tcga |>
    filter(`IDH status` == "WT")
  
    # Create survival labels
    survival_labels_idhwt <- survival_counts(input_df = program_model_tcga_idhwt, type = "IDH")
    p1 <- survfit2(Surv(duration, observed) ~ riskLabel, data = program_model_tcga_idhwt) |>
      ggsurvfit(linewidth = 1) +
      add_risktable(risktable_stats= "{n.risk} ({cum.event})",risktable_group='risktable_stats') +
      add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
      add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
      scale_ggsurvfit() +
      add_censor_mark(size = 4, alpha = 0.9) +
      annotate(geom="table", label=list(survival_labels_idhwt$table), x=max(program_model_tcga_idhwt$duration), y=1) +
      scale_color_manual(labels=c("High","Low"), values=c("#d73027","#4575b4")) +
      add_pvalue(location = "annotation", x = 1,y=0, hjust=0, size=5) +
      labs(title="Program model: TCGA | IDH wt", x="Time (days)")
    


    # create Ff fpr IDH wild type
    program_model_tcga_idhmut <- program_model_tcga |>
      filter(`IDH status` == "Mutant")
    
    # Create survival labels
    survival_labels_idhmut <- survival_counts(input_df = program_model_tcga_idhmut, type = "IDH")
    
    p2 <- survfit2(Surv(duration, observed) ~ riskLabel, data = program_model_tcga_idhmut) |>
      ggsurvfit(linewidth = 1) +
      add_risktable(risktable_stats= "{n.risk} ({cum.event})",risktable_group='risktable_stats') +
      add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
      add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
      scale_ggsurvfit() +
      add_censor_mark(size = 4, alpha = 0.9) +
      annotate(geom="table", label=list(survival_labels_idhmut$table), x=max(program_model_tcga_idhmut$duration), y=1) +
      scale_color_manual(labels=c("High","Low"), values=c("#d73027","#4575b4")) +
      add_pvalue(location = "annotation", x = 1,y=0, hjust=0, size=5) +
      labs(title="Program model: TCGA | IDH mut", x="Time (days)")


    # Create survival labels
    #survival_labels_all <- survival_counts(input_df = tcga_all, type = "IDH")
    
    p3 <- survfit2(Surv(duration, observed) ~ IDHStatus, data = tcga_all) |>
      ggsurvfit(linewidth = 1) +
      add_risktable(risktable_stats= "{n.risk} ({cum.event})",risktable_group='risktable_stats') +
      add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
      add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
      scale_ggsurvfit() +
      add_censor_mark(size = 4, alpha = 0.9) +
      annotate(geom="table", label=list(tcga_raw_counts), x=max(tcga_all$duration), y=1) +
      scale_color_manual(labels=c("IDHwt","IDHmut"), values=c("#d73027","#4575b4")) +
      add_pvalue(location = "annotation", x = 1,y=0, hjust=0, size=5) +
      labs(title="TCGA | IDH Status", x="Time (days)")
    
    
    #### MGMT
    ### Analysis specific to MGMT status based separation
    # create Ff fpr MGMT Unmthylated type
    program_model_tcga_mgmtunmet <- program_model_tcga |>
      filter(`MGMT promoter status` == "Unmethylated")
    
    # Create survival labels
    survival_labels_mgmtunmet <- survival_counts(input_df = program_model_tcga_mgmtunmet, type = "MGMT")
    m1 <- survfit2(Surv(duration, observed) ~ riskLabel, data = program_model_tcga_mgmtunmet) |>
      ggsurvfit(linewidth = 1) +
      add_risktable(risktable_stats= "{n.risk} ({cum.event})",risktable_group='risktable_stats') +
      add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
      add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
      scale_ggsurvfit() +
      add_censor_mark(size = 4, alpha = 0.9) +
      annotate(geom="table", label=list(survival_labels_mgmtunmet$table), x=max(program_model_tcga_mgmtunmet$duration), y=1) +
      scale_color_manual(labels=c("High","Low"), values=c("#d73027","#4575b4")) +
      add_pvalue(location = "annotation", x = 1,y=0, hjust=0, size=5) +
      labs(title="Program model: TCGA | MGMT Unmethylated", x="Time (days)")
    
    
    program_model_tcga_mgmtmet <- program_model_tcga |>
      filter(`MGMT promoter status` == "Methylated")
    
    survival_labels_mgmtmet <- survival_counts(input_df = program_model_tcga_mgmtmet, type = "MGMT")
    m2 <- survfit2(Surv(duration, observed) ~ riskLabel, data = program_model_tcga_mgmtmet) |>
      ggsurvfit(linewidth = 1) +
      add_risktable(risktable_stats= "{n.risk} ({cum.event})",risktable_group='risktable_stats') +
      add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
      add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
      scale_ggsurvfit() +
      add_censor_mark(size = 4, alpha = 0.9) +
      annotate(geom="table", label=list(survival_labels_mgmtmet$table), x=max(program_model_tcga_mgmtmet$duration), y=1) +
      scale_color_manual(labels=c("High","Low"), values=c("#d73027","#4575b4")) +
      add_pvalue(location = "annotation", x = 1,y=0, hjust=0, size=5) +
      labs(title="Program model: TCGA | MGMT Methylated", x="Time (days)")
    
    
    # Create survival labels
    #survival_labels_all <- survival_counts(input_df = tcga_all, type = "IDH")
    
    m3 <- survfit2(Surv(duration, observed) ~ MGMT_Status, data = tcga_all |> dplyr::rename("MGMT_Status" = "MGMT promoter status") |> filter(MGMT_Status != "unknown")) |>
      ggsurvfit(linewidth = 1) +
      add_risktable(risktable_stats= "{n.risk} ({cum.event})",risktable_group='risktable_stats') +
      add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
      add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
      scale_ggsurvfit() +
      add_censor_mark(size = 4, alpha = 0.9) +
      #annotate(geom="table", label=list(tcga_raw_counts), x=max(tcga_all$duration), y=1) +
      scale_color_manual(labels=c("Methylated","Unmethylated"), values=c("#d73027","#4575b4")) +
      add_pvalue(location = "annotation", x = 1,y=0, hjust=0, size=5) +
      labs(title="TCGA | MGMT Promoter Status", x="Time (days)")
    
    m3 + m1 + m2
    
    
    ## Comparison of IDH status vs SYGNAl prediction of risk
    merged1 <- program_model_tcga |>
      select(Patient_ID, duration, observed, riskLabel, `IDH status`, `MGMT promoter status`, WHO_CNS5_diagnosis) |>
      mutate(`IDH status` = if_else(`IDH status` == "WT (maf)", "WT", `IDH status`)) |>
      mutate(observed_risk = if_else(duration > 500, "Low", "High")) |>
      arrange(duration) |>
      mutate(correct_count_sygnal = if_else(riskLabel == observed_risk, TRUE, FALSE)) |>
      mutate(correct_count_idhstatus = case_when(`IDH status` == "WT" & observed_risk == "High" ~ TRUE,
                                                 `IDH status` == "Mutant" & observed_risk == "Low" ~ TRUE,
                                                 `IDH status` == "unknown"  ~ NA,
                                                 .default = FALSE))
    
    match_counts <- merged1 |>
      select(Patient_ID, observed_risk, correct_count_sygnal, correct_count_idhstatus) |>
      group_by(correct_count_sygnal, correct_count_idhstatus) |>
      summarise(count=n())
    
    labels <- sprintf("italic(R^2) == %.2f", match_counts$count)
    
    pp1 <- ggplot(merged1, aes(x=duration, fill=duration)) +
      geom_histogram(bins=30)
    
    pp2 <- ggplot(merged1, aes(x=factor(Patient_ID, levels=Patient_ID), y=duration, fill=observed_risk)) +
      geom_bar(stat="identity") +
      coord_flip() +
      facet_grid(correct_count_sygnal ~ correct_count_idhstatus,  scales = "fixed") +
      geom_text(x = 2000, y = 35, aes(label = label), data = labels,hjust = 0) +
      #facet_grid(`IDH status` ~ .,  scales = "fixed") +
      theme(axis.text.y = element_blank()) +
      labs()
    pp2
    
    
    ## Get Tail of distr
    

    ###### End analysis #####
    
    









d2 <- survfit2(Surv(duration, observed) ~ riskLabel + dataset, data = program_model)

d1 |> ggsurvfit(linewidth = 1)

d1 <- survfit2(Surv(duration, observed) ~ riskLabel + dataset, data = m1)

d1 |>
  ggsurvfit(linewidth = 1) +
  #add_confidence_interval() +
  #add_risktable() +
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  add_pvalue(location  = "annotation", x = 4000) + 
  facet_grid(. ~ riskLabel)
  

# Customized survival curves
ggsurvplot(d1, data = as.data.frame(m1),combine = T,
           linetype = c("solid", "solid","solid","solid","solid","solid","solid", "solid","solid","solid","solid","solid","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed"),
           #surv.median.line = "hv", # Add medians survival
           
           # Change legends: title & labels
           legend.title = "Risk",
           legend.labs = c("GR-High","TCGA-High","xCures-High","GR-Low","TCGA-Low","xCures-Low"),
           # Add p-value and tervals
           pval = TRUE,
           facet.by = "dataset",
           #conf.int = TRUE,
           # Add risk table
           risk.table = FALSE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           #color = "dataset",
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#fc8d59","#91bfdb","#998ec3","#fc8d59","#91bfdb","#998ec3","#fc8d59","#91bfdb","#998ec3","#fc8d59","#91bfdb","#998ec3","#fc8d59","#91bfdb","#998ec3","#fc8d59","#91bfdb","#998ec3","#fc8d59","#91bfdb","#998ec3","#fc8d59","#91bfdb","#998ec3"),
           ggtheme = theme_bw() # Change ggplot2 theme
)


# Customized survival curves
ggsurvplot(d1, data = as.data.frame(m1),combine = T,
           linetype = c("solid", "solid","solid","solid","solid","solid","solid", "solid","solid","solid","solid","solid","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed","dashed"),
           #surv.median.line = "hv", # Add medians survival
           
           # Change legends: title & labels
           legend.title = "Risk",
           legend.labs = c("High","Low"),
           # Add p-value and tervals
           pval = TRUE,
           facet.by = c("dataset","type"),
           #conf.int = TRUE,
           # Add risk table
           risk.table = FALSE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           #color = "dataset",
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#d73027","#4575b4","#998ec3","#d73027","#4575b4","#998ec3","#d73027","#4575b4","#998ec3","#d73027","#4575b4","#998ec3","#d73027","#4575b4","#998ec3","#d73027","#4575b4","#998ec3","#d73027","#4575b4","#998ec3","#d73027","#4575b4","#998ec3"),
           ggtheme = theme_bw() # Change ggplot2 theme
)






fit <- survfit( Surv(time = duration, event = observed) ~ riskLabel, data = gene_model)
ggsurvplot_facet(fit, gene_model, facet.by = "dataset",
                 palette = "jco", pval = TRUE)



  
  
  labs(x="Time (days)") +
  gghighlight::gghighlight(strata == as.character(1) | strata == as.character(2),
                           max_highlight = 3,
                           use_group_by = F,
                           unhighlighted_colour = "#DDDDDD",
                           calculate_per_facet = TRUE
  )