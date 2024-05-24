# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2024
# Email: sturkarslan@systemsbiology.org
#
# Date: 2024-05-17
#
# Script Name: new_subtype_plotting.R
#
# Script Description: Plots distribution of new GBM subtypes across states and programs
#                     as described by Zakharova_et_al.
#
# Script version:
# --------------------------------------------

library(jsonlite)

# Load new annotations
tcga2 <- read_csv("data/Zakharova_etal_Glioma_reclassification.csv") |>
  dplyr::select(Case_ID, WHO_CNS5_diagnosis)

# TCGA patients in the model
tcga_model <- read_json("../GBM-Model-052022/transcriptional_states.json", simplifyVector = T)

tcga_patients <- data.frame(patient_id = unique(unlist(tcga_model))) |>
  separate(patient_id, into = c("patient_id", "case_no"), sep="\\.", remove=F) |>
  dplyr::select(patient_id) |>
  unique()

# Combine all
tcga_all <- tcga_patients |>
  left_join(tcga2, by=c("patient_id" = "Case_ID")) |>
  mutate(WHO_CNS5_diagnosis = gsub(",", "\n", WHO_CNS5_diagnosis)) |>
  mutate(WHO_CNS5_diagnosis = if_else(WHO_CNS5_diagnosis == "NA. Grade NA.", NA, WHO_CNS5_diagnosis))



#### Enrichment of prognostic markers in States
# Read state hazard ratios
states_hr <- read_csv("../GBM-Model-052022/State_HR_plot_data.csv") |>
  dplyr::rename("Patient_ID" = "...1") |>
  group_by(label) |>
  mutate(medianGuanScore = median(GuanScore)) # |>
#mutate(label = label - 1) ## States start with number zero so we correct for that



####### Subtype enrichment

states_survival <- states_hr |>
  left_join(tcga_all, by=c("Patient_ID" = "patient_id")) |>
  dplyr::rename("Subtype" = "WHO_CNS5_diagnosis")

subtypes <- states_survival|>
  group_by(Subtype) |>
  summarise(subtype_total=n())

subs <- unique(states_survival$Subtype)[!is.na(unique(states_survival$Subtype))]


states_vs_subtypes <- states_survival |>
  group_by(label,Subtype) |>
  summarise(count=n()) |>
  mutate(Subtype = if_else(is.na(Subtype), "Unknown", Subtype))


p <-
  ggplot(states_vs_subtypes,
         aes(
           x = factor(label, levels = unique(states_hr$label)),
           y = count,
           group = label,
           fill = factor(
             Subtype,
             levels = c(
               "Glioblastoma\n IDH-wildtype. Grade 4.",
               "Astrocytoma\n IDH-mutant. Grade 4.",
               "Oligodendroglioma\n IDH-mutant\n 1p/19q-codeleted. Grade 3.",
               "Unknown"
             )
           )
         )) +
  geom_bar(stat = "identity") +
  labs(x = "States", y = "Count", ) +
  theme(legend.text = element_text(size = 7)) +
  geom_hline(yintercept = 1.30103, lty = "dashed") +
  #scale_fill_brewer(type = "seq",palette = "YlGnBu")
  scale_fill_manual(
    name = "Subtype",
    values = c(
      "Glioblastoma\n IDH-wildtype. Grade 4." = "#fdae61",
      "Astrocytoma\n IDH-mutant. Grade 4." = "#abdda4",
      "Oligodendroglioma\n IDH-mutant\n 1p/19q-codeleted. Grade 3." = "#2b83ba",
      "Unknown" = "gray70"
    )
  )

ggsave(p, filename="../GBM-Model-052022/State_vs_New_Subtypes.pdf", device = "pdf",width = 6.61, height= 1.0, units = "in")
#ggsave(p, filename="../GBM-Model-052022/State_vs_Subtypes.pdf", device = "pdf", width=6.61, height = 1.71, units = "in")



## Plot for LGG
# Load LGG data
# TCGA patients in the model
lgg_model <- read_json("../LGG/LGG-Model-110723/transcriptional_states.json", simplifyVector = T)
lgg_risk <- read_csv("../LGG/LGG-Model-110723/guanSurvivalDf_TCGA_LGG.csv")

# Assign states to patients
lgg_model_states <- tibble()
for(i in 1:length(lgg_model)){
  tmp1 <- unlist(lgg_model[[i]])
  df1 <- tibble(patient_id = tmp1) |>
    mutate(label = i)
  lgg_model_states <- bind_rows(lgg_model_states, df1)
}

# Combine all
lgg_all <- lgg_patients |>
  left_join(tcga2, by=c("patient_id" = "Case_ID")) |>
  mutate(WHO_CNS5_diagnosis = gsub(",", "\n", WHO_CNS5_diagnosis)) |>
  mutate(WHO_CNS5_diagnosis = if_else(WHO_CNS5_diagnosis == "NA. Grade NA.", NA, WHO_CNS5_diagnosis))

#### Enrichment of prognostic markers in States
# Read state hazard ratios
states_hr_lgg <- lgg_model_states |>
  mutate(patient_id = gsub("-01","",patient_id)) |>
  dplyr::rename("Patient_ID" = "patient_id") |>
  left_join(lgg_risk, by=c("Patient_ID" = "bcr_patient_barcode")) |>
  group_by(label) |>
  mutate(medianGuanScore = median(GuanScore)) # |>
#mutate(label = label - 1) ## States start with number zero so we correct for that

####### Subtype enrichment

states_survival_lgg <- states_hr_lgg |>
  left_join(lgg_all, by=c("Patient_ID" = "patient_id")) |>
  dplyr::rename("Subtype" = "WHO_CNS5_diagnosis")

subtypes_lgg<- states_survival_lgg |>
  group_by(Subtype) |>
  summarise(subtype_total=n())

subs <- unique(states_survival$Subtype)[!is.na(unique(states_survival$Subtype))]


states_vs_subtypes_lgg <- states_survival_lgg |>
  group_by(label,Subtype) |>
  summarise(count=n()) |>
  mutate(Subtype = if_else(is.na(Subtype), "Unknown", Subtype))


p <-
  ggplot(states_vs_subtypes_lgg,
         aes(
           x = factor(label, levels = unique(states_hr$label)),
           y = count,
           group = label,
           fill = factor(
             Subtype,
             levels = c(
               "Glioblastoma\n IDH-wildtype. Grade 4.",
               "Astrocytoma\n IDH-mutant. Grade 4.",
               "Oligodendroglioma\n IDH-mutant\n 1p/19q-codeleted. Grade 3.",
               "Unknown"
             )
           )
         )) +
  geom_bar(stat = "identity") +
  labs(x = "States", y = "Count", ) +
  theme(legend.text = element_text(size = 7)) +
  geom_hline(yintercept = 1.30103, lty = "dashed") +
  #scale_fill_brewer(type = "seq",palette = "YlGnBu")
  scale_fill_manual(
    name = "Subtype",
    values = c(
      "Glioblastoma\n IDH-wildtype. Grade 4." = "#fdae61",
      "Astrocytoma\n IDH-mutant. Grade 4." = "#abdda4",
      "Oligodendroglioma\n IDH-mutant\n 1p/19q-codeleted. Grade 3." = "#2b83ba",
      "Unknown" = "gray70"
    )
  )

ggsave(p, filename="../GBM-Model-052022/State_vs_New_Subtypes.pdf", device = "pdf",width = 6.61, height= 1.0, units = "in")
#ggsave(p, filename="../GBM-Model-052022/State_vs_Subtypes.pdf", device = "pdf", width=6.61, height = 1.71, units = "in")











### Test enrichment of each TCGA mutation in States
sub_df <- data_frame()
for(my_subtype in subs){
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
    
    res <- phyper(x,m1,n,k1,lower.tail = FALSE)
    
    tmp1 <- data_frame(State = my_state, Subtype = my_subtype, pvalue = res)
    sub_df <- bind_rows(sub_df, tmp1)
    
  }
  
}

## Sbtype enrichment
sub_plot_df <- sub_df |>
  mutate(padjust = p.adjust(pvalue)) |>
  mutate(minusLogpvalue = -log10(padjust))

p <- ggplot(survival, aes(x = factor(State, levels = unique(states_hr$label)), y=minusLogpvalue, group=State, fill=factor(Subtype))) +
  geom_bar(stat="identity") +
  labs(x="States", y="Count",) +
  theme(legend.text = element_text(size=7)) +
  geom_hline(yintercept = 1.30103, lty="dashed") +
  #scale_fill_brewer(type = "seq",palette = "YlGnBu")
  scale_fill_manual(name="Subtype",
                    values = c(
                      "Glioblastoma\n IDH-wildtype. Grade 4." = "#ffffbf",
                      "Astrocytoma\n IDH-mutant. Grade 4." = "#abdda4",
                      "Oligodendroglioma\n IDH-mutant\n 1p/19q-codeleted. Grade 3." = "#2b83ba",
                      "NA" = "#fdae61"
                    )
  )


ggsave(p, filename="../GBM-Model-052022/State_vs_Subtypes.pdf", device = "pdf", width=6.61, height = 1.71, units = "in")

