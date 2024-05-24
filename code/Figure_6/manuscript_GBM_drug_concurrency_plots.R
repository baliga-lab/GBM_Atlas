# --------------------------------------------
# Author: Serdar Turkarslan
# Copyright (c) Institute for Systems Biology, 2024
# Email: sturkarslan@systemsbiology.org
#
# Date: 2024-05-21
#
# Script Name: manuscript_GBM_drug_concurrency_plots.R
#
# Script Description: Plots the concurrency boxplots/violin plots
#                    for PDGSCs and xCures 
#
# Script version:
# --------------------------------------------


library(here)
library(readr)
library(tidyverse)
library(magrittr)
library(ComplexHeatmap)
library(gridExtra)

## put all files in this directory
input_path <- here("data/")
##input_path <- here("output/HTS/manuscript/tmp/serdar/")

## Data needed for drug sensitivity comparison
## prune_IC50: IC50 scores (Parvi) with scores > 1e-3 set to 1e-3 for
##             scaling purposes.
## AUC_scale: AUC scores scaled between 0 and 1 with 0 being strong response
##            and 1 being poor response (same orientation as IC50)
## via_mean: mean viability at highest dose concentration
## reg_mean: regulon drug constrained activity (mean of all drug mapped regulons)
## prog_mean: program drug constrained activity (mean of all drug mapped programs)
## gof: goodness of fit metric from nlpr
## ll4q_IC50 - updated IC50s using James method of quantifying 
dat <- read_csv(paste0(input_path, "GBM_drug_sensitivity_pruned.csv"))

## if TRUE use LL4q else use Parvi's IC50s
do_LL4Q <- TRUE

################ choose which IC50s are used in analysis ##################################

if(do_LL4Q) {
   ## rename column names - using Parvi's IC50s (prune_IC50)
   dat %>% 
      dplyr::select(Compound, pid, ll4q_IC50, reg_mean) %>%
      dplyr::mutate(hts_act=ll4q_IC50, net_act=reg_mean) -> 
      mat
   
   sig_cut <- 30 ## FDR <= 0.05 significant number of concurrencies (>= sig_cut)
   
   ## count values of # of significant drugs (>=28) 
   study_sig_count <- read_csv(paste0(input_path, "binary_ll4q_screen_sig_histogram1e-6.csv"))
   
   ## p-values of likliehood of seeing x number of significant drugs (>28 concurrencies) out of total (63)
   study_pvals <- read_csv(paste0(input_path, "binary_ll4q_table_3_drug_screen_pvals1e-6.csv"))
   
   ## fdrs of likelihood a drug has > x number of concurrencies out of 43 pdgscs
   read_csv(paste0(input_path, "binary_ll4q_table_2_single_drug_fdrs1e-6.csv")) %>%  
      pivot_longer(-class) %>%   
      pivot_wider(names_from = class, values_from = value) ->
      drug_fdrs
   
} else {
   
   dat %>%
      dplyr::select(Compound, pid, prune_IC50, reg_mean) %>%
      dplyr::mutate(hts_act=prune_IC50, net_act=reg_mean) -> 
      mat
   
   sig_cut <- 28 ## FDR <= 0.05 significant number of concurrencies (>= sig_cut)
   
   ## count values of # of significant drugs (>=30) 
   study_sig_count <- read_csv(paste0(input_path, "binary_nocombo_screen_sig_histogram.csv"))
   
   ## p-values of likliehood of seeing x number of significant drugs (>28 concurrencies) out of total (63)
   study_pvals <- read_csv(paste0(input_path, "binary_nocombo_table_3_drug_screen_pvals.csv"))
   
   ## fdrs of likelihood a drug has > x number of concurrencies out of 43 pdgscs
   read_csv(paste0(input_path, "binary_nocombo_table_2_single_drug_fdrs.csv")) %>%  
      pivot_longer(-class) %>%   
      pivot_wider(names_from = class, values_from = value) ->
      drug_fdrs
   
}



## read in drug screen metadata
meta1 <- read_csv(paste0(input_path, "HTS_PDGSC_Drug_Ranking_Table.csv"))
meta1 %>%
   dplyr::mutate(Drug=toupper(Drug)) %>%
   dplyr::mutate(Drug = case_when(
      Drug == "CLONIDINE HCL" ~ "CLONIDINE HYDROCHLORIDE",
      Drug == "EPIRUBICIN HCL" ~ "EPIRUBICIN HYDROCHLORIDE",
      Drug == "FLUPHENAZINE HCL" ~ "FLUPHENAZINE HYDROCHLORIDE",
      Drug == "MITOXANTRONE HCL" ~ "MITOXANTRONE HYDROCHLORIDE",
      Drug == "TRIFLUOPERAZINE HCL" ~ "TRIFLUOPERAZINE",
      Drug == "TOPOTECAN HCL" ~ "TOPOTECAN HYDROCHLORIDE",
      TRUE ~ Drug)) ->
   meta

## cutoffs for response classification 
hts_cut1 <- 1e-6    ## max for all drugs
hts_cut2 <- 1e-07   ## max for Bortezomib only
net_cut1 <- 0       ## cutoff for gbmSYGNAL


## classify samples for concurrency
mat %>%
   dplyr::mutate(hts_class=ifelse(hts_act > hts_cut1, 
                                         "nonresponder", "responder")) %>%
   dplyr::mutate(net_class=ifelse(net_act <= net_cut1, 
                                         "nonresponder", "responder")) %>%
   dplyr::mutate(hts_class=ifelse(Compound=="BORTEZOMIB", 
                                  ifelse(hts_act > hts_cut2, 
                                       "nonresponder", "responder"),
                                  hts_class)) %>% 
   dplyr::mutate(concurrency=ifelse(hts_class==net_class, net_class, 
                                    "noncoherent")) %>%
   dplyr::mutate(concurrency=ifelse(is.na(hts_act), "missing",
                                    concurrency)) ->
   mat1

##write_csv(mat1, here("data/HTS/GBM_drug_sensitivity_responses.csv"))

## transform into table for heatmap
mapp <- pivot_wider(mat1[,c("Compound", "pid", "concurrency")], 
                    names_from="pid", values_from="concurrency")

## assign row names
heat_con <- as.matrix(mapp[,-1])
rownames(heat_con) <- mapp$Compound

## count number of concurrencies for significant drugs
counts <- apply(heat_con, 1, function(x) sum(x=="responder" | 
                                                x=="nonresponder"))
count_dt <- data.frame(Drug=names(counts), concurrency=counts)
rownames(count_dt) <- NULL
count_dt %<>% dplyr::arrange(desc(concurrency))

## number of drugs showing significant number of concurrencies 
sig_drugs <- sum(counts >= sig_cut)
cat("number of drugs >= ", sig_cut, ": ", sig_drugs, "\n\n")

tot_con <- sum(heat_con=="responder" | heat_con=="nonresponder", na.rm=TRUE)
cat("total concurrency: ", tot_con, "\n\n")

## order to put responders on top, nonresponders on bottm
sum_resp <- apply(heat_con, 1, function(x) sum(x=="responder", na.rm=TRUE))
sum_noresp <- apply(heat_con, 1, function(x) sum(x=="nonresponder", na.rm=TRUE))
ordy <- order(-sum_resp, sum_noresp)

## order by patients - number of concordances
sum_nna <- apply(heat_con, 2, function(x) sum(x=="responder" | 
                                                 x=="nonresponder"))
c_ordy <- order(-sum_nna)

## colum annotations
mat1 %>%
   dplyr::group_by(pid) %>%
   dplyr::count(concurrency) %>%
   pivot_wider(names_from=concurrency, values_from=n) %>%
   mutate(across(c(responder, noncoherent, nonresponder, missing), ~replace_na(.x, 0))) %>%
   arrange(match(pid, colnames(heat_con[ordy,c_ordy])))  ->
   col_bar_dat

col_annot <- columnAnnotation(
   coh=anno_barplot(col_bar_dat[,-1], which="row",
                    beside = FALSE,
                    gp = gpar(fill = c("grey" ,"blue","red", "purple"))),
   # responder=anno_simple(as.character(pd_resp), 
   #                         height= unit(1, c('mm')),
   #                         col=c("TRUE"="green", "FALSE"="yellow")),
   # nonrespond=anno_simple(as.character(pd_noresp), 
   #                         height= unit(1, c('mm')),
   #                    col=c("TRUE"="green", "FALSE"="yellow")),
   show_annotation_name=FALSE)

## calculate basic performance (total concurrency and sig genes)
gene_sigs <- apply(heat_con[ordy, c_ordy], 1, function(x) return(sum(x=="responder" | x=="nonresponder") >= sig_cut))
##write_csv(data.frame(drug=names(gene_sigs), response=gene_sigs), 
##          here("output/HTS/figures/drug_response_counts.csv"))

# left row annotations  
mat1 %>%
   dplyr::group_by(Compound) %>%
   dplyr::count(concurrency) %>%
   pivot_wider(names_from=concurrency, values_from=n) %>%
   mutate(across(c(responder, noncoherent, nonresponder, missing), ~replace_na(.x, 0))) %>%
   arrange(match(Compound, rownames(heat_con[ordy,])))  ->
   row_bar_dat

row_annot_left <- rowAnnotation(
   coh=anno_barplot(row_bar_dat[,-1], which="row",
                    beside = F,
                    numbers_gp = gpar(fontsize = 3),
                    gp = gpar(fill = c( "grey" ,"blue","red", "purple"))),
   signficance=anno_simple(as.character(gene_sigs), 
                           width= unit(2, c('mm')),
                           col=c("TRUE"="green", "FALSE"="yellow")),
   show_annotation_name=FALSE)

## create legend for significant annotation
sigs = Legend(title = "Significant",  at = c("TRUE", "FALSE"), 
              legend_gp = gpar(fill = c("green", "yellow")),
              labels = c("yes", "no"))

## right row annotations
meta %>%
   dplyr::filter(!duplicated(Drug)) %>%
   dplyr::select(Drug, recommendation, Class) %>%
   dplyr::filter(Drug %in% rownames(heat_con)) %>%
   arrange(match(Drug, rownames(heat_con[ordy,]))) -> 
   row_anno_right

library(paletteer)
class_names <- sort(unique(row_anno_right$Class))
row_color_list <- c(paletteer_d(`"ggthemes::Tableau_20"`, n = length(class_names)))
names(row_color_list) <- class_names

ha = HeatmapAnnotation(
  ##rec=row_anno_right$recommendation, 
  class = row_anno_right$Class,
  which="row",
  col = as.list(row_color_list))


heatie <- Heatmap(heat_con[ordy,c_ordy],
                  col=c(nonresponder="blue", responder="red", noncoherent="grey", missing="purple"),
                  top_annotation = col_annot,
                  left_annotation = row_annot_left,
                  right_annotation = ha,
                  ##row_order=ordy,
                  ##column_order=c_ordy,
                  heatmap_legend_param = list(
                     title = "coherence",
                     legends_gp = gpar(fontsize = 4)
                     ##grid_height = unit(3, "mm"),
                     ##grid_width = unit(3, "mm")
                  ),
                  show_row_dend = F,
                  show_column_dend = F,
                  rect_gp = gpar(col = "white", lwd = 0.4),
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8),
)
pdf(file="output/PDGSC_Concurrency_Heatmap_0906.pdf",height = 7.73, width = 10.31)
draw(heatie, annotation_legend_list = list(sigs),  merge_legend = TRUE)
dev.off()

################################## plot study p-values ####################################

max_count <- max(study_sig_count$count)
ggplot(study_sig_count, aes(x=count)) + geom_histogram(col="black", fill="cyan3", bins=30) +
   xlim(0,30) +
   xlab(paste0("Number of Drugs(>= ", sig_cut, " Concurrency)")) +
   ylab("Count") +
   geom_vline(xintercept=max_count, linetype="dashed", col="cyan3") +
   geom_vline(xintercept=sig_drugs, linetype="dashed", col="red") 

############################### plot single drug FDRs #######################################

ggplot(as.data.frame(drug_fdrs), aes(x=as.numeric(name), y=total)) + geom_point() +
   xlab("Number of concurrency events (43)") +
   ylab("False Discovery Rate") +
   geom_vline(xintercept=sig_cut, linetype="dotted", color="red") +
   geom_hline(yintercept=0.05, linetype="dotted", color="red") +
   scale_x_continuous(breaks=c(0, 10, 20, 30, 40))  +
   scale_y_continuous(breaks=c(0, 0.05,  0.25, 0.50, 0.75, 1)) +
   theme(axis.text.y = element_text(colour=c("black", "red", "black", "black", "black", "black"))) +
   ## add 1 to x-axis so label is adjacent not on top of line
   geom_text(aes(x=sig_cut+1, y=1, label=sig_cut), col="red") +
   theme(aspect.ratio=1/2)
   
############################### plot significance count for each drug ##########################

## I wasn't sure how to show only theme background below scatter plot and how to format test
## labels so they wouldn't overlap

count_dt <- count_dt %>% 
   dplyr::mutate(Drug=factor(Drug, levels=Drug, ordered=TRUE)) %>%
   dplyr::mutate(significant=concurrency >= sig_cut) 

ggplot(count_dt, aes(x=Drug, y=concurrency)) + geom_point() +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   scale_fill_manual(breaks=c(TRUE, FALSE), values=c("cyan3", "grey"), guide="none") +
   geom_label(aes(label=concurrency, fill=significant), size=2, color="white") +
   geom_vline(xintercept="BORTEZOMIB", linetype="dashed", color="cyan3") +
   xlab("") + ylab("Concurrency Count") +
   theme(aspect.ratio=1/3) 

library(ggpubr)

ggdotchart(count_dt, x = "Drug", y = "concurrency",
           color = "significant",                                # Color by groups
           palette = c("gray","#00AFBB"), # Custom color palette
           sorting = "descending", # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 0.5), # Change segment color and size
           #group = as.factor("significant"),                                # Order by groups
           dot.size = 8,rotate = F,                                 # Large dot size
           label = "concurrency",                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 11, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
) + 
  theme_cleveland() + geom_vline(xintercept = "MITOXANTRONE HYDROCHLORIDE", size=1, lty="dashed",color="#00AFBB")


ggplot(count_dt, aes(x=Drug, y=concurrency)) + 
  geom_bar(aes(fill = significant, size = concurrency), alpha = 0.5) +
  scale_color_manual(values = c("gray", "#00AFBB")) +
  scale_size(range = c(4, 11)) +  # Adjust the range of points size
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()

ggplot(count_dt, aes(x=Drug, y=concurrency)) + 
  geom_point(aes(color=significant), size=4) +
  #geom_area(aes(fill = significant), stat="identity") +
  scale_color_manual(values = c("gray", "#00AFBB")) +
  #scale_size(range = c(4, 11)) +  # Adjust the range of points size
  theme(axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust=1)) #+
  #coord_flip()


############################ plot IC50 ratios ##################################################

## all drugs #########################################################
matt <- mat1
drug <- "all Drugs"

## wilcoxon test (non-parameterc rank test)
wil_test <- wilcox.test(hts_act ~ net_class, data=matt) 
wil_pval <- signif(wil_test$p.value, 2)

## t.test (don't think this is valid with IC50 data)
t_test <- t.test(hts_act ~ net_class, data=matt)
t_pval <- signif(t_test$p.value, 2)

p1 <- ggplot(matt, aes(x=net_class, y=hts_act)) + 
   geom_violin(aes(color=net_class),  draw_quantiles=0.5) +
   scale_y_log10() + scale_color_manual(values=c("blue", "red")) +
   geom_jitter(width=0.1, height=0.0000001, alpha=0.25, aes(color=hts_class)) + 
   geom_hline(yintercept=hts_cut1, linetype="dotted") +
   xlab("predicted response") + ylab("log10 IC50") +
   ggtitle(paste(drug, "Wilcox test: ", wil_pval, "| All"))

#devtools::install_github("psyteachr/introdataviz")
library(introdataviz)
library(ggdist)
p1 <- ggplot(matt, aes(x = factor(net_class), y = hts_act, fill = factor(net_class))) +
  #stat_dist_halfeye(adjust=0.5, justification=-.2, .width=0, point_color = NA) +
  #geom_boxplot(width=.12, alpha=0.5) + 
  #stat_dots(aes(color=net_class))
  #introdataviz::geom_split_violin(alpha = .4, trim = FALSE,draw_quantiles = T, show.legend = T) +
  #geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  geom_jitter(width=0.1, height=0.0000001, alpha=0.25, aes(color=hts_class)) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  scale_x_discrete(name = "Condition", labels = c("Non-word", "Word")) +
  scale_y_continuous(name = "Reaction time (ms)",
                     breaks = seq(200, 800, 100), 
                     limits = c(200, 800)) +
  scale_fill_brewer(palette = "Dark2", name = "Language group") +
  theme_minimal()


p1 <- ggplot(matt, aes(x=net_class, y=hts_act)) + 
  geom_violin(aes(color=net_class),  draw_quantiles=0.5) +
  scale_y_log10() + scale_color_manual(values=c("blue", "red")) +
  geom_jitter(width=0.1, height=0.0000001, aes(color=hts_class)) + 
  geom_hline(yintercept=hts_cut1, linetype="dotted") +
  xlab("predicted response") + ylab("log10 IC50") +
  ggtitle(paste(drug, "Wilcox test: ", wil_pval))



## significant drugs ###################################################
sig_drug_names <- count_dt %>% dplyr::filter(significant==TRUE)
matt %>%
   dplyr::filter(Compound %in% sig_drug_names$Drug) ->
   mat_test

## wilcoxon test (non-parameterc rank test)
wil_test <- wilcox.test(hts_act ~ net_class, data=mat_test) 
wil_pval <- signif(wil_test$p.value, 2)

## t.test (don't think this is valid with IC50 data)
t_test <- t.test(hts_act ~ net_class, data=mat_test)
t_pval <- signif(t_test$p.value, 2)

p2 <- ggplot(mat_test, aes(x=net_class, y=hts_act)) + 
   geom_violin(aes(color=net_class),  draw_quantiles=0.5) +
   scale_y_log10() + scale_color_manual(values=c("blue", "red")) +
   geom_jitter(width=0.1, height=0.0000001, aes(color=hts_class)) + 
   geom_hline(yintercept=hts_cut1, linetype="dotted") +
   xlab("predicted response") + ylab("log10 IC50") +
   ggtitle(paste(drug, "Wilcox test: ", wil_pval, "| Significant"))

grid.arrange(p1, p2, nrow=1)

## drug specific #########################################################
pdf(file="output/Wilcox_test_VioloinPlots_0906.pdf")

for(drug in unique(matt$Compound)){
  cat(drug, "\n")
  
  matt %>%
    dplyr::filter(Compound==drug) ->
    mat_test
  
  ## wilcoxon test (non-parameterc rank test)
  wil_test <- try(wilcox.test(hts_act ~ net_class, data=mat_test))
  wil_pval <- try(signif(wil_test$p.value, 2))
  
  ## t.test (don't think this is valid with IC50 data)
  #t_test <- t.test(hts_act ~ net_class, data=mat_test)
 # t_pval <- signif(t_test$p.value, 2)
  
  p <- ggplot(mat_test, aes(x=net_class, y=hts_act)) + 
    geom_violin(aes(color=net_class),  draw_quantiles=0.5) +
    scale_y_log10() + scale_color_manual(values=c("blue", "red")) +
    geom_jitter(width=0.1, height=0.0000001, aes(color=hts_class)) + 
    geom_hline(yintercept=hts_cut1, linetype="dotted") +
    xlab("predicted response") + ylab("log10 IC50") +
    ggtitle(paste(drug, "Wilcox test: ", wil_pval))
  
  
  print(p)
  
  rm(wil_pval)
  rm(wil_test)
  
}
dev.off()
#drug = "METHOTREXATE"



