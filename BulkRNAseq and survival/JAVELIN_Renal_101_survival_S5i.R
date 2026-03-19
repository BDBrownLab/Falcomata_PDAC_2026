library('tidyverse')
library('readxl')
library("writexl")
library(dplyr)
library(broom)
library(survival)
library(survminer)


setwd('../yourpath/')
filename <- "../yourpath/41591_2020_1044_MOESM3_ESM.xlsx"
md_df <- read_xlsx(filename, sheet = "S11_Clinical_data",skip = 1)
counts_df <- read_xlsx(filename, sheet = "S13_Gene_expression_TPM", skip =1)
counts_df %>%
  pivot_longer(
    cols = colnames(counts_df)[2:length(colnames(counts_df))],
    names_to = "ID",
    values_to = "Expression"
  ) -> counts_df

df_motzer <- merge(counts_df, md_df, by="ID")


#SERPINE1
df_motzer_sub <- df_motzer %>%
  filter(HUGO == "SERPINE1") %>%
  filter(TRT01P == "Avelumab+Axitinib")
df_motzer_sub$PFS_P <- as.numeric(df_motzer_sub$PFS_P)
df_motzer_sub$PFS_P_CNSR <- as.numeric(df_motzer_sub$PFS_P_CNSR)

cutoff <- median(df_motzer_sub$Expression)
df_motzer_sub$group <- ifelse(df_motzer_sub$Expression > (cutoff*1), "SERPINE1_hi", "SERPINE1_low")
df_motzer_sub$daysToEvent <- df_motzer_sub$PFS_P
df_motzer_sub$event <- (df_motzer_sub$PFS_P_CNSR-1)*-1
df_motzer_sub$event

surv_obj <- Surv(time = df_motzer_sub$daysToEvent, event = df_motzer_sub$event)
km_fit <- survfit(surv_obj ~ group, data = df_motzer_sub)

ggsurvplot(km_fit, 
           data = df_motzer_sub,
           pval = TRUE,                 
           risk.table = FALSE,           
           xlab = "Months",
           ylab = "PFS",
           title = "SERPINE1",
           legend.labs = c("high", "low"),
           palette = c("#21027D", "#A499C9"),
           size = 0.3,
           censor.size = 2
)  
setwd("../yourpath/")
ggsave(paste0("../yourpath/Motzer2020.KM.SERPINE1.png"), plot = last_plot(), dpi = 500, width = 3.5, height = 4)
ggsave(paste0("../yourpath/Motzer2020.KM.SERPINE1.pdf"), plot = last_plot(), dpi = 500, width = 3.5, height = 4)


#SERPINB2
df_motzer_sub <- df_motzer %>%
  filter(HUGO == "SERPINB2") %>%
  filter(TRT01P == "Avelumab+Axitinib")
df_motzer_sub$PFS_P <- as.numeric(df_motzer_sub$PFS_P)
df_motzer_sub$PFS_P_CNSR <- as.numeric(df_motzer_sub$PFS_P_CNSR)

cutoff <- median(df_motzer_sub$Expression)
df_motzer_sub$group <- ifelse(df_motzer_sub$Expression > (cutoff*1), "high", "low")
df_motzer_sub$daysToEvent <- df_motzer_sub$PFS_P
df_motzer_sub$event <- (df_motzer_sub$PFS_P_CNSR-1)*-1
df_motzer_sub$event

surv_obj <- Surv(time = df_motzer_sub$daysToEvent, event = df_motzer_sub$event)
km_fit <- survfit(surv_obj ~ group, data = df_motzer_sub)

ggsurvplot(km_fit, 
           data = df_motzer_sub,
           pval = TRUE,                 
           risk.table = FALSE,           
           xlab = "Months",
           ylab = "PFS",
           title = "SERPINB2",
           legend.labs = c("high", "low"),
           palette = c("#A72F28", "#CF8C86"),
           size = 0.3,
           censor.size = 2
)  
setwd("../yourpath/")
ggsave(paste0("../yourpath/Motzer2020.KM.SERPINB2.png"), plot = last_plot(), dpi = 500, width = 3.5, height = 4)
ggsave(paste0("../yourpath/Motzer2020.KM.SERPINB2.pdf"), plot = last_plot(), dpi = 500, width = 3.5, height = 4)

