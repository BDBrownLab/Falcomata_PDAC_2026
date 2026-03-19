library('tidyverse')
library('readxl')
library("writexl")
theme_set(theme_minimal(base_family = "Arial"))

library(extrafont)
font_import()
loadfonts(device = "pdf")
library(ggbreak)

#Load data
setwd('../yourpath/')
df_merged <- read.csv('../yourpath/XenaHubs.TcgaTargetGtex.normalized.metadata.csv')
df_merged %>%
  filter(X_study %in% c("TCGA")) -> df_merged
df_merged$Cancer.Type <- df_merged$detailed_category

df_merged %>%
  filter(!X_sample_type %in% c("Solid Tissue Normal")) -> df_merged

cancer_mapping <- c(
  "Brain Lower Grade Glioma" = "LGG", 
  "Liver Hepatocellular Carcinoma" = "LIHC",
  "Cervical & Endocervical Cancer" = "CESC", 
  "Lung Adenocarcinoma" = "LUAD",
  "Colon Adenocarcinoma" = "COAD", 
  "Acute Myeloid Leukemia" = "LAML",
  "Breast Invasive Carcinoma" = "BRCA", 
  "Testicular Germ Cell Tumor" = "TGCT",
  "Sarcoma" = "SARC",
  "Kidney Papillary Cell Carcinoma" = "KIRP",
  "Stomach Adenocarcinoma" = "STAD", 
  "Prostate Adenocarcinoma" = "PRAD",
  "Esophageal Carcinoma" = "ESCA",
  "Skin Cutaneous Melanoma" = "SKCM",
  "Head & Neck Squamous Cell Carcinoma" = "HNSC",
  "Glioblastoma Multiforme" = "GBM",
  "Kidney Clear Cell Carcinoma" = "KIRC",
  "Uterine Corpus Endometrioid Carcinoma" = "UCEC",
  "Thymoma" = "THYM",
  "Bladder Urothelial Carcinoma" = "BLCA",
  "Lung Squamous Cell Carcinoma" = "LUSC",
  "Thyroid Carcinoma" = "THCA",
  "Mesothelioma" = "MESO",
  "Rectum Adenocarcinoma" = "READ",
  "Pancreatic Adenocarcinoma" = "PAAD",
  "Ovarian Serous Cystadenocarcinoma" = "OV",
  "Pheochromocytoma & Paraganglioma" = "PCPG",
  "Uveal Melanoma" = "UVM",
  "Uterine Carcinosarcoma" = "UCS",
  "Kidney Chromophobe" = "KICH",
  "Diffuse Large B-Cell Lymphoma" = "DLBC",
  "Adrenocortical Cancer" = "ACC",
  "Cholangiocarcinoma" = "CHOL"
)

df_merged <- df_merged %>%
  mutate(Cancer.Type = recode(Cancer.Type, !!!cancer_mapping))

df_merged <- df_merged %>%
  pivot_longer(
    cols = c(SERPINE1, SERPINB2),
    names_to = "Gene.Name",
    values_to = "RSEM_norm_count"
  )

df_merged$daysToEvent <- df_merged$OS.time
df_merged$event <- df_merged$OS

df_merged <- df_merged[!is.na(df_merged$Cancer.Type),]

###SERPINB2-----
library(dplyr)
library(broom)
library(survival) 

calculate_hazard_ratio <- function(cancer_type, gene_name, df) {
  df_filtered <- df %>%
    filter(Cancer.Type == cancer_type) %>%
    filter(Gene.Name == gene_name)
  
  duplicated_case_ids <- df_filtered$sample[duplicated(df_filtered$sample)]
  df_filtered <- df_filtered %>%
    filter(!sample %in% duplicated_case_ids)
  
  if (nrow(df_filtered) < 10) {
    return(data.frame(Cancer.Type = cancer_type, Gene.Name = gene_name, HR = NA, p.value = NA))
  }
  
  cutoff <- median(df_filtered$RSEM_norm_count, na.rm = TRUE)
  df_filtered$group <- ifelse(df_filtered$RSEM_norm_count > (cutoff), "SERPINE1_hi", "SERPINE1_low")
  df_filtered$group <- factor(df_filtered$group, levels = c("SERPINE1_low","SERPINE1_hi"))
  model <- coxph(Surv(daysToEvent, event) ~ group, data = df_filtered)
  summary <- tidy(model)
  
  hazard_ratio <- exp(summary$estimate[1])
  p_value <- summary$p.value[1]
  
  return(data.frame(Cancer.Type = cancer_type, Gene.Name = gene_name, HR = hazard_ratio, p.value = p_value))
}



cancer_types <- unique(df_merged$Cancer.Type)
results <- lapply(cancer_types, function(cancer) calculate_hazard_ratio(cancer, "SERPINB2", df_merged))

results_df <- do.call(rbind, results)
results_df <- results_df[order(-results_df$HR), , drop = FALSE]
results_df$significant <- ifelse(results_df$p.value < 0.05, "significant", "not significant")
print(results_df)

results_df$Cancer.Type <- factor(results_df$Cancer.Type, levels = results_df$Cancer.Type[order(results_df$HR)])
ggplot(results_df, aes(x = HR, y = Cancer.Type, fill=significant)) +
  geom_vline(xintercept=1, color="black")+
  geom_segment(aes(x = 1, xend = HR, y = Cancer.Type, yend = Cancer.Type)) +
  geom_point(shape = 21, size = 3) +
  labs(
    x = "HR",
    title = "SERPINB2"
  ) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=c("lightgrey", "#C87D79")) +
  theme(axis.text = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.line.x.bottom = element_line(linewidth = 0.4),
        axis.line.y.left = element_line(linewidth = 0.4),
        text = element_text(family = "Arial"),
        legend.text = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
setwd('../yourpath/')
ggsave(paste0("../yourpath/HR.coxph.SERPINB2.median.png"), plot = last_plot(), dpi = 500, width = 3.5, height = 7)
ggsave(paste0("../yourpath/HR.coxph.SERPINB2.median.pdf"), plot = last_plot(), dpi = 500, width = 3.5, height = 7)

###SERPINE1-----
library(dplyr)
library(broom)
library(survival) 

calculate_hazard_ratio <- function(cancer_type, gene_name, df) {
  df_filtered <- df %>%
    filter(Cancer.Type == cancer_type) %>%
    filter(Gene.Name == gene_name)
  
  duplicated_case_ids <- df_filtered$sample[duplicated(df_filtered$sample)]
  df_filtered <- df_filtered %>%
    filter(!sample %in% duplicated_case_ids)
  
  if (nrow(df_filtered) < 10) {
    return(data.frame(Cancer.Type = cancer_type, Gene.Name = gene_name, HR = NA, p.value = NA))
  }
  
  cutoff <- median(df_filtered$RSEM_norm_count, na.rm = TRUE)
  df_filtered$group <- ifelse(df_filtered$RSEM_norm_count > (cutoff), "SERPINE1_hi", "SERPINE1_low")
  df_filtered$group <- factor(df_filtered$group, levels = c("SERPINE1_low","SERPINE1_hi"))
  model <- coxph(Surv(daysToEvent, event) ~ group, data = df_filtered)
  summary <- tidy(model)
  
  hazard_ratio <- exp(summary$estimate[1])
  p_value <- summary$p.value[1]
  
  return(data.frame(Cancer.Type = cancer_type, Gene.Name = gene_name, HR = hazard_ratio, p.value = p_value))
}

cancer_types <- unique(df_merged$Cancer.Type)
results <- lapply(cancer_types, function(cancer) calculate_hazard_ratio(cancer, "SERPINE1", df_merged))

results_df <- do.call(rbind, results)
results_df <- results_df[order(-results_df$HR), , drop = FALSE]
results_df$significant <- ifelse(results_df$p.value < 0.05, "significant", "not significant")
print(results_df)

results_df$Cancer.Type <- factor(results_df$Cancer.Type, levels = results_df$Cancer.Type[order(results_df$HR)])
ggplot(results_df, aes(x = HR, y = Cancer.Type, fill=significant)) +
  geom_vline(xintercept=1, color="black")+
  
  geom_segment(aes(x = 1, xend = HR, y = Cancer.Type, yend = Cancer.Type)) +
  geom_point(shape = 21, size = 3) +
  labs(
    x = "HR",
    title = "SERPINE1"
  ) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=c("lightgrey", "#7F61A6")) +
  theme(axis.text = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.line.x.bottom = element_line(linewidth = 0.4),
        axis.line.y.left = element_line(linewidth = 0.4),
        text = element_text(family = "Arial"),
        legend.text = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,),
        panel.grid = element_blank(),
        panel.border = element_blank(),      
        panel.spacing = unit(0, "lines"), 
        panel.background = element_blank() )
setwd('../yourpath/')
ggsave(paste0("../yourpath/HR.coxph.SERPINE1.median.png"), plot = last_plot(), dpi = 500, width = 5.5, height = 7)
ggsave(paste0("../yourpath/HR.coxph.SERPINE1.median.pdf"), plot = last_plot(), dpi = 500, width = 5.5, height = 7)



         