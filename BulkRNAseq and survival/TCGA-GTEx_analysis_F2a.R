library('tidyverse')
library('readxl')
library("writexl")
library("DESeq2")
library("dplyr")
theme_set(theme_minimal(base_family = "Arial"))

#Load TCGA data
setwd('../yourpath/')
df_merged <- read.csv('../yourpath/250122.XenaHubs.TcgaTargetGtex.normalized.metadata.csv')
df_merged %>%
  filter(X_study %in% c("TCGA", "GTEX")) -> df_merged
df_merged$Cancer.Type <- df_merged$detailed_category
df_merged <- df_merged %>%
  pivot_longer(
    cols = c(SERPINE1, SERPINB2),
    names_to = "Gene.Name",
    values_to = "RSEM_norm_count"
  )
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
df_merged <- df_merged[!is.na(df_merged$Cancer.Type),]

#-------------------SERPINE1-------------------



df_TCGA_sub <- df_merged %>%
  filter(Gene.Name == "SERPINE1") %>% 
  filter(X_study == "TCGA") %>%
  filter(!X_sample_type == "Solid Tissue Normal") %>%
  select(sample,RSEM_norm_count, Cancer.Type)
df_TCGA_sub$Cancer.Type <- as.factor(df_TCGA_sub$Cancer.Type)

df_normal_sub <- df_merged %>%
  filter(Gene.Name == "SERPINE1") %>% 
  filter(X_sample_type %in% c("Solid Tissue Normal", "Normal Tissue", 'Cell Line')) %>%
  select(sample,RSEM_norm_count, X_primary_site)

# Mapping
mapping_df <- data.frame(
  Cancer.Type = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", 
                  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", 
                  "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", 
                  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS"),
  X_primary_site = c("Adrenal Gland", "Bladder", "Breast", "Cervix Uteri", "Bile duct", "Colon", 
           "Blood", "Esophagus", "Brain", "Head and Neck region", "Kidney", "Kidney",  "Kidney", "Bone Marrow", "Brain", "Liver",
           "Lung", "Lung",  "Ovary", "Pancreas", "Paraganglia", "Prostate", "Colon", 
           "Soft tissue,Bone", "Skin", "Stomach", "Testis", "Thyroid", "Thymus", "Uterus", "Uterus"),
  stringsAsFactors = FALSE
)

results <- data.frame(cancer.type = character(), 
                      normal.tissue = character(),
                      cancer.RSEM.norm.count.mean = numeric(),
                      normal.RSEM.norm.count.mean = numeric(),
                      log2FC = numeric(),
                      p_value = numeric(), 
                      stringsAsFactors = FALSE)

for (cancer_type in unique(df_TCGA_sub$Cancer.Type)) {
  print(paste("Processing cancer type:", cancer_type))
  tissue_type <- mapping_df %>% filter(Cancer.Type == cancer_type) %>% pull(X_primary_site)
  if (length(tissue_type) == 0) {
    print(paste("No matching tissue type found for cancer type:", cancer_type))
    next
  }
  cancer_RSEM_norm_count <- df_TCGA_sub %>% filter(Cancer.Type == cancer_type) %>% pull(RSEM_norm_count)
  healthy_RSEM_norm_count <- df_normal_sub %>% filter(X_primary_site == tissue_type) %>% pull(RSEM_norm_count)
  
  if (length(healthy_RSEM_norm_count) == 0) {
    print(paste("No matching normal data for tissue type:", tissue_type))
    next
  }
  
  mean_cancer_RSEM_norm_count <- mean(cancer_RSEM_norm_count, na.rm = TRUE)
  mean_healthy_RSEM_norm_count <- mean(healthy_RSEM_norm_count, na.rm = TRUE)
  log2FC <- mean_cancer_RSEM_norm_count - mean_healthy_RSEM_norm_count
  wilcox_test <- wilcox.test(cancer_RSEM_norm_count, healthy_RSEM_norm_count, 
                             alternative = "two.sided", exact = FALSE)
  p_val <- wilcox_test$p.value
  
  results <- rbind(results, data.frame(
    cancer.type = cancer_type,
    normal.tissue = tissue_type,
    cancer.RSEM.norm.count.mean = mean_cancer_RSEM_norm_count,
    normal.RSEM.norm.count.mean = mean_healthy_RSEM_norm_count,
    log2FC = log2FC,
    p_value = p_val
  ))
}
results$p_adj <- p.adjust(results$p_value, method = "BH")  
print(results)

results$significant <- ifelse(results$p_adj < 0.05, "significant", "not significant")
results <- results[order(-results$log2FC), , drop = FALSE]
results$cancer.type <- factor(results$cancer.type, levels = results$cancer.type[order(-results$log2FC)])
ggplot(results, aes(x = cancer.type, y = log2FC, fill=significant)) +
  geom_hline(yintercept=1, color="lightgrey", linetype="dashed")+
  geom_hline(yintercept=-1, color="lightgrey", linetype="dashed")+
  geom_hline(yintercept=0, color="black")+
  geom_segment(aes(x = cancer.type, xend = cancer.type, y = 0, yend = log2FC)) +
  geom_point(shape = 21, size = 3) +
  labs(
    y = expression(Log[2]*FC ~ "(norm. count)"),
    title = "SERPINE1"
  ) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=c("lightgrey", "#7F61A6")) +
  theme(axis.text = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.line.x.bottom = element_line(linewidth = 0.4),
        axis.line.y.left = element_line(linewidth = 0.4),
        text = element_text(family = "Arial"),
        legend.text = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
setwd("../yourpath/")
ggsave(paste0("../yourpath/lfc.SERPINE1.png"), plot = last_plot(), dpi = 500, width = 7.5, height = 3)
ggsave(paste0("../yourpath/lfc.SERPINE1.pdf"), plot = last_plot(), dpi = 500, width = 7.5, height = 3)


#-------------------SERPINB2-------------------


df_TCGA_sub <- df_merged %>%
  filter(Gene.Name == "SERPINB2") %>% 
  filter(X_study == "TCGA") %>%
  filter(!X_sample_type == "Solid Tissue Normal") %>%
  select(sample,RSEM_norm_count, Cancer.Type)
df_TCGA_sub$Cancer.Type <- as.factor(df_TCGA_sub$Cancer.Type)

df_normal_sub <- df_merged %>%
  filter(Gene.Name == "SERPINB2") %>% 
  filter(X_sample_type %in% c("Solid Tissue Normal", "Normal Tissue", 'Cell Line')) %>%
  select(sample,RSEM_norm_count, X_primary_site)

# Mapping
mapping_df <- data.frame(
  Cancer.Type = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", 
                  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", 
                  "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", 
                  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS"),
  X_primary_site = c("Adrenal Gland", "Bladder", "Breast", "Cervix Uteri", "Bile duct", "Colon", 
                     "Blood", "Esophagus", "Brain", "Head and Neck region", "Kidney", "Kidney",  "Kidney", "Bone Marrow", "Brain", "Liver",
                     "Lung", "Lung",  "Ovary", "Pancreas", "Paraganglia", "Prostate", "Colon", 
                     "Soft tissue,Bone", "Skin", "Stomach", "Testis", "Thyroid", "Thymus", "Uterus", "Uterus"),
  stringsAsFactors = FALSE
)

results <- data.frame(cancer.type = character(), 
                      normal.tissue = character(),
                      cancer.RSEM.norm.count.mean = numeric(),
                      normal.RSEM.norm.count.mean = numeric(),
                      log2FC = numeric(),
                      p_value = numeric(), 
                      stringsAsFactors = FALSE)

for (cancer_type in unique(df_TCGA_sub$Cancer.Type)) {
  print(paste("Processing cancer type:", cancer_type))
  tissue_type <- mapping_df %>% filter(Cancer.Type == cancer_type) %>% pull(X_primary_site)
  if (length(tissue_type) == 0) {
    print(paste("No matching tissue type found for cancer type:", cancer_type))
    next
  }
  cancer_RSEM_norm_count <- df_TCGA_sub %>% filter(Cancer.Type == cancer_type) %>% pull(RSEM_norm_count)
  healthy_RSEM_norm_count <- df_normal_sub %>% filter(X_primary_site == tissue_type) %>% pull(RSEM_norm_count)
  
  if (length(healthy_RSEM_norm_count) == 0) {
    print(paste("No matching normal data for tissue type:", tissue_type))
    next
  }
  
  mean_cancer_RSEM_norm_count <- mean(cancer_RSEM_norm_count, na.rm = TRUE)
  mean_healthy_RSEM_norm_count <- mean(healthy_RSEM_norm_count, na.rm = TRUE)
  log2FC <- mean_cancer_RSEM_norm_count - mean_healthy_RSEM_norm_count
  wilcox_test <- wilcox.test(cancer_RSEM_norm_count, healthy_RSEM_norm_count, 
                             alternative = "two.sided", exact = FALSE)
  p_val <- wilcox_test$p.value
  
  results <- rbind(results, data.frame(
    cancer.type = cancer_type,
    normal.tissue = tissue_type,
    cancer.RSEM.norm.count.mean = mean_cancer_RSEM_norm_count,
    normal.RSEM.norm.count.mean = mean_healthy_RSEM_norm_count,
    log2FC = log2FC,
    p_value = p_val
  ))
}
results$p_adj <- p.adjust(results$p_value, method = "BH")  
print(results)

results$significant <- ifelse(results$p_adj < 0.05, "significant", "not significant")
results <- results[order(-results$log2FC), , drop = FALSE]
results$cancer.type <- factor(results$cancer.type, levels = results$cancer.type[order(-results$log2FC)])
ggplot(results, aes(x = cancer.type, y = log2FC, fill=significant)) +
  geom_hline(yintercept=1, color="lightgrey", linetype="dashed")+
  geom_hline(yintercept=-1, color="lightgrey", linetype="dashed")+
  geom_hline(yintercept=0, color="black")+
  geom_segment(aes(x = cancer.type, xend = cancer.type, y = 0, yend = log2FC)) +
  geom_point(shape = 21, size = 3) +
  labs(
    y = expression(Log[2]*FC ~ "(norm. count)"),
    title = "SERPINB2"
  ) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=c("lightgrey", "#C87D79")) +
  theme(axis.text = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.line.x.bottom = element_line(linewidth = 0.4),
        axis.line.y.left = element_line(linewidth = 0.4),
        text = element_text(family = "Arial"),
        legend.text = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
setwd("../yourpath/")
ggsave(paste0("../yourpath/lfc.SERPINB2.png"), plot = last_plot(), dpi = 500, width = 7.5, height = 3)
ggsave(paste0("../yourpath/lfc.SERPINB2.pdf"), plot = last_plot(), dpi = 500, width = 7.5, height = 3)

