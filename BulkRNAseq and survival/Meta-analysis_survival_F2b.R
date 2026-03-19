library(ggplot2)
library(readxl)
library(writexl)
library(dplyr)
library(broom)
library(survival)
library(survminer)
library(tidyr)
library(purrr)
library(metafor)
library(meta)
library(readr)

#Load data
setwd('../yourpath/')
df_merged <- read.csv('../yourpath/250625.XenaHubs.TcgaTargetGtex.normalized.metadata.csv')
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

df_merged$daysToEvent <- df_merged$OS.time
df_merged$event <- df_merged$OS
TCGA_PAAD_sub <- df_merged %>%
  filter(Cancer.Type == "PAAD") %>%
  dplyr::select(sample, SERPINE1, SERPINB2, SERPINF2,event, daysToEvent)
TCGA_PAAD_sub$study <- "TCGA-PAAD"

setwd('../yourpath/')

#Moffitt2015
#death=1 censor=0
#time to event in days
Moffitt2015_data <- read_csv('../yourpath/GSE71729_Moffitt2015_gene_by_sample.csv')
Moffitt2015_md <- read_csv('../yourpath/Moffitt2015_GSE71729_sample_annotations.csv')
Moffitt2015_data$sample <- sapply(strsplit(Moffitt2015_data$sample, "_"), `[`, 1)
Moffitt2015_md$sample <- Moffitt2015_md$geo_accession
Moffitt2015 <- merge(Moffitt2015_data, Moffitt2015_md, by="sample", how="inner")
Moffitt2015$month <- as.numeric(Moffitt2015$survival_months)
Moffitt2015$daysToEvent <- Moffitt2015$month * 30.45
Moffitt2015$event <- Moffitt2015$death_event_1death_0censor
Moffitt2015_sub <- as_tibble(Moffitt2015) %>%
  dplyr::select(sample, SERPINE1, SERPINF2, SERPINB2, daysToEvent, event)
Moffitt2015_sub$study <- "Moffitt2015"

#Puleo2018
#death=1 censor=0
#time to event in days
Puleo2018_data <- read_csv('../yourpath/E-MTAB-6134_Puleo2018_gene_by_sample.csv')
Puleo2018_md <- read_tsv('../yourpath/Puleo2018/E-MTAB-6134.sdrf.txt')
Puleo2018_md$sample <- Puleo2018_md$`Source Name`
Puleo2018 <- merge(Puleo2018_data, Puleo2018_md, by="sample", how="inner")
Puleo2018$month <- as.numeric(Puleo2018$`Characteristics[os.delay]`)
Puleo2018$daysToEvent <- Puleo2018$month * 30.45
Puleo2018$event <- Puleo2018$`Characteristics[os.event]`
Puleo2018_sub <- as_tibble(Puleo2018) %>%
  dplyr::select(sample, SERPINE1, SERPINF2, SERPINB2, daysToEvent, event)
Puleo2018_sub$study <- "Puleo2018"

#Yang2016
#death=1 censor=0
#time to event in days
Yang2016_data <- read_csv('../yourpath/GSE62452_Yang2016_gene_by_sample.csv')
Yang2016_md <- read_csv('../yourpath/Yang2016_GSE62452_sample_annotations.csv')
Yang2016_data$sample <- sapply(strsplit(Yang2016_data$sample, "_"), `[`, 1)
Yang2016_md$sample <- Yang2016_md$geo_accession
Yang2016 <- merge(Yang2016_data, Yang2016_md, by="sample", how="inner")
Yang2016$month <- as.numeric(Yang2016$survival_months)
Yang2016$daysToEvent <- Yang2016$month * 30.45
Yang2016$event <- Yang2016$survival_status
Yang2016_sub <- as_tibble(Yang2016) %>%
  dplyr::select(sample, SERPINE1, SERPINF2, SERPINB2, daysToEvent, event)
Yang2016_sub$study <- "Yang2016"

#Ohara2023
#death=1 censor=0
#time to event in days
Yuuki2023_data <- read_csv('../yourpath/GSE224564_Yuuki2023_gene_by_sample.csv')
Yuuki2023_md <- read_csv('../yourpath/Yuuki2023_GSE223909_sample_annotations.csv')
Yuuki2023_md$sample <- Yuuki2023_md$geo_accession
Yuuki2023 <- merge(Yuuki2023_data, Yuuki2023_md, by="sample", how="inner")
Yuuki2023$month <- as.numeric(Yuuki2023$month)
Yuuki2023$daysToEvent <- Yuuki2023$month * 30.45
Yuuki2023$event <- Yuuki2023$cancer_death
Yuuki2023_sub <- as_tibble(Yuuki2023) %>%
  dplyr::select(sample, SERPINE1, SERPINF2, SERPINB2, daysToEvent, event)
Yuuki2023_sub$study <- "Ohara2023"

#Zhang2012
#death=1 censor=0
#time to event in days
Zhang2012_data <- read_csv('../yourpath/GSE28735_Zhang2012_gene_by_sample.csv')
Zhang2012_md <- read_csv('../yourpath/Zhang2012_GSE28735_sample_annotations.csv')
Zhang2012_md$sample <- Zhang2012_md$geo_accession
Zhang2012 <- merge(Zhang2012_data, Zhang2012_md, by="sample", how="inner")
Zhang2012$month <- as.numeric(Zhang2012$survival_month)
Zhang2012$daysToEvent <- Zhang2012$month * 30.45
Zhang2012$event <- Zhang2012$cancer_death
Zhang2012_sub <- as_tibble(Zhang2012) %>%
  dplyr::select(sample, SERPINE1, SERPINF2, SERPINB2, daysToEvent, event)
Zhang2012_sub$study <- "Zhang2012"

#ICGC-PACA-AU
#death=1, censor=0
#time to event in days
ICGCPACAAU_data <- read_csv('../yourpath/ICGC-PACA-AU_gene_by_sample.csv')
ICGCPACAAU_md <- read_csv('../yourpath/ICGC-PACA-AU_sample_annotations.csv')
ICGCPACAAU <- merge(ICGCPACAAU_data, ICGCPACAAU_md, by="sample", how="inner")
ICGCPACAAU <- ICGCPACAAU %>%
  filter(Sample.type == "Primary tumour")
ICGCPACAAU$daysToEvent <- ICGCPACAAU$survivalA
ICGCPACAAU$event <- ICGCPACAAU$censorA.0yes.1no
ICGCPACAAU_sub <- as_tibble(ICGCPACAAU) %>%
  dplyr::select(sample, SERPINE1, SERPINF2, SERPINB2, daysToEvent, event)
ICGCPACAAU_sub$study <- "ICGC-PACA-AU"

#ICGC-PACA-AU
#death=1, censor=0
#time to event in days
ICGCPACACA_data <- read_csv('../yourpath/ICGC-PACA-CA_gene_by_sample.csv')
ICGCPACACA_md <- read_csv('../yourpath/ICGC-PACA-CA_sample_annotations.csv')
ICGCPACACA <- merge(ICGCPACACA_data, ICGCPACACA_md, by="sample", how="inner")
ICGCPACACA <- ICGCPACACA %>%
  filter(specimen_type %in% c("Primary tumour - solid tissue", "Primary tumour - other"))
ICGCPACACA$daysToEvent <- ICGCPACACA$survivalA
ICGCPACACA$event <- ICGCPACACA$censorA.0yes.1no
ICGCPACACA_sub <- as_tibble(ICGCPACACA) %>%
  dplyr::select(sample, SERPINE1, SERPINF2, SERPINB2, daysToEvent, event)
ICGCPACACA_sub$study <- "ICGC-PACA-CA"

data <- rbind(Moffitt2015_sub, Puleo2018_sub, Yang2016_sub,Yuuki2023_sub,Zhang2012_sub,TCGA_PAAD_sub,ICGCPACACA_sub,ICGCPACAAU_sub)
data <- na.omit(data)

################### SERPINE1 quartiles----
dat2 <- data %>%
  mutate(
    event = as.integer(trimws(as.character(event)))
  ) %>%
  filter(
    !is.na(daysToEvent),
    !is.na(event),
    !is.na(SERPINE1),
    !is.na(study)
  ) %>%
  group_by(study) %>%
  mutate(
    q25 = quantile(SERPINE1, 0.25, na.rm = TRUE),
    q75 = quantile(SERPINE1, 0.75, na.rm = TRUE),
    serpine1_group = case_when(
      SERPINE1 <= q25 ~ "low",
      SERPINE1 >= q75 ~ "high",
      TRUE ~ NA_character_
    ),
    serpine1_group = factor(serpine1_group, levels = c("low","high"))
  ) %>%
  ungroup() %>%
  filter(serpine1_group %in% c("high","low"))

eligible <- dat2 %>%
  filter(!is.na(serpine1_group)) %>%
  group_by(study, serpine1_group) %>%
  summarise(n = n(), events = sum(event == 1), .groups = "drop_last") %>%
  pivot_wider(
    names_from = serpine1_group,
    values_from = c(n, events),
    values_fill = 0
  ) %>%
  filter(events_low > 0 & events_high > 0) %>%
  distinct(study) %>%
  pull(study)

dat3 <- dat2 %>% filter(study %in% eligible)

per_study <- dat3 %>%
  group_by(study) %>%
  group_modify(~{
    fit <- survival::coxph(Surv(daysToEvent, event) ~ serpine1_group, data = .x)
    s <- summary(fit)
    ci <- confint(fit)  
    tibble(
      n = nrow(.x),
      events = sum(.x$event == 1),
      hr = unname(exp(coef(fit)[1])),
      logHR = unname(coef(fit)[1]),
      se = unname(s$coef[1, "se(coef)"]),
      ci_l = unname(exp(ci[1])),
      ci_u = unname(exp(ci[2])),
      p = unname(s$coef[1, "Pr(>|z|)"])
    )
  }) %>%
  ungroup() %>%
  arrange(hr)

m_serpine1 <- metagen(
  TE= logHR,seTE= se,studlab= study,data= per_study,sm= "HR",
  comb.fixed = FALSE,comb.random= TRUE,method.tau = "REML",hakn= TRUE)

summary(m_serpine1)

pdf("../yourpath/forestplot.SERPINE1.quartiles.pdf", width = 7, height = 3.5)

forest(m_serpine1,
       comb.fixed=F,comb.random=T,prediction=F,backtransf=T,leftcols = c("studlab"),rightcols=c("effect","ci"),print.pval=F,
       print.tau2 = F, print.I2 = F, print.pval.Q = F,col.square = "#9584C1",col.square.lines = "#9584C1",
       col.diamond = "#71489A",col.diamond.lines = "black",)

dev.off()

################### SERPINB2 quartiles----
dat2 <- data %>%
  mutate(
    event = as.integer(trimws(as.character(event)))
  ) %>%
  filter(
    !is.na(daysToEvent),
    !is.na(event),
    !is.na(SERPINB2),
    !is.na(study)
  ) %>%
  group_by(study) %>%
  mutate(
    q25 = quantile(SERPINB2, 0.25, na.rm = TRUE),
    q75 = quantile(SERPINB2, 0.75, na.rm = TRUE),
    serpinb2_group = case_when(
      SERPINB2 <= q25 ~ "low",
      SERPINB2 >= q75 ~ "high",
      TRUE ~ NA_character_
    ),
    serpinb2_group = factor(serpinb2_group, levels = c("low","high"))
  ) %>%
  ungroup() %>%
  filter(serpinb2_group %in% c("high","low"))

eligible <- dat2 %>%
  filter(!is.na(serpinb2_group)) %>%
  group_by(study, serpinb2_group) %>%
  summarise(n = n(), events = sum(event == 1), .groups = "drop_last") %>%
  pivot_wider(
    names_from = serpinb2_group,
    values_from = c(n, events),
    values_fill = 0
  ) %>%
  filter(events_low > 0 & events_high > 0) %>%
  distinct(study) %>%
  pull(study)

dat3 <- dat2 %>% filter(study %in% eligible)

per_study <- dat3 %>%
  group_by(study) %>%
  group_modify(~{
    fit <- survival::coxph(Surv(daysToEvent, event) ~ serpinb2_group, data = .x)
    s   <- summary(fit)
    ci  <- confint(fit)  
    tibble(
      n = nrow(.x),
      events = sum(.x$event == 1),
      hr = unname(exp(coef(fit)[1])),
      logHR = unname(coef(fit)[1]),
      se = unname(s$coef[1, "se(coef)"]),
      ci_l = unname(exp(ci[1])),
      ci_u= unname(exp(ci[2])),
      p = unname(s$coef[1, "Pr(>|z|)"])
    )
  }) %>%
  ungroup() %>%
  arrange(hr)

m_serpine1 <- metagen(
  TE= logHR,seTE= se,studlab= study,data= per_study,sm= "HR",
  comb.fixed = FALSE,comb.random= TRUE,method.tau = "REML",hakn= TRUE)

summary(m_serpine1)

pdf("../yourpath/forestplot.SERPINB2.quartiles.pdf", width = 7, height = 3.5)

forest(m_serpine1,
       comb.fixed=F,comb.random=T,prediction=F,backtransf=T,leftcols = c("studlab"),rightcols=c("effect","ci"),print.pval=F,
       print.tau2 = F, print.I2 = F, print.pval.Q = F,col.square = "#D07371",col.square.lines = "#D07371",
       col.diamond = "#CB4665",col.diamond.lines = "black",)

dev.off()
