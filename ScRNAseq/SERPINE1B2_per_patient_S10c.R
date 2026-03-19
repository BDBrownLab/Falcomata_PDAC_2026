library('tidyverse')
library('readxl')
library("writexl")
library("DESeq2")
library("dplyr")
theme_set(theme_minimal(base_family = "Arial"))

library(extrafont)
font_import()
loadfonts(device = "pdf")

setwd("../yourpath/")
barplot_data <- read_csv("../yourpath.csv")

barplot_data_long <- barplot_data %>%
  filter(n_cancer_cells_in_patient_id >= 100) %>%
  select(patient_id, n_SERPINE1_pos, n_SERPINB2_pos, n_SERPINE1_and_SERPINB2_pos) %>%
  pivot_longer(cols = c("n_SERPINE1_pos","n_SERPINB2_pos","n_SERPINE1_and_SERPINB2_pos"),
               names_to = c("serpin.status"),
               values_to = "value",
               values_drop_na = TRUE)

barplot_data_long$serpin.status <- factor(barplot_data_long$serpin.status, levels = c("n_SERPINE1_pos","n_SERPINB2_pos","n_SERPINE1_and_SERPINB2_pos"))

barplot_data$pct_double_pos_of_pos <- barplot_data$n_SERPINE1_and_SERPINB2_pos / (barplot_data$n_SERPINE1_and_SERPINB2_pos + barplot_data$n_SERPINE1_pos + barplot_data$n_SERPINB2_pos)

patient_order <- barplot_data %>%
  filter(n_cancer_cells_in_patient_id >= 100) %>%
  select(patient_id, pct_double_pos_of_pos) %>%
  arrange(desc(pct_double_pos_of_pos)) %>%
  pull(patient_id)

barplot_data_long2 <- barplot_data_long %>%
  mutate(patient_id = factor(patient_id, levels = patient_order))

ggplot(barplot_data_long2, aes(fill=serpin.status, y=value, x=patient_id)) + 
  geom_bar(position="fill", stat="identity", colour = NA, linewidth = 0,width = 1) +
  theme_classic(base_size = 15) +
  labs(x = "Proportion", 
       y = "") +
  scale_fill_manual(values = c("#7F61A6","#C87D79", "#E39D3F")) +
  theme(
    axis.text = element_text(size = 17, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.line.x.bottom = element_line(linewidth = 0.4),
    axis.line.y.left = element_line(linewidth = 0.4),
    text = element_text(family = "Arial"),
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  )
setwd("../yourpath/")
ggsave(paste0("../yourpath/human.scRNAseq.barplot.booleanE1B2.patient.png"), plot = last_plot(), dpi = 500, width = 7, height = 6)
ggsave(paste0("../yourpath/scRNAseq.barplot.booleanE1B2.patient.pdf"), plot = last_plot(), dpi = 500, width = 7, height = 6)
