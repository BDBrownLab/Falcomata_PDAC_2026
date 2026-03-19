library('tidyverse')
library('readxl')
library("writexl")
library("dplyr")
theme_set(theme_minimal(base_family = "Arial"))

library(extrafont)
font_import()
loadfonts(device = "pdf")

#Barplot SERPINE1+/-SERPINB2+/- subtypes
setwd("../yourpath/")
barplot_data <- read_csv("../yourpath.csv")

barplot_data_long <- barplot_data %>%
  pivot_longer(cols = c("SERPINE1+SERPINB2+","SERPINE1-SERPINB2+","SERPINE1+SERPINB2-","SERPINE1-SERPINB2-"),
               names_to = c("serpin.status"),
               values_to = "value",
               values_drop_na = TRUE)

barplot_data_long$anno_sub <- factor(barplot_data_long$anno_sub)
levels(barplot_data_long$anno_sub) <- c("basal-like", "classical", "exocrine-like")
barplot_data_long$serpin.status <- factor(barplot_data_long$serpin.status, levels = c("SERPINE1+SERPINB2+", "SERPINE1+SERPINB2-", "SERPINE1-SERPINB2+","SERPINE1-SERPINB2-"))

ggplot(barplot_data_long, aes(fill=serpin.status, y=value, x=source)) + 
  geom_bar(position="fill", stat="identity", color = "black") +
  facet_wrap(~anno_sub) +
  theme_classic(base_size = 15) +
  labs(x = "Proportion", 
       y = "") +
  scale_fill_manual(values = c("#E39D3F","#7F61A6","#C87D79", 'lightgrey')) +
  theme(
    axis.text = element_text(size = 17, color = "black"),
    axis.text.x = element_text(angle = 90,size = 14, hjust = 1,vjust = 0.5),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.line.x.bottom = element_line(linewidth = 0.4),
    axis.line.y.left = element_line(linewidth = 0.4),
    text = element_text(family = "Arial"),
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  )
setwd("../yourpath/")
ggsave(paste0("../yourpath/human.scRNAseq.barplot.booleanE1B2.wNEG.study.png"), plot = last_plot(), dpi = 500, width = 8, height = 6)
ggsave(paste0("../yourpath/human.scRNAseq.barplot.booleanE1B2.wNEG.study.pdf"), plot = last_plot(), dpi = 500, width = 8, height = 6)

barplot_data_long %>%
  filter(serpin.status != "SERPINE1-SERPINB2-") -> barplot_data_long
ggplot(barplot_data_long, aes(fill=serpin.status, y=value, x=source)) + 
  geom_bar(position="fill", stat="identity", color = "black") +
  facet_wrap(~anno_sub) +
  theme_classic(base_size = 15) +
  labs(x = "Proportion", 
       y = "") +
  scale_fill_manual(values = c("#E39D3F","#7F61A6","#C87D79", 'lightgrey')) +
  theme(
    axis.text = element_text(size = 17, color = "black"),
    axis.text.x = element_text(angle = 90,size = 14, hjust = 1,vjust = 0.5),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.line.x.bottom = element_line(linewidth = 0.4),
    axis.line.y.left = element_line(linewidth = 0.4),
    text = element_text(family = "Arial"),
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  )
setwd("../yourpath/")
ggsave(paste0("../yourpath/human.scRNAseq.barplot.booleanE1B2.study.png"), plot = last_plot(), dpi = 500, width = 8, height = 6)
ggsave(paste0("../yourpath/human.scRNAseq.barplot.booleanE1B2.study.pdf"), plot = last_plot(), dpi = 500, width = 8, height = 6)
