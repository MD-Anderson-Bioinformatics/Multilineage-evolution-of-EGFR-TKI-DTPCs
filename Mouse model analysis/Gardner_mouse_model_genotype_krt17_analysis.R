library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)



##### Gardner mouse model genotype Krt17 expression analysis #####


#Read in data
gardner_geno_krt17 <- Gardner_mouse_model_KRT17_expression_by_genotype_and_condition

#Set sample order
gardner_geno_krt17$Sample <- factor(gardner_geno_krt17$Sample, levels = c("AT2_1mo", "Spc_PtenMT", "EPMT", "EPMT_MRD", "ERPT", "ERPMT_pool", "ERPMT_osiMRD", "ERPMT_doxMRD", "ERPMT_3moMRD", "ERPMT_MRDrestart"))

#Plot barchart
gardner_geno_krt17_bp <- ggplot(gardner_geno_krt17, aes(x = Sample, y = Krt17_pos_percentage, fill = Sample)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = signif(Krt17_pos_percentage, 3), vjust = -0.3)) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 12, color = "black"),
        axis.title.x = element_blank()) +
  scale_y_continuous("Percentage Krt17+", expand = c(0, 0), limits = c(0, 10)) +
  labs(title = "Gardner mouse genotype Krt17 expression")

gardner_geno_krt17_bp #Plot: 800 x 450








