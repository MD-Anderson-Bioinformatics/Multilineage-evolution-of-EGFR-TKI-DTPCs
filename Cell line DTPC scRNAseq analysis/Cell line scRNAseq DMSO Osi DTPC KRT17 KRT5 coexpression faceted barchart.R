library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(forcats)



##### HCC827 HCC4006 H1975 DMSO Osi DTPC KRT17/KRT5 co-expression faceted barchart #####

#Read in data
dmso_dtpc_krt17_krt5_f <- HCC827_HCC4006_H1975_DMSO_Osi_DTPC_KRT17_KRT5_co_expression_faceted_barchart_data

#Set plot order
dmso_dtpc_krt17_krt5_f$Cell_line <- factor(dmso_dtpc_krt17_krt5_f$Cell_line, levels = c("HCC827", "HCC4006", "H1975"))

#Draw barchart
dmso_dtpc_krt17_krt5_f_bar <- ggplot(dmso_dtpc_krt17_krt5_f, aes(x = dmso_dtpc_krt17_krt5_f$Condition, y = dmso_dtpc_krt17_krt5_f$Percentage, fill = dmso_dtpc_krt17_krt5_f$Type)) +
  facet_wrap(~ Cell_line) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right") +
  scale_fill_manual(values = c("KRT17+ KRT5-" = "red",
                               "KRT17+ KRT5+" = "green")) +
  scale_y_continuous("Fraction") +
  labs(title = "KRT17+ cell KRT5 co-expression", fill = "Co-expression")


dmso_dtpc_krt17_krt5_f_bar #Plot 550 x 450




