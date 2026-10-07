library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)


##### GSE121634 Heymach Erlotinib resistant bulk RNAseq analysis #####

###HCC827 ER

#Read in data
heymach_hcc827_er_krt17 <- Heymach_HCC827_erlotinib_resistant_bulk_rnaseq_KRT17_expression

#Set group order
heymach_hcc827_er_krt17$Group <- factor(x = heymach_hcc827_er_krt17$Group, levels = c("HCC827_Parental", "HCC827_ER1", "HCC827_ER3", "HCC827_ER6"))

#Draw faceted boxplot for KRT17 expression
heymach_hcc827_er_krt17_bp <- ggplot(heymach_hcc827_er_krt17, aes(x = Group, y = KRT17_3872, fill = Group)) +
                       geom_boxplot() + 
                       theme_bw() +
                       theme(panel.grid = element_blank(),
                             axis.text = element_text(size = 12, color = "black"),
                             axis.title.x = element_blank(),
                             axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                       scale_y_continuous("KRT17 expression (TPM)") +
                       scale_fill_manual(values = c("HCC827_Parental" = "lightblue",
                                                    "HCC827_ER1" = "blue",
                                                    "HCC827_ER3" = "blue",
                                                    "HCC827_ER6" = "blue")) +
                       stat_compare_means(method = "anova", label.y = 3.5) +
                       labs(title = "HCC827 Parental vs Erolitinib Resistant KRT17 expression")

heymach_hcc827_er_krt17_bp #Plot: 600 x 400



###HCC4006 ER

#Read in data
heymach_hcc4006_er_krt17 <- Heymach_HCC4006_erlotinib_resistant_bulk_rnaseq_KRT17_expression

#Set group order
heymach_hcc4006_er_krt17$Group <- factor(x = heymach_hcc4006_er_krt17$Group, levels = c("HCC4006_Parental", "HCC4006_ER1", "HCC4006_ER2", "HCC4006_ER3", "HCC4006_ER5", "HCC4006_ER6"))

#Draw faceted boxplot for KRT17 expression
heymach_hcc4006_er_krt17_bp <- ggplot(heymach_hcc4006_er_krt17, aes(x = Group, y = KRT17_3872, fill = Group)) +
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous("KRT17 expression (TPM)") +
  scale_fill_manual(values = c("HCC4006_Parental" = "lightgreen",
                               "HCC4006_ER1" = "green",
                               "HCC4006_ER2" = "green",
                               "HCC4006_ER3" = "green",
                               "HCC4006_ER5" = "green",
                               "HCC4006_ER6" = "green")) +
  stat_compare_means(method = "anova", label.y = 150) +
  labs(title = "HCC4006 Parental vs Erolitinib Resistant KRT17 expression")

heymach_hcc4006_er_krt17_bp #Plot: 600 x 400



##### GSE253742 Treatment naive vs osi MRD bulk RNAseq analysis (4/30/26) #####

#Read in data
pre_post_osi_krt17 <- Chinese_patient_pre_post_osi_bulk_RNAseq_GSE253742_KRT17_expression

#Set group order
pre_post_osi_krt17$Group <- factor(x = pre_post_osi_krt17$Group, levels = c("Treatment_naive", "Osi_MRD"))

#Draw faceted boxplot for KRT17 expression
pre_post_osi_krt17_bp <- ggplot(pre_post_osi_krt17, aes(x = Group, y = KRT17_3872, fill = Group)) +
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous("KRT17 expression (TPM)", limits = c(0, 500)) +
  scale_fill_manual(values = c("Treatment_naive" = "indianred1",
                               "Osi_MRD" = "green3")) +
  stat_compare_means(method = "wilcox.test", label.y = 480, label.x = 0.9) +
  labs(title = "GSE253742 Treatment naive vs \nOsi MRD KRT17 expression")

pre_post_osi_krt17_bp #Plot: 350 x 450



##### GSM5777996 Astrazeneca osi persister time course RNAseq analysis #####

###HCC827

#Read in data
az_hcc827_krt17 <- AstraZeneca_HCC827_osi_persister_time_course_bulk_RNAseq_GSM5777996

#Set group order
az_hcc827_krt17$Group <- factor(x = az_hcc827_krt17$Group, levels = c("HCC827_DMSO", "HCC827_acute_osi", "HCC827_Osi_DTPC", "HCC827_short_washout", "HCC827_long_washout"))

#Draw faceted boxplot for KRT17 expression
az_hcc827_krt17_bp <- ggplot(az_hcc827_krt17, aes(x = Group, y = KRT17_3872, fill = Group)) +
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous("KRT17 expression (TPM)", limits = c(0, 250)) +
  scale_fill_manual(values = c("HCC827_DMSO" = "indianred1",
                               "HCC827_acute_osi" = "green3",
                               "HCC827_Osi_DTPC" = "blue",
                               "HCC827_short_washout" = "yellow",
                               "HCC827_long_washout" = "magenta")) +
  stat_compare_means(method = "anova", label.y = 220) +
  labs(title = "GSM5777996 HCC827 Osi persister time course \nKRT17 expression")

az_hcc827_krt17_bp #Plot: 800 x 450


###HCC2935

#Read in data
az_HCC2935_krt17 <- AstraZeneca_HCC2935_osi_persister_time_course_bulk_RNAseq_GSM5777996

#Set group order
az_HCC2935_krt17$Group <- factor(x = az_HCC2935_krt17$Group, levels = c("HCC2935_DMSO", "HCC2935_acute_osi", "HCC2935_Osi_DTPC", "HCC2935_short_washout", "HCC2935_long_washout"))

#Draw faceted boxplot for KRT17 expression
az_HCC2935_krt17_bp <- ggplot(az_HCC2935_krt17, aes(x = Group, y = KRT17_3872, fill = Group)) +
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous("KRT17 expression (TPM)", limits = c(0, 250)) +
  scale_fill_manual(values = c("HCC2935_DMSO" = "indianred1",
                               "HCC2935_acute_osi" = "green3",
                               "HCC2935_Osi_DTPC" = "blue",
                               "HCC2935_short_washout" = "yellow",
                               "HCC2935_long_washout" = "magenta")) +
  stat_compare_means(method = "anova", label.y = 220) +
  labs(title = "GSM5777996 HCC2935 Osi persister time course \nKRT17 expression")

az_HCC2935_krt17_bp #Plot: 800 x 500



###H1975

#Read in data
az_H1975_krt17 <- AstraZeneca_H1975_osi_persister_time_course_bulk_RNAseq_GSM5777996

#Set group order
az_H1975_krt17$Group <- factor(x = az_H1975_krt17$Group, levels = c("H1975_DMSO", "H1975_acute_osi", "H1975_Osi_DTPC", "H1975_short_washout", "H1975_long_washout"))

#Draw faceted boxplot for KRT17 expression
az_H1975_krt17_bp <- ggplot(az_H1975_krt17, aes(x = Group, y = KRT17_3872, fill = Group)) +
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous("KRT17 expression (TPM)", limits = c(0, 500)) +
  scale_fill_manual(values = c("H1975_DMSO" = "indianred1",
                               "H1975_acute_osi" = "green3",
                               "H1975_Osi_DTPC" = "blue",
                               "H1975_short_washout" = "yellow",
                               "H1975_long_washout" = "magenta")) +
  stat_compare_means(method = "anova", label.y = 475) +
  labs(title = "GSM5777996 H1975 Osi persister time course \nKRT17 expression")

az_H1975_krt17_bp #Plot: 800 x 500



###PC9

#Read in data
az_PC9_krt17 <- AstraZeneca_PC9_osi_persister_time_course_bulk_RNAseq_GSM5777996

#Set group order
az_PC9_krt17$Group <- factor(x = az_PC9_krt17$Group, levels = c("PC9_DMSO", "PC9_acute_osi", "PC9_osi_DTPC", "PC9_short_washout", "PC9_long_washout"))

#Draw faceted boxplot for KRT17 expression
az_PC9_krt17_bp <- ggplot(az_PC9_krt17, aes(x = Group, y = KRT17_3872, fill = Group)) +
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous("KRT17 expression (TPM)", limits = c(0, 15)) +
  scale_fill_manual(values = c("PC9_DMSO" = "indianred1",
                               "PC9_acute_osi" = "green3",
                               "PC9_osi_DTPC" = "blue",
                               "PC9_short_washout" = "yellow",
                               "PC9_long_washout" = "magenta")) +
  stat_compare_means(method = "anova", label.y = 14) +
  labs(title = "GSM5777996 PC9 Osi persister time course \nKRT17 expression")

az_PC9_krt17_bp #Plot: 800 x 500












