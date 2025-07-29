library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(SeuratExtend)
library(SeuratDisk)
library(sctransform)
library(scplotter)
library(BPCells)
library(presto)
library(glmGamPoi)
library(scran)
library(ggpubr)
library(harmony)
library(DoubletFinder)
library(kableExtra)
library(CytoTRACE2)
library(HiClimR)
library(devtools)
library(forcats)
library(scCustomize)
library(scales)
library(patchwork)


##### Initial clinical cohort treatment naive KRT17 analysis #####

#Load in data
egfr_combined <- LoadSeuratRds("EGFR_data.Rds")

#List of treatment naive tumor samples
specimen_ids_to_subset <- c("Biopsy-1", "Lung-Tumor-10", "JH064",
                            "JH067", "JH067_ATAC", "JH095",
                            "JH128", "JH139", "JH145", "JH-217",
                            "JH-295", "JH-308", "JH-310", "JH114")


#Subset for TN samples only
egfr_epi_tn <- subset(egfr_combined, subset = Specimen_ID %in% specimen_ids_to_subset)

#Filter for epithelial cells
egfr_epi_tn <- subset(egfr_epi_tn, subset = General_Annot == "EPITHELIAL")

#Filter for InferCnv annotated tumor cells
egfr_epi_t_tn <- subset(egfr_epi_tn, subset = InfercnvR == "tumor")


#Summarize the number of epithelial tumor cells by Specimen_ID
egfr_epi_stats <- egfr_epi_t_tn@meta.data %>%
  group_by(Specimen_ID) %>%
  summarise(Count = n())

write.table(egfr_epi_stats, "C:/Users/bbmorris1/Desktop/MDACC EGFR mutant scRNAseq treatment naive sample epithelial cell count by Specimen ID.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Filter for treatment naive tumor samples with >20 tumor cells
specimen_ids_to_subset <- c("Biopsy-1", "Lung-Tumor-10", "JH064",
                            "JH067", "JH095", "JH128", "JH139")

#Subset for epithelial TN > 20 cells only
egfr_epi_tn <- subset(egfr_combined, subset = Specimen_ID %in% specimen_ids_to_subset)

#Filter for epithelial cells
egfr_epi_tn <- subset(egfr_epi_tn, subset = General_Annot == "EPITHELIAL")

#Filter for InferCnv annotated tumor cells
egfr_epi_t_tn <- subset(egfr_epi_tn, subset = InfercnvR == "tumor")

#Change active identity to Specimen ID
Idents(egfr_epi_t_tn) <- "Specimen_ID"

#Draw KRT17 expression violin plot across treatment naive samples
krt17_tn_vp <- VlnPlot(egfr_epi_t_tn, features = "KRT17", group.by = "Specimen_ID") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("Biopsy-1" = "#F8766D",
                               "JH064" = "#E9842C",
                               "JH067" = '#D69100',
                               "JH095" = '#BC9D00',
                               "JH128" = '#9CA700',
                               "JH139" = '#6FB000',
                               "Lung-Tumor-10" = '#00BD61')) +
  labs(title = "Treatment Naive KRT17 expression")


krt17_tn_vp



###Annotate KRT17+ vs KRT17- cells by Specimen_ID

#Label KRT17+ cells
egfr_epi_t_tn$KRT17_exp <- ifelse(FetchData(egfr_epi_t_tn, vars = "KRT17") > 0, 
                                  "KRT17+", "KRT17-")

head(egfr_epi_t_tn)

#Summarize the number of cells types by Specimen_ID
egfr_epi_t_tn_krt17_stats <- table(egfr_epi_t_tn@meta.data$Specimen_ID, egfr_epi_t_tn@meta.data$KRT17_exp)
egfr_epi_t_tn_krt17_stats

#Write data to table
write.table(egfr_epi_t_tn_krt17_stats, "C:/Users/bbmorris1/Desktop/MDACC EGFR mutant scRNAseq treatment naive KRT17 expression counts.txt", sep = "\t", row.names = TRUE, quote = FALSE)


###Percent KRT17 expressing barchart

#Read in data
egfr_epi_t_tn_krt17_d <- MDACC_EGFR_mutant_scRNAseq_treatment_naive_KRT17_expression_counts_with_percentage

#Round to 2 decimal places for plotting
egfr_epi_t_tn_krt17_d$KRT17_perc_pos_round <- round(egfr_epi_t_tn_krt17_d$KRT17_perc_pos, 2)
egfr_epi_t_tn_krt17_d$KRT17_perc_neg_round <- round(egfr_epi_t_tn_krt17_d$KRT17_perc_neg, 2)

#Draw barchart summarizing % of KRT17+ cells in each biopsy
egfr_epi_t_tn_krt17_bar <- ggplot(egfr_epi_t_tn_krt17_d, aes(x = egfr_epi_t_tn_krt17_d$Specimen_ID, y = egfr_epi_t_tn_krt17_d$KRT17_perc_pos_round, fill = egfr_epi_t_tn_krt17_d$Specimen_ID)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = KRT17_perc_pos_round), vjust = -0.5) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("Biopsy-1" = "#F8766D",
                               "JH064" = "#E9842C",
                               "JH067" = '#D69100',
                               "JH095" = '#BC9D00',
                               "JH128" = '#9CA700',
                               "JH139" = '#6FB000',
                               "Lung-Tumor-10" = '#00BD61')) +
  scale_y_continuous("Percent KRT17+", limits = c(0, 100)) +
  labs(title = "Treatment Naive Percent KRT17 expressing")


egfr_epi_t_tn_krt17_bar


###KRT17/KRT5 co-expression analysis by Sample_ID

#Label KRT5+ cells
egfr_epi_t_tn$KRT5_exp <- ifelse(FetchData(egfr_epi_t_tn, vars = "KRT5") > 0, 
                                 "KRT5+", "KRT5-")

head(egfr_epi_t_tn)


###Count the number of KRT17+ and KRT5+ cells for all MRD samples
egfr_epi_t_tn_krt17_krt5_counts <- egfr_epi_t_tn@meta.data %>%
  mutate(sample_id = Specimen_ID) %>%  
  group_by(Specimen_ID, KRT17_exp, KRT5_exp) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  pivot_wider(names_from = c(KRT17_exp, KRT5_exp), 
              values_from = cell_count,
              values_fill = 0)

print(egfr_epi_t_tn_krt17_krt5_counts)

write.table(egfr_epi_t_tn_krt17_krt5_counts, "C:/Users/bbmorris1/Desktop/MDACC EGFR mutant treatment naive tumor cell KRT17 KRT5 co expression cell counts.txt", sep = '\t', row.names = FALSE, quote = FALSE)


###KRT17/KRT5 co-expression barchart

#Read in data
egfr_epi_t_tn_krt17_krt5 <- MDACC_EGFR_mutant_scRNAseq_treatment_naive_KRT17_KRT5_expression_with_percentage

#Draw barchart
egfr_epi_t_tn_krt17_krt5_bar <- ggplot(egfr_epi_t_tn_krt17_krt5, aes(x = egfr_epi_t_tn_krt17_krt5$Specimen_ID, y = egfr_epi_t_tn_krt17_krt5$KRT17_pos_KRT5_neg_perc)) +
  geom_bar(stat = "identity", fill = "red") +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Percentage", limits = c(0, 100)) +
  labs(title = "Treatment Naive KRT17+ KRT5+ co-expression")


egfr_epi_t_tn_krt17_krt5_bar


#Stitch plots together for visualization
krt17_tn_vp | egfr_epi_t_tn_krt17_bar | egfr_epi_t_tn_krt17_krt5_bar #1400 x 400



##### Initial clinical cohort Osi MRD KRT17 analysis #####

#Load in data
egfr_combined <- LoadSeuratRds("EGFR_data.Rds")

#List of Osi MRD tumor samples
specimen_ids_to_subset <- c("NSTAR1-TumorA", "NSTAR1-TumorB",
                            "Lung-Tumor-1", "Lung-Tumor-2",
                            "Lung-Tumor-8", "JH050", "JH053", "JH099",
                            "JH105", "JH105_ATAC", "JH109", "JH125",
                            "JH138", "JH143", "JH144", "JH-297",
                            "JH-316T", "Lung-Tumor-7")


#Subset for Osi MRD only
egfr_epi_omrd <- subset(egfr_combined, subset = Specimen_ID %in% specimen_ids_to_subset)

#Subset for epithelial only
egfr_epi_omrd <- subset(egfr_epi_omrd, subset = General_Annot == "EPITHELIAL")

#Subset for InferCnv tumor cells
egfr_epi_t_omrd <- subset(egfr_epi_omrd, subset = InfercnvR == "tumor")


#Summarize the number of tumor cells by Specimen_ID
egfr_epi_ormd_stats <- egfr_epi_t_omrd@meta.data %>%
  group_by(Specimen_ID) %>%
  summarise(Count = n())

write.table(egfr_epi_ormd_stats, "C:/Users/bbmorris1/Desktop/MDACC EGFR mutant scRNAseq Osi MRD sample epithelial tumor cell count by Specimen ID.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Filter for Osi MRD tumor samples with > 20 tumor cells (exclude JH105_ATAC due to NE component & possible NE progression)
specimen_ids_to_subset <- c("JH-297",
                            "JH143",
                            "Lung-Tumor-1",
                            "Lung-Tumor-7",
                            "Lung-Tumor-8",
                            "NSTAR1-TumorA")

#Subset for Osi MRD only
egfr_epi_omrd <- subset(egfr_combined, subset = Specimen_ID %in% specimen_ids_to_subset)

#Subset for epithelial only
egfr_epi_omrd <- subset(egfr_epi_omrd, subset = General_Annot == "EPITHELIAL")

#Subset for InferCnv tumor cells
egfr_epi_t_omrd <- subset(egfr_epi_omrd, subset = InfercnvR == "tumor")


#Change active identity to Specimen ID
Idents(egfr_epi_t_omrd) <- "Specimen_ID"

#Draw KRT17 expression violin plot across Osi MRD samples
krt17_epi_t_mrd_vp <- VlnPlot(egfr_epi_t_omrd, features = "KRT17", group.by = "Specimen_ID") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH-297" = "#00B813",
                               "JH143" = "#00C08E",
                               "Lung-Tumor-1" = "#00C0B4",
                               "Lung-Tumor-7" = "#00BDD4",
                               "Lung-Tumor-8" = "#00A7FF",
                               "NSTAR1-TumorA" = "#7F96FF")) +
  labs(title = "Osi MRD KRT17 expression")

krt17_epi_t_mrd_vp


###Annotate KRT17+ vs KRT17- cells by Specimen ID

#Label KRT17+ cells
egfr_epi_t_omrd$KRT17_exp <- ifelse(FetchData(egfr_epi_t_omrd, vars = "KRT17") > 0, 
                                    "KRT17+", "KRT17-")

head(egfr_epi_t_omrd)

#Summarize the number of cells types by Specimen_ID
egfr_epi_ormd_krt17_stats <- table(egfr_epi_t_omrd@meta.data$Specimen_ID, egfr_epi_t_omrd@meta.data$KRT17_exp)
egfr_epi_ormd_krt17_stats

#Write data to file
write.table(egfr_epi_ormd_krt17_stats, "C:/Users/bbmorris1/Desktop/MDACC EGFR mutant scRNAseq osimertinib MRD KRT17 expression stats.txt", sep = "\t", row.names = TRUE, quote = FALSE)


###Summary barchart for KRT17+%

#Read in data
egfr_epi_t_ormd_krt17_d <- MDACC_EGFR_mutant_scRNAseq_osimertinib_MRD_KRT17_expression_stats_with_percentage

#Round to 2 decimal places for plotting
egfr_epi_t_ormd_krt17_d$KRT17_perc_pos_round <- round(egfr_epi_t_ormd_krt17_d$KRT17_perc_pos, 2)
egfr_epi_t_ormd_krt17_d$KRT17_perc_neg_round <- round(egfr_epi_t_ormd_krt17_d$KRT17_perc_neg, 2)

#Draw barchart summarizing % of KRT17+ cells in each biopsy
egfr_epi_t_ormd_krt17_bar <- ggplot(egfr_epi_t_ormd_krt17_d, aes(x = egfr_epi_t_ormd_krt17_d$Specimen_ID, y = egfr_epi_t_ormd_krt17_d$KRT17_perc_pos_round, fill = egfr_epi_t_ormd_krt17_d$Specimen_ID)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = KRT17_perc_pos_round), vjust = -0.5) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("JH-297" = "#00B813",
                               "JH143" = "#00C08E",
                               "Lung-Tumor-1" = "#00C0B4",
                               "Lung-Tumor-7" = "#00BDD4",
                               "Lung-Tumor-8" = "#00A7FF",
                               "NSTAR1-TumorA" = "#7F96FF")) +
  scale_y_continuous("Percent KRT17+", limits = c(0, 100)) +
  labs(title = "Osimertinib MRD")


egfr_epi_t_ormd_krt17_bar


###KRT17/KRT5 co-expression analysis by Specimen_ID

#Label KRT5+ cells
egfr_epi_t_omrd$KRT5_exp <- ifelse(FetchData(egfr_epi_t_omrd, vars = "KRT5") > 0, 
                                   "KRT5+", "KRT5-")

head(egfr_epi_t_omrd)


###Count the number of KRT17+ and KRT5+ cells for all MRD samples
egfr_epi_t_omrd_krt17_krt5_counts <- egfr_epi_t_omrd@meta.data %>%
  mutate(sample_id = Specimen_ID) %>%  
  group_by(Specimen_ID, KRT17_exp, KRT5_exp) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  pivot_wider(names_from = c(KRT17_exp, KRT5_exp), 
              values_from = cell_count,
              values_fill = 0)

print(egfr_epi_t_omrd_krt17_krt5_counts)

write.table(egfr_epi_t_omrd_krt17_krt5_counts, "C:/Users/bbmorris1/Desktop/MDACC EGFR mut scRNAseq Osi MRD tumor cell KRT17 KRT5 co expression cell counts.txt", sep = '\t', row.names = FALSE)



###KRT17/KRT5 co-expression barchart

#Read in data
egfr_epi_t_ormd_krt17_krt5 <- MDACC_EGFR_mut_scRNAseq_Osi_MRD_KRT17_KRT5_co_expression_barchart_data

#Draw barchart
egfr_epi_t_tn_krt17_krt5_bar <- ggplot(egfr_epi_t_ormd_krt17_krt5, aes(x = egfr_epi_t_ormd_krt17_krt5$Specimen_ID, y = egfr_epi_t_ormd_krt17_krt5$Percentage, fill = egfr_epi_t_ormd_krt17_krt5$Type)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("KRT17+ KRT5-" = "red",
                               "KRT17+ KRT5+" = "green")) +
  scale_y_continuous("Fraction") +
  labs(title = "Osi MRD KRT17+ KRT5+ co-expression")


egfr_epi_t_tn_krt17_krt5_bar


#Stitch plots together for visualization
krt17_epi_t_mrd_vp | egfr_epi_t_ormd_krt17_bar | egfr_epi_t_tn_krt17_krt5_bar #1400 x 400



##### Initial clinical cohort Osi progression (PD) KRT17 analysis #####

#Load in data
egfr_combined <- LoadSeuratRds("EGFR_data.Rds")

#List of Osi progression tumor samples
specimen_ids_to_subset <- c("Biopsy-2", "JH033", "JH038", "JH054",
                            "JH086", "JH104", "JH115", "JH158", "JH-290",
                            "JH-304", "JH-305", "JH093", "JH-300")

#Subset for osi progression samples only
egfr_epi_opd <- subset(egfr_combined, subset = Specimen_ID %in% specimen_ids_to_subset)

#Subset for epithelial cells only
egfr_epi_opd <- subset(egfr_epi_opd, subset = General_Annot == "EPITHELIAL")

#Subset for InferCNV tumor cells only
egfr_epi_t_opd <- subset(egfr_epi_opd, subset = InfercnvR == "tumor")


#Summarize the number of epithelial cells by Specimen_ID
egfr_epi_opd_stats <- egfr_epi_t_opd@meta.data %>%
  group_by(Specimen_ID) %>%
  summarise(Count = n())

write.table(egfr_epi_opd_stats, "C:/Users/bbmorris1/Desktop/MDACC EGFR mutant scRNAseq Osi progression sample epithelial tumor cell count by Specimen ID.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Filter for Osi progression tumor samples with > 20 tumor cells
specimen_ids_to_subset <- c("JH033", "JH038",
                            "JH104", "JH-304", "JH-305")

#Subset for osi progression samples only
egfr_epi_opd <- subset(egfr_combined, subset = Specimen_ID %in% specimen_ids_to_subset)

#Subset for epithelial cells only
egfr_epi_opd <- subset(egfr_epi_opd, subset = General_Annot == "EPITHELIAL")

#Subset for InferCNV tumor cells only
egfr_epi_t_opd <- subset(egfr_epi_opd, subset = InfercnvR == "tumor")


#Change active identity to Specimen ID
Idents(egfr_epi_t_opd) <- "Specimen_ID"

#Draw KRT17 expression violin plot across osi progression samples
krt17_epi_t_opd <- VlnPlot(egfr_epi_t_opd, features = "KRT17", group.by = "Specimen_ID") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH-304" = "#FF62BF", 
                               "JH-305" = "#FF6A9A")) +
  labs(title = "Osi PD KRT17 expression")

krt17_epi_t_opd


###Annotate KRT17+ vs KRT17- cells by Specimen ID

#Label KRT17+ cells
egfr_epi_t_opd$KRT17_exp <- ifelse(FetchData(egfr_epi_t_opd, vars = "KRT17") > 0, 
                                   "KRT17+", "KRT17-")

head(egfr_epi_t_opd)

#Summarize the number of cells types by Specimen_ID
egfr_epi_opd_krt17_stats <- table(egfr_epi_t_opd@meta.data$Specimen_ID, egfr_epi_t_opd@meta.data$KRT17_exp)
egfr_epi_opd_krt17_stats

#Write data to file
write.table(egfr_epi_opd_krt17_stats, "C:/Users/bbmorris1/Desktop/MDACC EGFR mutant scRNAseq osimertinib PD KRT17 expression stats.txt", sep = "\t", row.names = TRUE, quote = FALSE)


###Summary barchart for KRT17+%

#Read in data
egfr_epi_t_opd_krt17_d <- MDACC_EGFR_mutant_scRNAseq_osimertinib_PD_KRT17_expression_stats_with_percentage

#Round to 2 decimal places for plotting
egfr_epi_t_opd_krt17_d$KRT17_perc_pos_round <- round(egfr_epi_t_opd_krt17_d$KRT17_perc_pos, 2)

#Draw barchart summarizing % of KRT17+ cells in each biopsy
egfr_epi_t_opd_krt17_bar <- ggplot(egfr_epi_t_opd_krt17_d, aes(x = egfr_epi_t_opd_krt17_d$Specimen_ID, y = egfr_epi_t_opd_krt17_d$KRT17_perc_pos_round, fill = egfr_epi_t_opd_krt17_d$Specimen_ID)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = KRT17_perc_pos_round), vjust = -0.5) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH-304" = "#FF62BF", 
                               "JH-305" = "#FF6A9A")) +
  scale_y_continuous("Percent KRT17+", limits = c(0, 100)) +
  labs(title = "Osimertinib PD")


egfr_epi_t_opd_krt17_bar


###KRT17/KRT5 co-expression analysis by Specimen_ID

#Label KRT5+ cells
egfr_epi_t_opd$KRT5_exp <- ifelse(FetchData(egfr_epi_t_opd, vars = "KRT5") > 0, 
                                  "KRT5+", "KRT5-")

head(egfr_epi_t_opd)


###Count the number of KRT17+ and KRT5+ cells for all MRD samples
egfr_epi_t_opd_krt17_krt5_counts <- egfr_epi_t_opd@meta.data %>%
  mutate(sample_id = Specimen_ID) %>%  
  group_by(Specimen_ID, KRT17_exp, KRT5_exp) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  pivot_wider(names_from = c(KRT17_exp, KRT5_exp), 
              values_from = cell_count,
              values_fill = 0)

print(egfr_epi_t_opd_krt17_krt5_counts)

write.table(egfr_epi_t_opd_krt17_krt5_counts, "C:/Users/bbmorris1/Desktop/MDACC EGFR mut scRNAseq Osi PD tumor cell KRT17 KRT5 co expression cell counts.txt", sep = '\t', row.names = FALSE)



###KRT17/KRT5 co-expression barchart

#Read in data
egfr_epi_t_opd_krt17_krt5 <- MDACC_EGFR_mut_scRNAseq_Osi_PD_tumor_cell_KRT17_KRT5_co_expression_barchart

#Draw barchart
egfr_epi_t_opd_krt17_krt5_bar <- ggplot(egfr_epi_t_opd_krt17_krt5, aes(x = egfr_epi_t_opd_krt17_krt5$Specimen_ID, y = egfr_epi_t_opd_krt17_krt5$Percentage, fill = egfr_epi_t_opd_krt17_krt5$Type)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("KRT17+ KRT5-" = "red",
                               "KRT17+ KRT5+" = "green")) +
  scale_y_continuous("Fraction") +
  labs(title = "Osi PD KRT17+ KRT5+ co-expression")


egfr_epi_t_opd_krt17_krt5_bar


#Stitch plots together for visualization
krt17_epi_t_opd | egfr_epi_t_opd_krt17_bar | egfr_epi_t_opd_krt17_krt5_bar #1400 x 400


##### Initial clinical cohort treatment naive, osimertinib MRD, and osimertinib progression %KRT17 expressing boxplot #####

#Read in data
tn_mrd_pd_krt17_perc <- MDACC_EGFR_mut_scRNAseq_percent_KRT17_expressing_naive_mrd_pd_barchart

#Facet order
tn_mrd_pd_krt17_perc$Type <- factor(tn_mrd_pd_krt17_perc$Type, levels = c("TN", "MRD", "PD"))

#Define groups for statistical comparisons
my_comparisons <- list(c("TN", "MRD"),
                       c("MRD", "PD"),
                       c("TN", "PD"))

#Draw boxplot
tn_mrd_pd_krt17_perc_p <- ggplot(tn_mrd_pd_krt17_perc, aes(x = tn_mrd_pd_krt17_perc$Type, y = tn_mrd_pd_krt17_perc$KRT17_perc_pos, fill = tn_mrd_pd_krt17_perc$Type)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  scale_y_continuous("Percentage KRT17+") +
  theme_bw() +
  theme(text = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(title = "Percentage KRT17+")


tn_mrd_pd_krt17_perc_p



##### Initial clinical cohort osi progression (PD) samples MET & SCLC gene expression plots #####

#Draw MET expression violin plot across osi progression samples
met_epi_t_opd <- VlnPlot(egfr_epi_t_opd, features = "MET", group.by = "Specimen_ID") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH-304" = "#FF62BF", 
                               "JH-305" = "#FF6A9A")) +
  labs(title = "Osi PD MET expression")

met_epi_t_opd #500x300


#Draw ASCL1 expression violin plot across osi progression samples
ascl1_epi_t_opd <- VlnPlot(egfr_epi_t_opd, features = "ASCL1", group.by = "Specimen_ID") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH-304" = "#FF62BF", 
                               "JH-305" = "#FF6A9A")) +
  labs(title = "Osi PD ASCL1 expression")

ascl1_epi_t_opd


#Draw TP63 expression violin plot across osi progression samples
tp63_epi_t_opd <- VlnPlot(egfr_epi_t_opd, features = "TP63", group.by = "Specimen_ID") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH-304" = "#FF62BF", 
                               "JH-305" = "#FF6A9A")) +
  labs(title = "Osi PD TP63 expression")

tp63_epi_t_opd



##### Initial clinical cohort treatment naive clinical data tile plot #####

#Create dataframe
tn_egfr <- data.frame(Specimen_ID = c("Biopsy-1", "JH064", "JH067", "JH095", "JH128", "JH139", "Lung-Tumor-10"),
                      EGFR_mut = c("L858R", "Exon20ins", "Exon19_del", "L858R", "L858R", "Exon19_del", "Exon20ins"),
                      Biopsy_location = c("Lung", "Lymph Node", "Liver", "Liver", "Liver", "Lung", "Lung"),
                      Pathology = c("Adenocarcinoma", "Adenocarcinoma", "Mucinous adenocarcinoma with NE features", "Unknown", "Moderately differentiated adenocarcinoma", "Moderately differentiated adenocarcinoma", "Well differentiated adenocarcinoma"))


#Set sample order
tn_egfr$Specimen_ID <- factor(tn_egfr$Specimen_ID, levels = c("Biopsy-1", "JH064", "JH067", "JH095", "JH128", "JH139", "Lung-Tumor-10"))

#Draw tile plot
tn_egfr_p <- ggplot(tn_egfr, aes(y = Specimen_ID)) +
  geom_tile(aes(x = "EGFR Mut", fill = EGFR_mut), color = "black") +
  scale_fill_manual(name = "EGFR Mut",
                    values = c("L858R" = "lightblue",
                               "Exon19_del" = "blue",
                               "Exon20ins" = "orange"),
                    guide = guide_legend(order = 1)) +
  new_scale_fill() + 
  geom_tile(aes(x = "Biopsy Location", fill = Biopsy_location), color = "black") +
  scale_fill_manual(name = "Biopsy Location",
                    values = c("Lung" = "dodgerblue",
                               "Lymph Node" = "khaki1",
                               "Liver" = "salmon"),
                    guide = guide_legend(order = 2)) +
  new_scale_fill() +
  geom_tile(aes(x = "Pathology", fill = Pathology), color = "black") +
  scale_fill_manual(name = "Pathology",
                    values = c("Adenocarcinoma" = "green3",
                               "Mucinous adenocarcinoma with NE features" = "lightseagreen",
                               "Moderately differentiated adenocarcinoma" = "seagreen1",
                               "Well differentiated adenocarcinoma" = "forestgreen",
                               "Unknown" = "gray50"),
                    guide = guide_legend(order = 3)) +
  #scale_x_discrete(limits = c("EGFR Mut", "Biopsy Location", "Pathology")) +
  scale_x_discrete(limits = c("Pathology", "Biopsy Location", "EGFR Mut"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)) +
  coord_flip()


tn_egfr_p #Plot 650 x 500



##### Initial clinical cohort osimertinib MRD clinical data tile plot #####

#Create dataframe
omrd_egfr <- data.frame(Specimen_ID = c("JH-297", "JH143", "Lung-Tumor-1", "Lung-Tumor-7", "Lung-Tumor-8", "NSTAR1-TumorA"),
                        EGFR_mut = c("L858R", "L858R", "L858R", "L858R", "Exon19_del", "L858R"),
                        Biopsy_location = c("Lung", "Lung", "Lung", "Lung", "Lung", "Lung"),
                        Pathology = c("Well differentiated adenocarcinoma", "Moderately differentiated adenocarcinoma", "Adenocarcinoma", "Adenocarcinoma", "Adenocarcinoma", "Adenocarcinoma"))


#Set sample order
omrd_egfr$Specimen_ID <- factor(omrd_egfr$Specimen_ID, levels = c("JH-297", "JH143", "Lung-Tumor-1", "Lung-Tumor-7", "Lung-Tumor-8", "NSTAR1-TumorA"))

#Draw tile plot
omrd_egfr_p <- ggplot(omrd_egfr, aes(y = Specimen_ID)) +
  geom_tile(aes(x = "EGFR Mut", fill = EGFR_mut), color = "black") +
  scale_fill_manual(name = "EGFR Mut",
                    values = c("L858R" = "lightblue",
                               "Exon19_del" = "blue"),
                    guide = guide_legend(order = 1)) +
  new_scale_fill() + 
  geom_tile(aes(x = "Biopsy Location", fill = Biopsy_location), color = "black") +
  scale_fill_manual(name = "Biopsy Location",
                    values = c("Lung" = "dodgerblue"),
                    guide = guide_legend(order = 2)) +
  new_scale_fill() +
  geom_tile(aes(x = "Pathology", fill = Pathology), color = "black") +
  scale_fill_manual(name = "Pathology",
                    values = c("Adenocarcinoma" = "green3",
                               "Mucinous adenocarcinoma with NE features" = "lightseagreen",
                               "Moderately differentiated adenocarcinoma" = "seagreen1",
                               "Well differentiated adenocarcinoma" = "forestgreen"),
                    guide = guide_legend(order = 3)) +
  scale_x_discrete(limits = c("Pathology", "Biopsy Location", "EGFR Mut"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)) +
  coord_flip()


omrd_egfr_p #Plot 650 x 500



##### Initial clinical cohort osimertinib resistant clinical data tile plot #####

#Create dataframe
ord_egfr <- data.frame(Specimen_ID = c("JH-304", "JH-305", "JH033", "JH038", "JH104"),
                       EGFR_mut = c("Exon19_del", "L858R", "Exon19_del", "Exon19_del", "L858R"),
                       Biopsy_location = c("Liver", "Lung", "Pleura", "Liver", "Liver"),
                       Pathology = c("Carcinoma with squamous features", "Adenocarcinoma", "Small cell lung cancer", "Adenocarcinoma", "Adenocarcinoma"),
                       Resistance_Mechanism = c("Unknown", "Unknown", "SCLC transformation", "MET amplification", "Unknown"))


#Set sample order
ord_egfr$Specimen_ID <- factor(ord_egfr$Specimen_ID, levels = c("JH-304", "JH-305", "JH033", "JH038", "JH104"))

#Draw tile plot
ord_egfr_p <- ggplot(ord_egfr, aes(y = Specimen_ID)) +
  geom_tile(aes(x = "EGFR Mut", fill = EGFR_mut), color = "black") +
  scale_fill_manual(name = "EGFR Mut",
                    values = c("L858R" = "lightblue",
                               "Exon19_del" = "blue"),
                    guide = guide_legend(order = 1)) +
  new_scale_fill() + 
  geom_tile(aes(x = "Biopsy Location", fill = Biopsy_location), color = "black") +
  scale_fill_manual(name = "Biopsy Location",
                    values = c("Lung" = "dodgerblue",
                               "Pleura" = "khaki1",
                               "Liver" = "salmon"),
                    guide = guide_legend(order = 2)) +
  new_scale_fill() +
  geom_tile(aes(x = "Pathology", fill = Pathology), color = "black") +
  scale_fill_manual(name = "Pathology",
                    values = c("Adenocarcinoma" = "green3",
                               "Carcinoma with squamous features" = "red",
                               "Small cell lung cancer" = "purple"),
                    guide = guide_legend(order = 3)) +
  new_scale_fill() +
  geom_tile(aes(x = "Resistance Mechanism", fill = Resistance_Mechanism), color = "black") +
  scale_fill_manual(name = "Resistance Mechanism",
                    values = c("Unknown" = "gray50",
                               "SCLC transformation" = "purple",
                               "MET amplification" = "chartreuse"),
                    guide = guide_legend(order = 4)) +
  scale_x_discrete(limits = c("Resistance Mechanism", "Pathology", "Biopsy Location", "EGFR Mut"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)) +
  coord_flip()


ord_egfr_p #Plot: 650 x 500

