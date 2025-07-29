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



##### HCC827 DMSO and Osi DTPC Lito Cell Cycle Analysis #####

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_annotated.Rds")

#Set SCT data as default assay
DefaultAssay(object = hcc827_combined) <- "SCT"

#Read in cell cycle genes table from Piro Lito
cc_lito <- Cell_cycle_genes_Lito

#Set G0 features
g0_features <- list(c(cc_lito$G0_oki))

#Add module score for G0
hcc827_combined <- AddModuleScore(hcc827_combined, features = g0_features,
                                  name = "G0")

#Set G1S features
g1s_features <- list(c(cc_lito$G1S))

#Add module score for G1S
hcc827_combined <- AddModuleScore(hcc827_combined, features = g1s_features,
                                  name = "G1S")

#Set S features
s_features <- list(c(cc_lito$S))

#Add module score for S
hcc827_combined <- AddModuleScore(hcc827_combined, features = s_features,
                                  name = "S")

#Set G2M features
g2m_features <- list(c(cc_lito$G2M))

#Add module score for G2M
hcc827_combined <- AddModuleScore(hcc827_combined, features = g2m_features,
                                  name = "G2M")

#Set M features
m_features <- list(c(cc_lito$M))

#Add module score for M
hcc827_combined <- AddModuleScore(hcc827_combined, features = m_features,
                                  name = "M")

#Set MG1 features
mg1_features <- list(c(cc_lito$MG1))

#Add module score for MG1
hcc827_combined <- AddModuleScore(hcc827_combined, features = mg1_features,
                                  name = "MG1")

###Rename metadata columns (removing "1" appended to every cell cycle entry)

#Copy the G01 column to G0
hcc827_combined@meta.data$G0 <- hcc827_combined@meta.data$G01

#Remove the old G01 column
hcc827_combined@meta.data$G01 <- NULL

#Copy the G1S1 column to G1S
hcc827_combined@meta.data$G1S <- hcc827_combined@meta.data$G1S1

#Remove the old G1S1 column
hcc827_combined@meta.data$G1S1 <- NULL

#Copy the S1 column to S
hcc827_combined@meta.data$S <- hcc827_combined@meta.data$S1

#Remove the old S column
hcc827_combined@meta.data$S1 <- NULL


#Copy the G2M1 column to G2M
hcc827_combined@meta.data$G2M <- hcc827_combined@meta.data$G2M1

#Remove the old G2M column
hcc827_combined@meta.data$G2M1 <- NULL


#Copy the M1 column to M
hcc827_combined@meta.data$M <- hcc827_combined@meta.data$M1

#Remove the old G2M column
hcc827_combined@meta.data$M1 <- NULL


#Copy the MG11 column to MG1
hcc827_combined@meta.data$MG1 <- hcc827_combined@meta.data$MG11

#Remove the old G2M column
hcc827_combined@meta.data$MG11 <- NULL




#Compare all cell cycle scores and assign the highest as cell cycle phase
hcc827_combined@meta.data$cell_cycle_expanded <- apply(hcc827_combined@meta.data[, c("G0", "G1S", "S", "G2M", "M", "MG1")], 1, function(x) {
  if (x["G0"] > x["G1S"] && x["G0"] > x["S"] && x["G0"] > x["G2M"] && x["G0"] > x["M"] && x["G0"] > x["MG1"]) {
    return("G0")
  }
  else if (x["G1S"] > x["G0"] && x["G1S"] > x["S"] && x["G1S"] > x["G2M"] && x["G1S"] > x["M"] && x["G1S"] > x["MG1"]) {
    return("G1S")
  }
  else if (x["S"] > x["G0"] && x["S"] > x["G1S"] && x["S"] > x["G2M"] && x["S"] > x["M"] && x["S"] > x["MG1"]) {
    return("S")
  }
  else if (x["G2M"] > x["G0"] && x["G2M"] > x["G1S"] && x["G2M"] > x["S"] && x["G2M"] > x["M"] && x["G2M"] > x["MG1"]) {
    return("G2M")
  }
  else if (x["M"] > x["G0"] && x["M"] > x["G1S"] && x["M"] > x["S"] && x["M"] > x["G2M"] && x["M"] > x["MG1"]) {
    return("M")
  }
  else if (x["MG1"] > x["G0"] && x["MG1"] > x["G1S"] && x["MG1"] > x["S"] && x["MG1"] > x["G2M"] && x["MG1"] > x["M"]) {
    return("MG1")
  }
  else {
    return(NA)  # Optional: Assign NA if no clear maximum
  }
})


#Save annotated dataset
SaveSeuratRds(hcc827_combined, "HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###UMAP

#Draw UMAP plot where cells are false colored by cell cycle state
hcc827_cc_umap <- DimPlot(hcc827_combined, reduction = "umap", group.by = "cell_cycle_expanded") +
  scale_color_manual(values = c("G1S" = "lightgreen", 
                                "S" = "blue", 
                                "G2M" = "purple",
                                "G0" = "red",
                                "M" = "magenta",
                                "MG1" = "gold")) +
  labs(title = "HCC827 Combined Lito Cell Cycle")

hcc827_cc_umap #550 x 400


#Read in percentage table
hcc827_cc_lito_perc <- HCC827_DMSO_HCC827_Osi_DTPC_lito_cellcyclescoring_percentage

#Set factor levels
hcc827_cc_lito_perc$Phase = factor(hcc827_cc_lito_perc$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
hcc827_cc_lito_perc_p <- ggplot(hcc827_cc_lito_perc, aes(x = Sample, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(label = paste0(hcc827_cc_lito_perc$Percentage, "%"), size = 3.5,
            color = "black", position = position_fill(vjust = 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC827 DMSO vs. Osi DTPC Lito Cell Cycle Phase", 
       y = "Fraction of cells")

hcc827_cc_lito_perc_p


###Cluster specific cell cycle states

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

#Create a table summarizing the number of cells in each phase
cell_cycle_counts_seurat <- table(hcc827_combined$cell_cycle_expanded, hcc827_combined$seurat_clusters)
cell_cycle_seurat_df <- as.data.frame(cell_cycle_counts_seurat)

write.table(cell_cycle_seurat_df, "HCC827 combined lito cellcyclescoring counts by seurat cluster.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Read in percentage table
hcc827_cc_lito_perc_seurat <- HCC827_combined_lito_cellcyclescoring_percentage_by_seurat_cluster

#Set factor levels
hcc827_cc_lito_perc_seurat$Phase = factor(hcc827_cc_lito_perc_seurat$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
hcc827_cc_lito_perc_seurat_p <- ggplot(hcc827_cc_lito_perc_seurat, aes(x = Cluster, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC827 Cluster Lito Cell Cycle Phase", 
       y = "Fraction of cells")

hcc827_cc_lito_perc_seurat_p


#Create a table summarizing the number of cells in each phase
cell_cycle_counts_seurat_sample <- table(hcc827_combined$cell_cycle_expanded, hcc827_combined$seurat_clusters, hcc827_combined$orig.ident)
cell_cycle_counts_seurat_sample_df <- as.data.frame(cell_cycle_counts_seurat_sample)

write.table(cell_cycle_counts_seurat_sample_df, "HCC827 DMSO and HCC827 Osi DTPC lito cellcyclescoring counts by seurat cluster.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###HCC827 DMSO Lito Cell Cycle by Seurat Cluster

#Read in percentage table
hcc827_dmso_cc_lito_perc_seurat <- HCC827_DMSO_lito_cell_cycle_scoring_percentage_by_seurat_cluster

#Set factor levels
hcc827_dmso_cc_lito_perc_seurat$Phase = factor(hcc827_dmso_cc_lito_perc_seurat$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
hcc827_dmso_cc_lito_perc_seurat_p <- ggplot(hcc827_dmso_cc_lito_perc_seurat, aes(x = Cluster, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC827 DMSO Cluster Lito Cell Cycle Phase", 
       y = "Fraction of cells")

hcc827_dmso_cc_lito_perc_seurat_p #600x300


###HCC827 Osi DTPC Lito Cell Cycle by Seurat Cluster

#Read in percentage table
hcc827_osi_dtpc_cc_lito_perc_seurat <- HCC827_Osi_DTPC_lito_cell_cycle_scoring_percentage_by_seurat_cluster

#Set factor levels
hcc827_osi_dtpc_cc_lito_perc_seurat$Phase = factor(hcc827_osi_dtpc_cc_lito_perc_seurat$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
hcc827_osi_dtpc_cc_lito_perc_seurat_p <- ggplot(hcc827_osi_dtpc_cc_lito_perc_seurat, aes(x = Cluster, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC827 Osi DTPC Cluster Lito Cell Cycle Phase", 
       y = "Fraction of cells")

hcc827_osi_dtpc_cc_lito_perc_seurat_p #600x300





##### HCC4006 DMSO and Osi DTPC Lito Cell Cycle Analysis #####

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_annotated.Rds")

#Set SCT data as default assay
DefaultAssay(object = hcc4006_combined) <- "SCT"

#Read in cell cycle genes table from Piro Lito
cc_lito <- Cell_cycle_genes_Lito

#Set G0 features
g0_features <- list(c(cc_lito$G0_oki))

#Add module score for G0
hcc4006_combined <- AddModuleScore(hcc4006_combined, features = g0_features,
                                   name = "G0")

#Set G1S features
g1s_features <- list(c(cc_lito$G1S))

#Add module score for G1S
hcc4006_combined <- AddModuleScore(hcc4006_combined, features = g1s_features,
                                   name = "G1S")

#Set S features
s_features <- list(c(cc_lito$S))

#Add module score for S
hcc4006_combined <- AddModuleScore(hcc4006_combined, features = s_features,
                                   name = "S")

#Set G2M features
g2m_features <- list(c(cc_lito$G2M))

#Add module score for G2M
hcc4006_combined <- AddModuleScore(hcc4006_combined, features = g2m_features,
                                   name = "G2M")

#Set M features
m_features <- list(c(cc_lito$M))

#Add module score for M
hcc4006_combined <- AddModuleScore(hcc4006_combined, features = m_features,
                                   name = "M")

#Set MG1 features
mg1_features <- list(c(cc_lito$MG1))

#Add module score for MG1
hcc4006_combined <- AddModuleScore(hcc4006_combined, features = mg1_features,
                                   name = "MG1")

###Rename metadata columns (removing "1" appended to every cell cycle entry)

#Copy the G01 column to G0
hcc4006_combined@meta.data$G0 <- hcc4006_combined@meta.data$G01

#Remove the old G01 column
hcc4006_combined@meta.data$G01 <- NULL

#Copy the G1S1 column to G1S
hcc4006_combined@meta.data$G1S <- hcc4006_combined@meta.data$G1S1

#Remove the old G1S1 column
hcc4006_combined@meta.data$G1S1 <- NULL

#Copy the S1 column to S
hcc4006_combined@meta.data$S <- hcc4006_combined@meta.data$S1

#Remove the old S column
hcc4006_combined@meta.data$S1 <- NULL


#Copy the G2M1 column to G2M
hcc4006_combined@meta.data$G2M <- hcc4006_combined@meta.data$G2M1

#Remove the old G2M column
hcc4006_combined@meta.data$G2M1 <- NULL


#Copy the M1 column to M
hcc4006_combined@meta.data$M <- hcc4006_combined@meta.data$M1

#Remove the old G2M column
hcc4006_combined@meta.data$M1 <- NULL


#Copy the MG11 column to MG1
hcc4006_combined@meta.data$MG1 <- hcc4006_combined@meta.data$MG11

#Remove the old G2M column
hcc4006_combined@meta.data$MG11 <- NULL




#Compare all cell cycle scores and assign the highest as cell cycle phase
hcc4006_combined@meta.data$cell_cycle_expanded <- apply(hcc4006_combined@meta.data[, c("G0", "G1S", "S", "G2M", "M", "MG1")], 1, function(x) {
  if (x["G0"] > x["G1S"] && x["G0"] > x["S"] && x["G0"] > x["G2M"] && x["G0"] > x["M"] && x["G0"] > x["MG1"]) {
    return("G0")
  }
  else if (x["G1S"] > x["G0"] && x["G1S"] > x["S"] && x["G1S"] > x["G2M"] && x["G1S"] > x["M"] && x["G1S"] > x["MG1"]) {
    return("G1S")
  }
  else if (x["S"] > x["G0"] && x["S"] > x["G1S"] && x["S"] > x["G2M"] && x["S"] > x["M"] && x["S"] > x["MG1"]) {
    return("S")
  }
  else if (x["G2M"] > x["G0"] && x["G2M"] > x["G1S"] && x["G2M"] > x["S"] && x["G2M"] > x["M"] && x["G2M"] > x["MG1"]) {
    return("G2M")
  }
  else if (x["M"] > x["G0"] && x["M"] > x["G1S"] && x["M"] > x["S"] && x["M"] > x["G2M"] && x["M"] > x["MG1"]) {
    return("M")
  }
  else if (x["MG1"] > x["G0"] && x["MG1"] > x["G1S"] && x["MG1"] > x["S"] && x["MG1"] > x["G2M"] && x["MG1"] > x["M"]) {
    return("MG1")
  }
  else {
    return(NA)  # Optional: Assign NA if no clear maximum
  }
})


#Save annotated dataset
SaveSeuratRds(hcc4006_combined, "HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###UMAP

#Draw UMAP plot where cells are false colored by cell cycle state
hcc4006_cc_umap <- DimPlot(hcc4006_combined, reduction = "umap", group.by = "cell_cycle_expanded") +
  scale_color_manual(values = c("G1S" = "lightgreen", 
                                "S" = "blue", 
                                "G2M" = "purple",
                                "G0" = "red",
                                "M" = "magenta",
                                "MG1" = "gold")) +
  labs(title = "HCC4006 Combined Lito Cell Cycle")

hcc4006_cc_umap #550 x 400


#Create a table summarizing the number of cells in each phase
cell_cycle_counts <- table(hcc4006_combined$cell_cycle_expanded, hcc4006_combined$orig.ident)
cell_cycle_df <- as.data.frame(cell_cycle_counts)

write.table(cell_cycle_df, "HCC4006 DMSO HCC4006 Osi DTPC lito cellcyclescoring counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Read in percentage table
hcc4006_cc_lito_perc <- HCC4006_DMSO_HCC4006_Osi_DTPC_lito_cellcyclescoring_percentage

#Set factor levels
hcc4006_cc_lito_perc$Phase = factor(hcc4006_cc_lito_perc$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
hcc4006_cc_lito_perc_p <- ggplot(hcc4006_cc_lito_perc, aes(x = Sample, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(label = paste0(hcc4006_cc_lito_perc$Percentage, "%"), size = 3.5,
            color = "black", position = position_fill(vjust = 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC4006 DMSO vs. Osi DTPC Lito Cell Cycle Phase", 
       y = "Fraction of cells")

hcc4006_cc_lito_perc_p #450 x 600


###Cluster specific cell cycle states

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

#Create a table summarizing the number of cells in each phase
cell_cycle_counts_seurat <- table(hcc4006_combined$cell_cycle_expanded, hcc4006_combined$seurat_clusters)
cell_cycle_seurat_df <- as.data.frame(cell_cycle_counts_seurat)

write.table(cell_cycle_seurat_df, "HCC4006 combined lito cellcyclescoring counts by seurat cluster.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Combined (DMSO + Osi DTPC) seurat cluster Lito cell cycle states

#Read in percentage table
hcc4006_cc_lito_perc_seurat <- HCC4006_combined_lito_cellcyclescoring_percentage_by_seurat_cluster

#Set factor levels
hcc4006_cc_lito_perc_seurat$Phase = factor(hcc4006_cc_lito_perc_seurat$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
hcc4006_cc_lito_perc_seurat_p <- ggplot(hcc4006_cc_lito_perc_seurat, aes(x = Cluster, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC4006 Cluster Lito Cell Cycle Phase", 
       y = "Fraction of cells")

hcc4006_cc_lito_perc_seurat_p



#Create a table summarizing the number of cells in each phase
cell_cycle_counts_seurat_sample <- table(hcc4006_combined$cell_cycle_expanded, hcc4006_combined$seurat_clusters, hcc4006_combined$orig.ident)
cell_cycle_counts_seurat_sample_df <- as.data.frame(cell_cycle_counts_seurat_sample)

write.table(cell_cycle_counts_seurat_sample_df, "HCC4006 DMSO and HCC4006 Osi DTPC lito cellcyclescoring counts by seurat cluster.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###HCC4006 DMSO Lito Cell Cycle by Seurat Cluster

#Read in percentage table
hcc4006_dmso_cc_lito_perc_seurat <- HCC4006_DMSO_lito_cell_cycle_scoring_percentage_by_seurat_cluster

#Set factor levels
hcc4006_dmso_cc_lito_perc_seurat$Phase = factor(hcc4006_dmso_cc_lito_perc_seurat$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
hcc4006_dmso_cc_lito_perc_seurat_p <- ggplot(hcc4006_dmso_cc_lito_perc_seurat, aes(x = Cluster, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC4006 DMSO Cluster Lito Cell Cycle Phase", 
       y = "Fraction of cells")

hcc4006_dmso_cc_lito_perc_seurat_p #600x300



###HCC4006 Osi DTPC Lito Cell Cycle by Seurat Cluster
#Read in percentage table
hcc4006_osi_dtpc_cc_lito_perc_seurat <- HCC4006_Osi_DTPC_lito_cell_cycle_scoring_percentage_by_seurat_cluster

#Set factor levels
hcc4006_osi_dtpc_cc_lito_perc_seurat$Phase = factor(hcc4006_osi_dtpc_cc_lito_perc_seurat$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
hcc4006_osi_dtpc_cc_lito_perc_seurat_p <- ggplot(hcc4006_osi_dtpc_cc_lito_perc_seurat, aes(x = Cluster, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC4006 Osi DTPC Cluster Lito Cell Cycle Phase", 
       y = "Fraction of cells")

hcc4006_osi_dtpc_cc_lito_perc_seurat_p #600x300



##### H1975 DMSO and Osi DTPC Lito Cell Cycle Analysis #####

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_annotated.Rds")

#Set SCT data as default assay
DefaultAssay(object = h1975_combined) <- "SCT"

#Read in cell cycle genes table from Piro Lito
cc_lito <- Cell_cycle_genes_Lito

#Set G0 features
g0_features <- list(c(cc_lito$G0_oki))

#Add module score for G0
h1975_combined <- AddModuleScore(h1975_combined, features = g0_features,
                                 name = "G0")

#Set G1S features
g1s_features <- list(c(cc_lito$G1S))

#Add module score for G1S
h1975_combined <- AddModuleScore(h1975_combined, features = g1s_features,
                                 name = "G1S")

#Set S features
s_features <- list(c(cc_lito$S))

#Add module score for S
h1975_combined <- AddModuleScore(h1975_combined, features = s_features,
                                 name = "S")

#Set G2M features
g2m_features <- list(c(cc_lito$G2M))

#Add module score for G2M
h1975_combined <- AddModuleScore(h1975_combined, features = g2m_features,
                                 name = "G2M")

#Set M features
m_features <- list(c(cc_lito$M))

#Add module score for M
h1975_combined <- AddModuleScore(h1975_combined, features = m_features,
                                 name = "M")

#Set MG1 features
mg1_features <- list(c(cc_lito$MG1))

#Add module score for MG1
h1975_combined <- AddModuleScore(h1975_combined, features = mg1_features,
                                 name = "MG1")

###Rename metadata columns (removing "1" appended to every cell cycle entry)

#Copy the G01 column to G0
h1975_combined@meta.data$G0 <- h1975_combined@meta.data$G01

#Remove the old G01 column
h1975_combined@meta.data$G01 <- NULL

#Copy the G1S1 column to G1S
h1975_combined@meta.data$G1S <- h1975_combined@meta.data$G1S1

#Remove the old G1S1 column
h1975_combined@meta.data$G1S1 <- NULL

#Copy the S1 column to S
h1975_combined@meta.data$S <- h1975_combined@meta.data$S1

#Remove the old S column
h1975_combined@meta.data$S1 <- NULL


#Copy the G2M1 column to G2M
h1975_combined@meta.data$G2M <- h1975_combined@meta.data$G2M1

#Remove the old G2M column
h1975_combined@meta.data$G2M1 <- NULL


#Copy the M1 column to M
h1975_combined@meta.data$M <- h1975_combined@meta.data$M1

#Remove the old G2M column
h1975_combined@meta.data$M1 <- NULL


#Copy the MG11 column to MG1
h1975_combined@meta.data$MG1 <- h1975_combined@meta.data$MG11

#Remove the old G2M column
h1975_combined@meta.data$MG11 <- NULL




#Compare all cell cycle scores and assign the highest as cell cycle phase
h1975_combined@meta.data$cell_cycle_expanded <- apply(h1975_combined@meta.data[, c("G0", "G1S", "S", "G2M", "M", "MG1")], 1, function(x) {
  if (x["G0"] > x["G1S"] && x["G0"] > x["S"] && x["G0"] > x["G2M"] && x["G0"] > x["M"] && x["G0"] > x["MG1"]) {
    return("G0")
  }
  else if (x["G1S"] > x["G0"] && x["G1S"] > x["S"] && x["G1S"] > x["G2M"] && x["G1S"] > x["M"] && x["G1S"] > x["MG1"]) {
    return("G1S")
  }
  else if (x["S"] > x["G0"] && x["S"] > x["G1S"] && x["S"] > x["G2M"] && x["S"] > x["M"] && x["S"] > x["MG1"]) {
    return("S")
  }
  else if (x["G2M"] > x["G0"] && x["G2M"] > x["G1S"] && x["G2M"] > x["S"] && x["G2M"] > x["M"] && x["G2M"] > x["MG1"]) {
    return("G2M")
  }
  else if (x["M"] > x["G0"] && x["M"] > x["G1S"] && x["M"] > x["S"] && x["M"] > x["G2M"] && x["M"] > x["MG1"]) {
    return("M")
  }
  else if (x["MG1"] > x["G0"] && x["MG1"] > x["G1S"] && x["MG1"] > x["S"] && x["MG1"] > x["G2M"] && x["MG1"] > x["M"]) {
    return("MG1")
  }
  else {
    return(NA)  # Optional: Assign NA if no clear maximum
  }
})


#Save annotated dataset
SaveSeuratRds(h1975_combined, "H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###UMAP

#Draw UMAP plot where cells are false colored by cell cycle state
h1975_cc_umap <- DimPlot(h1975_combined, reduction = "umap", group.by = "cell_cycle_expanded") +
  scale_color_manual(values = c("G1S" = "lightgreen", 
                                "S" = "blue", 
                                "G2M" = "purple",
                                "G0" = "red",
                                "M" = "magenta",
                                "MG1" = "gold")) +
  labs(title = "H1975 Combined Lito Cell Cycle")

h1975_cc_umap #550 x 400


#Create a table summarizing the number of cells in each phase
cell_cycle_counts <- table(h1975_combined$cell_cycle_expanded, h1975_combined$orig.ident)
cell_cycle_df <- as.data.frame(cell_cycle_counts)

write.table(cell_cycle_df, "H1975 DMSO H1975 Osi DTPC lito cellcyclescoring counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Read in percentage table
h1975_cc_lito_perc <- H1975_DMSO_H1975_Osi_DTPC_lito_cellcyclescoring_percentage

#Set factor levels
h1975_cc_lito_perc$Phase = factor(h1975_cc_lito_perc$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
h1975_cc_lito_perc_p <- ggplot(h1975_cc_lito_perc, aes(x = Sample, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(label = paste0(h1975_cc_lito_perc$Percentage, "%"), size = 3.5,
            color = "black", position = position_fill(vjust = 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "H1975 DMSO vs. Osi DTPC Lito Cell Cycle Phase", 
       y = "Fraction of cells")

h1975_cc_lito_perc_p


###Cluster specific cell cycle states

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

#Create a table summarizing the number of cells in each phase
cell_cycle_counts_seurat <- table(h1975_combined$cell_cycle_expanded, h1975_combined$seurat_clusters)
cell_cycle_seurat_df <- as.data.frame(cell_cycle_counts_seurat)

write.table(cell_cycle_seurat_df, "H1975 combined lito cellcyclescoring counts by seurat cluster.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Combined (DMSO + Osi DTPC) seurat cluster Lito cell cycle states

#Read in percentage table
h1975_cc_lito_perc_seurat <- H1975_combined_lito_cellcyclescoring_percentage_by_seurat_cluster

#Set factor levels
h1975_cc_lito_perc_seurat$Phase = factor(h1975_cc_lito_perc_seurat$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
h1975_cc_lito_perc_seurat_p <- ggplot(h1975_cc_lito_perc_seurat, aes(x = Cluster, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "H1975 Cluster Lito Cell Cycle Phase", 
       y = "Fraction of cells")

h1975_cc_lito_perc_seurat_p



#Create a table summarizing the number of cells in each phase
cell_cycle_counts_seurat_sample <- table(h1975_combined$cell_cycle_expanded, h1975_combined$seurat_clusters, h1975_combined$orig.ident)
cell_cycle_counts_seurat_sample_df <- as.data.frame(cell_cycle_counts_seurat_sample)

write.table(cell_cycle_counts_seurat_sample_df, "H1975 DMSO and H1975 Osi DTPC lito cellcyclescoring counts by seurat cluster.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###H1975 DMSO Lito Cell Cycle by Seurat Cluster

#Read in percentage table
h1975_dmso_cc_lito_perc_seurat <- H1975_DMSO_lito_cell_cycle_scoring_percentage_by_seurat_cluster

#Set factor levels
h1975_dmso_cc_lito_perc_seurat$Phase = factor(h1975_dmso_cc_lito_perc_seurat$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
h1975_dmso_cc_lito_perc_seurat_p <- ggplot(h1975_dmso_cc_lito_perc_seurat, aes(x = Cluster, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "H1975 DMSO Cluster Lito Cell Cycle Phase", 
       y = "Fraction of cells")

h1975_dmso_cc_lito_perc_seurat_p #600x300


###H1975 Osi DTPC Lito Cell Cycle by Seurat Cluster

#Read in percentage table
h1975_osi_dtpc_cc_lito_perc_seurat <- H1975_Osi_DTPC_lito_cell_cycle_scoring_percentage_by_seurat_cluster

#Set factor levels
h1975_osi_dtpc_cc_lito_perc_seurat$Phase = factor(h1975_osi_dtpc_cc_lito_perc_seurat$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot stacked bar chart
h1975_osi_dtpc_cc_lito_perc_seurat_p <- ggplot(h1975_osi_dtpc_cc_lito_perc_seurat, aes(x = Cluster, y = Percentage, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "H1975 Osi DTPC Cluster Lito Cell Cycle Phase", 
       y = "Fraction of cells")

h1975_osi_dtpc_cc_lito_perc_seurat_p #600x300









