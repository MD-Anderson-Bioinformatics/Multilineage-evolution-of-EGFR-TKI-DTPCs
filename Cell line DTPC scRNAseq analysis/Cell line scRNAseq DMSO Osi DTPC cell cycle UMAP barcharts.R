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
library(patchwork)


##### HCC827 DMSO and HCC827 Osi DTPC Lito Cell Cycle UMAP #####

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

###Lito cell cycle faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc827_combined, "umap"))
umap_data$lito_cc <- hcc827_combined@meta.data$cell_cycle_expanded
umap_data$orig.ident <- hcc827_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
hcc827_lito_cc_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = lito_cc)) +
  geom_point() +
  scale_color_manual(values = c("G1S" = "lightgreen", 
                                "S" = "blue", 
                                "G2M" = "purple",
                                "G0" = "red",
                                "M" = "magenta",
                                "MG1" = "gold")) +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "Cell Cycle") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"),
        legend.position = "none")

hcc827_lito_cc_umap_f #Plot 600 x 300


####Lito cell cycle barchart (DMSO vs Osi DTPC)

#Read in percentage table
hcc827_cc_lito_perc <- HCC827_DMSO_HCC827_Osi_DTPC_lito_cellcyclescoring_percentage

#Round values to one decimal place for plotting
hcc827_cc_lito_perc$Percentage_round <- round(hcc827_cc_lito_perc$Percentage, 1)

#Set factor levels
hcc827_cc_lito_perc$Phase = factor(hcc827_cc_lito_perc$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot bar chart
hcc827_cc_lito_perc_p <- ggplot(hcc827_cc_lito_perc, aes(x = Sample, y = Percentage_round, fill = Phase)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(label = paste0(hcc827_cc_lito_perc$Percentage_round), size = 3.5,
            color = "black", position = position_dodge(width = 0.9),
            vjust = -0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_y_continuous(limits = c(0, 65)) + 
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC827 Lito Cell Cycle Phase", 
       y = "Percentage of cells")

hcc827_cc_lito_perc_p #Plot 550 x 300



####Lito cell cycle barchart (DMSO only by seurat cluster)

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
  labs(title = "HCC827 DMSO", 
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
  labs(title = "HCC827 Osi DTPC", 
       y = "Fraction of cells")

hcc827_osi_dtpc_cc_lito_perc_seurat_p #600x300



##### HCC4006 DMSO and HCC4006 Osi DTPC Lito Cell Cycle UMAP #####

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###Lito cell cycle faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc4006_combined, "umap"))
umap_data$lito_cc <- hcc4006_combined@meta.data$cell_cycle_expanded
umap_data$orig.ident <- hcc4006_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
hcc4006_lito_cc_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = lito_cc)) +
  geom_point() +
  scale_color_manual(values = c("G1S" = "lightgreen", 
                                "S" = "blue", 
                                "G2M" = "purple",
                                "G0" = "red",
                                "M" = "magenta",
                                "MG1" = "gold")) +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "Cell Cycle") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"),
        legend.position = "none")

hcc4006_lito_cc_umap_f #Plot 600 x 300


####Lito cell cycle barchart (DMSO vs Osi DTPC)

#Read in percentage table
hcc4006_cc_lito_perc <- HCC4006_DMSO_HCC4006_Osi_DTPC_lito_cellcyclescoring_percentage

#Round values to one decimal place for plotting
hcc4006_cc_lito_perc$Percentage_round <- round(hcc4006_cc_lito_perc$Percentage, 1)

#Set factor levels
hcc4006_cc_lito_perc$Phase = factor(hcc4006_cc_lito_perc$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot bar chart
hcc4006_cc_lito_perc_p <- ggplot(hcc4006_cc_lito_perc, aes(x = Sample, y = Percentage_round, fill = Phase)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(label = paste0(hcc4006_cc_lito_perc$Percentage_round), size = 3.5,
            color = "black", position = position_dodge(width = 0.9),
            vjust = -0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_y_continuous(limits = c(0, 65)) + 
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "HCC4006 Lito Cell Cycle Phase", 
       y = "Percentage of cells")

hcc4006_cc_lito_perc_p #Plot 550 x 300


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
  labs(title = "HCC4006 DMSO", 
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
  labs(title = "HCC4006 Osi DTPC", 
       y = "Fraction of cells")

hcc4006_osi_dtpc_cc_lito_perc_seurat_p #600x300



##### H1975 DMSO and H1975 Osi DTPC Lito Cell Cycle UMAP #####

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###Lito cell cycle faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(h1975_combined, "umap"))
umap_data$lito_cc <- h1975_combined@meta.data$cell_cycle_expanded
umap_data$orig.ident <- h1975_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
h1975_lito_cc_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = lito_cc)) +
  geom_point() +
  scale_color_manual(values = c("G1S" = "lightgreen", 
                                "S" = "blue", 
                                "G2M" = "purple",
                                "G0" = "red",
                                "M" = "magenta",
                                "MG1" = "gold")) +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "Cell Cycle") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"),
        legend.position = "none")

h1975_lito_cc_umap_f #Plot 600 x 300


####Lito cell cycle barchart (DMSO vs Osi DTPC)

#Read in percentage table
h1975_cc_lito_perc <- H1975_DMSO_H1975_Osi_DTPC_lito_cellcyclescoring_percentage

#Round values to one decimal place for plotting
h1975_cc_lito_perc$Percentage_round <- round(h1975_cc_lito_perc$Percentage, 1)

#Set factor levels
h1975_cc_lito_perc$Phase = factor(h1975_cc_lito_perc$Phase, levels=c("G1S", "S", "G2M", "M", "MG1", "G0"))

#Plot bar chart
h1975_cc_lito_perc_p <- ggplot(h1975_cc_lito_perc, aes(x = Sample, y = Percentage_round, fill = Phase)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(label = paste0(h1975_cc_lito_perc$Percentage_round), size = 3.5,
            color = "black", position = position_dodge(width = 0.9),
            vjust = -0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = 10)) +
  scale_y_continuous(limits = c(0, 80)) + 
  scale_fill_manual(values = c("G0" = "red", 
                               "G1S" = "lightgreen", 
                               "S" = "blue",
                               "G2M" = "purple",
                               "MG1" = "gold",
                               "M" = "magenta")) +
  labs(title = "H1975 Lito Cell Cycle Phase", 
       y = "Percentage of cells")

h1975_cc_lito_perc_p #Plot 550 x 300


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
  labs(title = "H1975 DMSO", 
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
  labs(title = "H1975 Osi DTPC", 
       y = "Fraction of cells")

h1975_osi_dtpc_cc_lito_perc_seurat_p #600x300






