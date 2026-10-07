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
library(ggpubr)
library(forcats)



##### HCC827 DMSO & Osi DTPC Seurat transcriptional cluster plot #####

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


#Plot UMAP with seurat cluster pseudocoloring
hcc827_sc_dim <- DimPlot(hcc827_combined, group.by = "seurat_clusters", 
                         split.by = "orig.ident",
                         pt.size = 1) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "Seurat clusters") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

hcc827_sc_dim



##### HCC4006 DMSO & Osi DTPC Seurat transcriptional cluster plot #####

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

#Plot UMAP with seurat cluster pseudocoloring
hcc4006_sc_dim <- DimPlot(hcc4006_combined, group.by = "seurat_clusters", 
                          split.by = "orig.ident",
                          pt.size = 1) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "Seurat clusters") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank())

hcc4006_sc_dim



##### H1975 DMSO & Osi DTPC Seurat transcriptional cluster plot #####

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


#Plot UMAP with seurat cluster pseudocoloring
h1975_sc_dim <- DimPlot(h1975_combined, group.by = "seurat_clusters", 
                        split.by = "orig.ident",
                        pt.size = 1) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "Seurat clusters") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank())

h1975_sc_dim





