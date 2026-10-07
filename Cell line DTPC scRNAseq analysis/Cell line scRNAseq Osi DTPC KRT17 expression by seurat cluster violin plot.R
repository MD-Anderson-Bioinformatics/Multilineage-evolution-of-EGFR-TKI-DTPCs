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
library(ggpubr)



##### HCC827 DMSO and HCC827 Osi DTPC integrated analysis #####

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

#Subset DTPCs
hcc827_dtpc <- subset(hcc827_combined, subset = orig.ident == "HCC827_Osi_DTPC")

#Draw violin plot for DTPC KRT17 expression by seurat cluster
krt17_vln_p <- VlnPlot(hcc827_dtpc, features = "KRT17") +
  theme(legend.position = "none")

krt17_vln_p


##### HCC4006 DMSO and HCC4006 Osi DTPC integrated analysis #####

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

#Subset DTPCs
hcc4006_dtpc <- subset(hcc4006_combined, subset = orig.ident == "HCC4006_Osi_DTPC")

#Draw violin plot for DTPC KRT17 expression by seurat cluster
krt17_vln_p <- VlnPlot(hcc4006_dtpc, features = "KRT17") +
  theme(legend.position = "none")

krt17_vln_p


##### H1975 DMSO and H1975 Osi DTPC integrated analysis #####

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

#Subset DTPCs
h1975_dtpc <- subset(h1975_combined, subset = orig.ident == "H1975_Osi_DTPC")

#Draw violin plot for DTPC KRT17 expression by seurat cluster
krt17_vln_p <- VlnPlot(h1975_combined, features = "KRT17") +
  theme(legend.position = "none")

krt17_vln_p


