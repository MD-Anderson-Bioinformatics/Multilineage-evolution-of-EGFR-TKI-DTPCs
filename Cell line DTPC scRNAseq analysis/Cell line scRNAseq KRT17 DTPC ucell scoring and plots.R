library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(UCell)
library(SeuratData)
library(SeuratObject)
library(SeuratExtend)
library(SeuratDisk)
library(ggpubr)
library(forcats)



##### HCC827 KRT17+ Osi DTPC integrated UCell scoring analysis #####

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


#Define genes of interest
genes_of_interest <- c("KRT17")

#Extract expression data for the three genes
expression_data <- FetchData(hcc827_combined, vars = genes_of_interest)

#Convert to dataframe
expression_data <- as.data.frame(expression_data)

#Define thresholds for gene expression (adjust as needed)
threshold <- 0.1
expression_data <- expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT17 > threshold, 1, 0),
    ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc827_combined <- AddMetaData(hcc827_combined, expression_data$ExpressionCategory, col.name = "ExpressionCategory")

#Build custom color key for ExpressionCategory
existing_column <- "ExpressionCategory"

label_map <- c(
  "0" = "KRT17-",
  "1" = "KRT17+")


#Assign ExpressionCategory data to separate metadata variable
metadata <- hcc827_combined@meta.data

#Match surface target co-expression labels to ExpressionCategory binary values using key
metadata$ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                   label_map[metadata[[existing_column]]], 
                                   "Other")  # Default "Other" if value doesn't match

#Add new KRT17 labels back to metadata
hcc827_combined <- AddMetaData(hcc827_combined, metadata$ExpressionLabel, col.name = "ExpressionLabel")

#Check new metadata entry
head(hcc827_combined@meta.data)


###Aberrant basaloid scoring (Wang)

#Read in aberrant basaloid genes
wang_ab_basaloid_genes <- Wang_aberrant_basaloid_signature

#Set aberrant basaloid features
wang_ab_basaloid_genes_features <- list(c(wang_ab_basaloid_genes$Gene))

#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=wang_ab_basaloid_genes_features, name="ab_basaloid_wang")


#Subset for only Osi DTPCs
hcc827_dtpc <- subset(hcc827_combined, orig.ident == "HCC827_Osi_DTPC")


###Aberrant basaloid violin plot
hcc827_krt17_dtpc_ab_basaloid_vp <- VlnPlot(hcc827_dtpc, features = "signature_1ab_basaloid_wang", group.by = "ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "HCC827 Osi DTPC \nAberrant Basaloid")

hcc827_krt17_dtpc_ab_basaloid_vp #Plot: 350 x 500


###Wong Embryonic Stem Cell Core

#Read in Wong ESC genes
wong_esc_genes <- Wong_embryonic_stem_cell_core

#Set Wong ESC gene features
wong_esc_genes_features <- list(c(wong_esc_genes$Gene))

#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=wong_esc_genes_features, name="wong_esc")


#Subset for only Osi DTPCs
hcc827_dtpc <- subset(hcc827_combined, orig.ident == "HCC827_Osi_DTPC")


###Wong Embryonic Stem Cell core violin plot
hcc827_krt17_dtpc_wong_esc_vp <- VlnPlot(hcc827_dtpc, features = "signature_1wong_esc", group.by = "ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "HCC827 Osi DTPC \nWong ESC")

hcc827_krt17_dtpc_wong_esc_vp #Plot: 350 x 500



###EMT

#Read in Hallmark EMT genes
emt_genes <- Hallmark_EMT_genes

#Set Hallmark EMT gene features
emt_genes_features <- list(c(emt_genes$Gene))

#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=emt_genes_features, name="emt")


#Subset for only Osi DTPCs
hcc827_dtpc <- subset(hcc827_combined, orig.ident == "HCC827_Osi_DTPC")


###Hallmark EMT genes violin plot
hcc827_krt17_dtpc_emt_vp <- VlnPlot(hcc827_dtpc, features = "signature_1emt", group.by = "ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "HCC827 Osi DTPC \nEMT")

hcc827_krt17_dtpc_emt_vp #Plot: 350 x 500



##### HCC4006 KRT17+ Osi DTPC integrated UCell scoring analysis #####

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


#Define genes of interest
genes_of_interest <- c("KRT17")

#Extract expression data for the three genes
expression_data <- FetchData(hcc4006_combined, vars = genes_of_interest)

#Convert to dataframe
expression_data <- as.data.frame(expression_data)

#Define thresholds for gene expression (adjust as needed)
threshold <- 0.1
expression_data <- expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT17 > threshold, 1, 0),
    ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc4006_combined <- AddMetaData(hcc4006_combined, expression_data$ExpressionCategory, col.name = "ExpressionCategory")

#Build custom color key for ExpressionCategory
existing_column <- "ExpressionCategory"

label_map <- c(
  "0" = "KRT17-",
  "1" = "KRT17+")


#Assign ExpressionCategory data to separate metadata variable
metadata <- hcc4006_combined@meta.data

#Match surface target co-expression labels to ExpressionCategory binary values using key
metadata$ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                   label_map[metadata[[existing_column]]], 
                                   "Other")  # Default "Other" if value doesn't match

#Add new KRT17 labels back to metadata
hcc4006_combined <- AddMetaData(hcc4006_combined, metadata$ExpressionLabel, col.name = "ExpressionLabel")

#Check new metadata entry
head(hcc4006_combined@meta.data)


###Aberrant basaloid scoring (Wang)

#Read in aberrant basaloid genes
wang_ab_basaloid_genes <- Wang_aberrant_basaloid_signature

#Set aberrant basaloid features
wang_ab_basaloid_genes_features <- list(c(wang_ab_basaloid_genes$Gene))

#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=wang_ab_basaloid_genes_features, name="ab_basaloid_wang")


#Subset for only Osi DTPCs
hcc4006_dtpc <- subset(hcc4006_combined, orig.ident == "HCC4006_Osi_DTPC")


###Aberrant basaloid violin plot
hcc4006_krt17_dtpc_ab_basaloid_vp <- VlnPlot(hcc4006_dtpc, features = "signature_1ab_basaloid_wang", group.by = "ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "HCC4006 Osi DTPC \nAberrant Basaloid")

hcc4006_krt17_dtpc_ab_basaloid_vp #Plot: 350 x 500



###Wong Embryonic Stem Cell Core

#Read in Wong ESC genes
wong_esc_genes <- Wong_embryonic_stem_cell_core

#Set Wong ESC gene features
wong_esc_genes_features <- list(c(wong_esc_genes$Gene))

#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=wong_esc_genes_features, name="wong_esc")


#Subset for only Osi DTPCs
hcc4006_dtpc <- subset(hcc4006_combined, orig.ident == "HCC4006_Osi_DTPC")


###Wong Embryonic Stem Cell core violin plot
hcc4006_krt17_dtpc_wong_esc_vp <- VlnPlot(hcc4006_dtpc, features = "signature_1wong_esc", group.by = "ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "HCC4006 Osi DTPC \nWong ESC")

hcc4006_krt17_dtpc_wong_esc_vp #Plot: 350 x 500


###EMT

#Read in Hallmark EMT genes
emt_genes <- Hallmark_EMT_genes

#Set Hallmark EMT gene features
emt_genes_features <- list(c(emt_genes$Gene))

#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=emt_genes_features, name="emt")


#Subset for only Osi DTPCs
hcc4006_dtpc <- subset(hcc4006_combined, orig.ident == "HCC4006_Osi_DTPC")


###Hallmark EMT genes violin plot
hcc4006_krt17_dtpc_emt_vp <- VlnPlot(hcc4006_dtpc, features = "signature_1emt", group.by = "ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "HCC4006 Osi DTPC \nEMT")

hcc4006_krt17_dtpc_emt_vp #Plot: 350 x 500



##### H1975 KRT17+ Osi DTPC integrated UCell scoring analysis #####

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


#Define genes of interest
genes_of_interest <- c("KRT17")

#Extract expression data for the three genes
expression_data <- FetchData(h1975_combined, vars = genes_of_interest)

#Convert to dataframe
expression_data <- as.data.frame(expression_data)

#Define thresholds for gene expression (adjust as needed)
threshold <- 0.1
expression_data <- expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT17 > threshold, 1, 0),
    ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
h1975_combined <- AddMetaData(h1975_combined, expression_data$ExpressionCategory, col.name = "ExpressionCategory")

#Build custom color key for ExpressionCategory
existing_column <- "ExpressionCategory"

label_map <- c(
  "0" = "KRT17-",
  "1" = "KRT17+")


#Assign ExpressionCategory data to separate metadata variable
metadata <- h1975_combined@meta.data

#Match surface target co-expression labels to ExpressionCategory binary values using key
metadata$ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                   label_map[metadata[[existing_column]]], 
                                   "Other")  # Default "Other" if value doesn't match

#Add new KRT17 labels back to metadata
h1975_combined <- AddMetaData(h1975_combined, metadata$ExpressionLabel, col.name = "ExpressionLabel")

#Check new metadata entry
head(h1975_combined@meta.data)


###Aberrant basaloid scoring (Wang)

#Read in aberrant basaloid genes
wang_ab_basaloid_genes <- Wang_aberrant_basaloid_signature

#Set aberrant basaloid features
wang_ab_basaloid_genes_features <- list(c(wang_ab_basaloid_genes$Gene))

#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=wang_ab_basaloid_genes_features, name="ab_basaloid_wang")


#Subset for only Osi DTPCs
h1975_dtpc <- subset(h1975_combined, orig.ident == "H1975_Osi_DTPC")


###Aberrant basaloid violin plot
H1975_krt17_dtpc_ab_basaloid_vp <- VlnPlot(h1975_dtpc, features = "signature_1ab_basaloid_wang", group.by = "ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "H1975 Osi DTPC \nAberrant Basaloid")

H1975_krt17_dtpc_ab_basaloid_vp #Plot: 350 x 500



###Wong Embryonic Stem Cell Core

#Read in Wong ESC genes
wong_esc_genes <- Wong_embryonic_stem_cell_core

#Set Wong ESC gene features
wong_esc_genes_features <- list(c(wong_esc_genes$Gene))

#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=wong_esc_genes_features, name="wong_esc")


#Subset for only Osi DTPCs
h1975_dtpc <- subset(h1975_combined, orig.ident == "H1975_Osi_DTPC")


###Wong Embryonic Stem Cell core violin plot (6/17/2025 update; no points)
H1975_krt17_dtpc_wong_esc_vp <- VlnPlot(h1975_dtpc, features = "signature_1wong_esc", group.by = "ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "H1975 Osi DTPC \nWong ESC")

H1975_krt17_dtpc_wong_esc_vp #Plot: 350 x 500



###EMT

#Read in Hallmark EMT genes
emt_genes <- Hallmark_EMT_genes

#Set Hallmark EMT gene features
emt_genes_features <- list(c(emt_genes$Gene))

#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=emt_genes_features, name="emt")


#Subset for only Osi DTPCs
h1975_dtpc <- subset(h1975_combined, orig.ident == "H1975_Osi_DTPC")


###Hallmark EMT genes violin plot
h1975_krt17_dtpc_emt_vp <- VlnPlot(h1975_dtpc, features = "signature_1emt", group.by = "ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "H1975 Osi DTPC \nEMT")

h1975_krt17_dtpc_emt_vp #Plot: 350 x 500



