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
library(sctransform)
library(ggpubr)
library(forcats)



##### HCC827 KRT17+ Osi DTPC lung lineage marker co-expression analysis scoring analysis #####

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


#Subset for only Osi DTPCs
hcc827_dtpc <- subset(hcc827_combined, orig.ident == "HCC827_Osi_DTPC")


###NKX2-1 violin plot
hcc827_krt17_dtpc_nkx2_1_vp <- VlnPlot(hcc827_dtpc, features = "NKX2-1", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "HCC827 Osi DTPC \nNKX2-1") 

hcc827_krt17_dtpc_nkx2_1_vp


###HOPX violin plot
hcc827_krt17_dtpc_hopx_vp <- VlnPlot(hcc827_dtpc, features = "HOPX", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "HCC827 Osi DTPC \nHOPX") 

hcc827_krt17_dtpc_hopx_vp


###TP63 violin plot
hcc827_krt17_dtpc_tp63_vp <- VlnPlot(hcc827_dtpc, features = "TP63", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "HCC827 Osi DTPC \nTP63") 

hcc827_krt17_dtpc_tp63_vp


###SOX4 violin plot
hcc827_sox4_dtpc_tp63_vp <- VlnPlot(hcc827_dtpc, features = "SOX4", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "HCC827 Osi DTPC \nSOX4") 

hcc827_sox4_dtpc_tp63_vp



##### HCC4006 KRT17+ Osi DTPC Osi DTPC lung lineage marker co-expression analysis scoring analysis #####

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



#Subset for only Osi DTPCs
hcc4006_dtpc <- subset(hcc4006_combined, orig.ident == "HCC4006_Osi_DTPC")


###NKX2-1 violin plot
hcc4006_krt17_dtpc_nkx2_1_vp <- VlnPlot(hcc4006_dtpc, features = "NKX2-1", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "HCC4006 Osi DTPC \nNKX2-1") 

hcc4006_krt17_dtpc_nkx2_1_vp


###HOPX violin plot
hcc4006_krt17_dtpc_hopx_vp <- VlnPlot(hcc4006_dtpc, features = "HOPX", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "HCC4006 Osi DTPC \nHOPX") 

hcc4006_krt17_dtpc_hopx_vp


###TP63 violin plot
hcc4006_krt17_dtpc_tp63_vp <- VlnPlot(hcc4006_dtpc, features = "TP63", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "HCC4006 Osi DTPC \nTP63") 

hcc4006_krt17_dtpc_tp63_vp


###SOX4 violin plot
hcc4006_krt17_dtpc_sox4_vp <- VlnPlot(hcc4006_dtpc, features = "SOX4", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "HCC4006 Osi DTPC \nSOX4") 

hcc4006_krt17_dtpc_sox4_vp



##### H1975 KRT17+ Osi DTPC Osi DTPC lung lineage marker co-expression analysis scoring analysis #####

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


#Subset for only Osi DTPCs
h1975_dtpc <- subset(h1975_combined, orig.ident == "H1975_Osi_DTPC")


###NKX2-1 violin plot
h1975_krt17_dtpc_nkx2_1_vp <- VlnPlot(h1975_dtpc, features = "NKX2-1", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "H1975 Osi DTPC \nNKX2-1") 

h1975_krt17_dtpc_nkx2_1_vp


###HOPX violin plot
h1975_krt17_dtpc_hopx_vp <- VlnPlot(h1975_dtpc, features = "HOPX", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "H1975 Osi DTPC \nHOPX") 

h1975_krt17_dtpc_hopx_vp


###TP63 violin plot
h1975_krt17_dtpc_tp63_vp <- VlnPlot(h1975_dtpc, features = "TP63", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "H1975 Osi DTPC \nTP63") 

h1975_krt17_dtpc_tp63_vp



###SOX4 violin plot
h1975_krt17_dtpc_sox4_vp <- VlnPlot(h1975_dtpc, features = "SOX4", group.by = "ExpressionLabel") +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Expression") +
  labs(title = "H1975 Osi DTPC \nSOX4") 

h1975_krt17_dtpc_sox4_vp







