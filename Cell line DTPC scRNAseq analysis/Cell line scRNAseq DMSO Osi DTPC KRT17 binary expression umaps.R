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
library(forcats)



##### HCC827 DMSO & Osi DTPC KRT17 binary expression (KRT17+ vs KRT17-) dim plot #####

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

#Define custom color mapping (required for proper visualization)
colors <- c("KRT17-" = "gray",
            "KRT17+" = "red")

#Plot UMAP with KRT17 pseudocoloring
hcc827_krt17_dim <- DimPlot(hcc827_combined, group.by = "ExpressionLabel", 
                            split.by = "orig.ident",
                            pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT17 expression",
       color = "KRT17") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

hcc827_krt17_dim



##### HCC4006 DMSO & Osi DTPC KRT17 binary expression (KRT17+ vs KRT17-) dim plot #####

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

#Define custom color mapping (required for proper visualization)
colors <- c("KRT17-" = "gray",
            "KRT17+" = "red")

#Plot UMAP with KRT17 pseudocoloring
hcc4006_krt17_dim <- DimPlot(hcc4006_combined, group.by = "ExpressionLabel", 
                             split.by = "orig.ident",
                             pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT17 expression",
       color = "KRT17") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank())

hcc4006_krt17_dim



##### H1975 DMSO & Osi DTPC KRT17 binary expression (KRT17+ vs KRT17-) dim plot #####

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

#Define custom color mapping (required for proper visualization)
colors <- c("KRT17-" = "gray",
            "KRT17+" = "red")

#Plot UMAP with KRT17 pseudocoloring
h1975_krt17_dim <- DimPlot(h1975_combined, group.by = "ExpressionLabel", 
                           split.by = "orig.ident",
                           pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT17 expression",
       color = "KRT17") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank())

h1975_krt17_dim








