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
library(scplotter)
library(ggpubr)



##### HCC827 DMSO and HCC827 Osi DTPC KRT17 MAPK activity analysis #####

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###MAPK pathway activity scoring

#Read in mapk pathway activity score
mapk_act_genes <- MAPK_pathway_activity_score

#Set mapk pathway activity score features
mapk_act_genes_features <- list(c(mapk_act_genes$Gene))


#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=mapk_act_genes_features, name="MAPK_act")

###MAPK Pathway Activity Score violin plot
hcc827_mapk_vp <- VlnPlot(hcc827_combined, features = "signature_1MAPK_act", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "MAPK Pathway \nActivity Score")

hcc827_mapk_vp #Plot: 350 x 500


###KRT17 labeling

#Define genes of interest
genes_of_interest <- c("KRT17")

#Extract expression data for the three genes
KRT17_expression_data <- FetchData(hcc827_combined, vars = genes_of_interest)

#Convert to dataframe
KRT17_expression_data <- as.data.frame(KRT17_expression_data)

#Define thresholds for gene expression (adjust as needed)
threshold <- 0.1
KRT17_expression_data <- KRT17_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT17 > threshold, 1, 0),
    KRT17_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc827_combined <- AddMetaData(hcc827_combined, KRT17_expression_data$KRT17_ExpressionCategory, col.name = "KRT17_ExpressionCategory")

#Build custom color key for KRT17_ExpressionCategory
existing_column <- "KRT17_ExpressionCategory"

label_map <- c(
  "0" = "KRT17-",
  "1" = "KRT17+")


#Assign KRT17_ExpressionCategory data to separate metadata variable
metadata <- hcc827_combined@meta.data

#Match surface target co-expression labels to ExpressionCategory binary values using key
metadata$KRT17_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT17 labels back to metadata
hcc827_combined <- AddMetaData(hcc827_combined, metadata$KRT17_ExpressionLabel, col.name = "KRT17_ExpressionLabel")

#Check new metadata entry
head(hcc827_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT17-" = "gray",
            "KRT17+" = "red")

#Plot UMAP with KRT17 pseudocoloring
hcc827_KRT17_dim <- DimPlot(hcc827_combined, group.by = "KRT17_ExpressionLabel", 
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

hcc827_KRT17_dim #Plot: 800 x 400



###MAPK Pathway Activity Score violin plot
hcc827_krt17_mapk_vp <- VlnPlot(hcc827_combined, features = "signature_1MAPK_act", group.by = "orig.ident", split.by = "KRT17_ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "right") +
  scale_y_continuous("Activity") +
  labs(title = "KRT17 expression MAPK Pathway \nActivity Score")

hcc827_krt17_mapk_vp #Plot: 500 x 500



##### HCC4006 DMSO and HCC4006 Osi DTPC KRT17 MAPK activity analysis #####

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###MAPK pathway activity scoring

#Read in mapk pathway activity score
mapk_act_genes <- MAPK_pathway_activity_score

#Set mapk pathway activity score features
mapk_act_genes_features <- list(c(mapk_act_genes$Gene))


#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=mapk_act_genes_features, name="MAPK_act")

###MAPK Pathway Activity Score violin plot
hcc4006_mapk_vp <- VlnPlot(hcc4006_combined, features = "signature_1MAPK_act", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "MAPK Pathway \nActivity Score")

hcc4006_mapk_vp #Plot: 350 x 500


###KRT17 labeling

#Define genes of interest
genes_of_interest <- c("KRT17")

#Extract expression data for the three genes
KRT17_expression_data <- FetchData(hcc4006_combined, vars = genes_of_interest)

#Convert to dataframe
KRT17_expression_data <- as.data.frame(KRT17_expression_data)

#Define thresholds for gene expression (adjust as needed)
threshold <- 0.1
KRT17_expression_data <- KRT17_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT17 > threshold, 1, 0),
    KRT17_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc4006_combined <- AddMetaData(hcc4006_combined, KRT17_expression_data$KRT17_ExpressionCategory, col.name = "KRT17_ExpressionCategory")

#Build custom color key for KRT17_ExpressionCategory
existing_column <- "KRT17_ExpressionCategory"

label_map <- c(
  "0" = "KRT17-",
  "1" = "KRT17+")


#Assign KRT17_ExpressionCategory data to separate metadata variable
metadata <- hcc4006_combined@meta.data

#Match surface target co-expression labels to ExpressionCategory binary values using key
metadata$KRT17_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT17 labels back to metadata
hcc4006_combined <- AddMetaData(hcc4006_combined, metadata$KRT17_ExpressionLabel, col.name = "KRT17_ExpressionLabel")

#Check new metadata entry
head(hcc4006_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT17-" = "gray",
            "KRT17+" = "red")

#Plot UMAP with KRT17 pseudocoloring
hcc4006_KRT17_dim <- DimPlot(hcc4006_combined, group.by = "KRT17_ExpressionLabel", 
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

hcc4006_KRT17_dim #Plot: 800 x 400


###MAPK Pathway Activity Score violin plot
hcc4006_krt17_mapk_vp <- VlnPlot(hcc4006_combined, features = "signature_1MAPK_act", group.by = "orig.ident", split.by = "KRT17_ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "right") +
  scale_y_continuous("Activity") +
  labs(title = "KRT17 expression MAPK Pathway \nActivity Score")

hcc4006_krt17_mapk_vp #Plot: 700 x 500



##### H1975 DMSO and H1975 Osi DTPC integrated KRT17 MAPK activity analysis #####

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###MAPK pathway activity scoring

#Read in mapk pathway activity score
mapk_act_genes <- MAPK_pathway_activity_score

#Set mapk pathway activity score features
mapk_act_genes_features <- list(c(mapk_act_genes$Gene))


#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=mapk_act_genes_features, name="MAPK_act")


###MAPK Pathway Activity Score violin plot 
h1975_mapk_vp <- VlnPlot(h1975_combined, features = "signature_1MAPK_act", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "MAPK Pathway \nActivity Score")

h1975_mapk_vp #Plot: 350 x 500


###KRT17 labeling

#Define genes of interest
genes_of_interest <- c("KRT17")

#Extract expression data for the three genes
KRT17_expression_data <- FetchData(h1975_combined, vars = genes_of_interest)

#Convert to dataframe
KRT17_expression_data <- as.data.frame(KRT17_expression_data)

#Define thresholds for gene expression (adjust as needed)
threshold <- 0.1
KRT17_expression_data <- KRT17_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT17 > threshold, 1, 0),
    KRT17_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
h1975_combined <- AddMetaData(h1975_combined, KRT17_expression_data$KRT17_ExpressionCategory, col.name = "KRT17_ExpressionCategory")

#Build custom color key for KRT17_ExpressionCategory
existing_column <- "KRT17_ExpressionCategory"

label_map <- c(
  "0" = "KRT17-",
  "1" = "KRT17+")


#Assign KRT17_ExpressionCategory data to separate metadata variable
metadata <- h1975_combined@meta.data

#Match surface target co-expression labels to ExpressionCategory binary values using key
metadata$KRT17_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT17 labels back to metadata
h1975_combined <- AddMetaData(h1975_combined, metadata$KRT17_ExpressionLabel, col.name = "KRT17_ExpressionLabel")

#Check new metadata entry
head(h1975_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT17-" = "gray",
            "KRT17+" = "red")

#Plot UMAP with KRT17 pseudocoloring
h1975_KRT17_dim <- DimPlot(h1975_combined, group.by = "KRT17_ExpressionLabel", 
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

h1975_KRT17_dim #Plot: 800 x 400


###MAPK Pathway Activity Score violin plot 
h1975_krt17_mapk_vp <- VlnPlot(h1975_combined, features = "signature_1MAPK_act", group.by = "orig.ident", split.by = "KRT17_ExpressionLabel", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "right") +
  scale_y_continuous("Activity") +
  labs(title = "KRT17 expression MAPK Pathway \nActivity Score")

h1975_krt17_mapk_vp #Plot: 700 x 500






