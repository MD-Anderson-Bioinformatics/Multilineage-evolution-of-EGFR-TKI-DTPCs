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
library(scico)



##### HCC827 KRT17+ Osi DTPC reclustering #####

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


#Subset for only Osi DTPCs
hcc827_dtpc <- subset(hcc827_combined, orig.ident == "HCC827_Osi_DTPC")

#Subset for only KRT17+ osi DTPCs
hcc827_krt17_dtpc <- subset(hcc827_dtpc, ExpressionLabel == "KRT17+")

#Re-run SCTransform
hcc827_krt17_dtpc <- SCTransform(hcc827_krt17_dtpc, 
                                 vst.flavor = "v2", #specifying v2 for updated vst
                                 vars.to.regress = "percent.mt", 
                                 verbose = FALSE) #requires matrixStats < version 1.2


#Set SCT transformed/normalized values as VariableFeatures, required for PCA
VariableFeatures(hcc827_krt17_dtpc[["SCT"]]) <- rownames(hcc827_krt17_dtpc[["SCT"]]@scale.data)

#Run PCA
hcc827_krt17_dtpc <- RunPCA(hcc827_krt17_dtpc, verbose = FALSE)

#Run UMAP
hcc827_krt17_dtpc <- RunUMAP(hcc827_krt17_dtpc, reduction = "pca", dims = 1:30, verbose = FALSE)

#Find neighbors (SNN)
hcc827_krt17_dtpc <- FindNeighbors(hcc827_krt17_dtpc, reduction = "pca", dims = 1:30)

#Find clusters
hcc827_krt17_dtpc <- FindClusters(hcc827_krt17_dtpc, resolution = 0.8) #seurat default resolution


#Plot UMAP, false coloring by sample condition
hcc827_krt17_dtpc_reclust_umap <- DimPlot(hcc827_krt17_dtpc, reduction = "umap", group.by = "seurat_clusters",
                                          pt.size = 2) +
  theme(plot.title = element_text(color = "black", face = "bold", size = 12),
        legend.position = "right",
        legend.direction = "vertical",
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12)) +
  scale_y_continuous("UMAP2") +
  scale_x_continuous("UMAP1") +
  labs(title = "HCC827 KRT17+ Osi DTPC")

hcc827_krt17_dtpc_reclust_umap #450x350


###DE cluster markers gene analysis

#Set seurat cluster as identity
Idents(hcc827_krt17_dtpc) <- "seurat_clusters"

#Prep SCTFindMarkers
hcc827_krt17_dtpc %<>% PrepSCTFindMarkers(assay = "SCT")

#FindAllMarkers (requires normalizing library sizes first)
hcc827_krt17_dtpc_cluster_markers <- FindAllMarkers(hcc827_krt17_dtpc, assay = "SCT", recorrect_umi = FALSE)

#Check beginning of marker gene table
head(hcc827_krt17_dtpc_cluster_markers)

#Write marker gene table to file
write.table(hcc827_krt17_dtpc_cluster_markers, "C:/Users/bbmorris1/Desktop/HCC827 KRT17+ Osi DTPC reclustered seurat cluster marker genes.txt", sep = "\t", quote = FALSE)


#Filter for top 10 DE genes for each cluster
hcc827_krt17_dtpc_cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#Draw cluster marker heatmap
hcc827_krt17_dtpc_markers_h <- DoHeatmap(hcc827_krt17_dtpc, features = top10$gene, raster = FALSE) + 
  NoLegend() +
  theme(axis.text.y = element_text(color = "black", size = 4))

hcc827_krt17_dtpc_markers_h #Save also as PDF for gene name clarity

#Save as svg with ggsave for increased label clarity
ggsave(hcc827_krt17_dtpc_markers_h, file = "2025.6.18 HCC827 KRT17+ Osi DTPC reclustered cluster marker heatmap.svg", width = 4, height = 3)


##### HCC4006 KRT17+ Osi DTPC reclustering #####

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

#Subset for only KRT17+ osi DTPCs
hcc4006_krt17_dtpc <- subset(hcc4006_dtpc, ExpressionLabel == "KRT17+")

#Re-run SCTransform
hcc4006_krt17_dtpc <- SCTransform(hcc4006_krt17_dtpc, 
                                  vst.flavor = "v2", #specifying v2 for updated vst
                                  vars.to.regress = "percent.mt", 
                                  verbose = FALSE) #requires matrixStats < version 1.2


#Set SCT transformed/normalized values as VariableFeatures, required for PCA
VariableFeatures(hcc4006_krt17_dtpc[["SCT"]]) <- rownames(hcc4006_krt17_dtpc[["SCT"]]@scale.data)

#Run PCA
hcc4006_krt17_dtpc <- RunPCA(hcc4006_krt17_dtpc, verbose = FALSE)

#Run UMAP
hcc4006_krt17_dtpc <- RunUMAP(hcc4006_krt17_dtpc, reduction = "pca", dims = 1:30, verbose = FALSE)

#Find neighbors (SNN)
hcc4006_krt17_dtpc <- FindNeighbors(hcc4006_krt17_dtpc, reduction = "pca", dims = 1:30)

#Find clusters
hcc4006_krt17_dtpc <- FindClusters(hcc4006_krt17_dtpc, resolution = 0.8) #seurat default resolution


#Plot UMAP, false coloring by sample condition
hcc4006_krt17_dtpc_reclust_umap <- DimPlot(hcc4006_krt17_dtpc, reduction = "umap", group.by = "seurat_clusters",
                                           pt.size = 2) +
  theme(plot.title = element_text(color = "black", face = "bold", size = 12),
        legend.position = "right",
        legend.direction = "vertical",
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12)) +
  scale_y_continuous("UMAP2") +
  scale_x_continuous("UMAP1") +
  labs(title = "HCC4006 KRT17+ Osi DTPC")

hcc4006_krt17_dtpc_reclust_umap


###DE cluster markers gene analysis

#Set seurat cluster as identity
Idents(hcc4006_krt17_dtpc) <- "seurat_clusters"

#PrepSCTFindMarkers
hcc4006_krt17_dtpc %<>% PrepSCTFindMarkers(assay = "SCT")

#FindAllMarkers (requires normalizing library sizes first)
hcc4006_krt17_dtpc_cluster_markers <- FindAllMarkers(hcc4006_krt17_dtpc, assay = "SCT", recorrect_umi = FALSE)

#Check beginning of marker gene table
head(hcc4006_krt17_dtpc_cluster_markers)

#Write marker gene table to file
write.table(hcc4006_krt17_dtpc_cluster_markers, "C:/Users/bbmorris1/Desktop/HCC4006 KRT17+ Osi DTPC reclustered seurat cluster marker genes.txt", sep = "\t", quote = FALSE)


#Filter for top 10 genes for each cluster
hcc4006_krt17_dtpc_cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#Draw cluster marker heatmap
hcc4006_krt17_dtpc_markers_h <- DoHeatmap(hcc4006_krt17_dtpc, features = top10$gene, raster = FALSE) + 
  NoLegend() +
  theme(axis.text.y = element_text(color = "black", size = 4))

hcc4006_krt17_dtpc_markers_h #Save also as PDF for gene name clarity


##### H1975 KRT17+ Osi DTPC reclustering #####

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

#Subset for only osi DTPCs
h1975_dtpc <- subset(h1975_combined, orig.ident == "H1975_Osi_DTPC")

#Subset for only KRT17+ osi DTPCs
h1975_krt17_dtpc <- subset(h1975_dtpc, ExpressionLabel == "KRT17+")


#Re-run SCTransform
h1975_krt17_dtpc <- SCTransform(h1975_krt17_dtpc, 
                                vst.flavor = "v2", #specifying v2 for updated vst
                                vars.to.regress = "percent.mt", 
                                verbose = FALSE) #requires matrixStats < version 1.2


#Set SCT transformed/normalized values as VariableFeatures, required for PCA
VariableFeatures(h1975_krt17_dtpc[["SCT"]]) <- rownames(h1975_krt17_dtpc[["SCT"]]@scale.data)

#Run PCA
h1975_krt17_dtpc <- RunPCA(h1975_krt17_dtpc, verbose = FALSE)

#Run UMAP
h1975_krt17_dtpc <- RunUMAP(h1975_krt17_dtpc, reduction = "pca", dims = 1:30, verbose = FALSE)

#Find neighbors (SNN)
h1975_krt17_dtpc <- FindNeighbors(h1975_krt17_dtpc, reduction = "pca", dims = 1:30)

#Find clusters
h1975_krt17_dtpc <- FindClusters(h1975_krt17_dtpc, resolution = 0.8) #seurat default resolution


#Plot UMAP, false coloring by sample condition
h1975_krt17_dtpc_reclust_umap <- DimPlot(h1975_krt17_dtpc, reduction = "umap", group.by = "seurat_clusters",
                                         pt.size = 2) +
  theme(plot.title = element_text(color = "black", face = "bold", size = 12),
        legend.position = "right",
        legend.direction = "vertical",
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12)) +
  scale_y_continuous("UMAP2") +
  scale_x_continuous("UMAP1") +
  labs(title = "H1975 KRT17+ Osi DTPC")

h1975_krt17_dtpc_reclust_umap


###DE cluster markers gene analysis

#Set seurat cluster as identity
Idents(h1975_krt17_dtpc) <- "seurat_clusters"

#PrepSCTFindMarkers
h1975_krt17_dtpc %<>% PrepSCTFindMarkers(assay = "SCT")

#FindAllMarkers (requires normalizing library sizes first)
h1975_krt17_dtpc_cluster_markers <- FindAllMarkers(h1975_krt17_dtpc, assay = "SCT", recorrect_umi = FALSE)

#Check beginning of marker gene table
head(h1975_krt17_dtpc_cluster_markers)

#Write marker gene table to file
write.table(h1975_krt17_dtpc_cluster_markers, "C:/Users/bbmorris1/Desktop/H1975 KRT17+ Osi DTPC reclustered seurat cluster marker genes.txt", sep = "\t", quote = FALSE)


#Filter for top 10 genes for clusters
h1975_krt17_dtpc_cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#Draw cluster marker heatmap
h1975_krt17_dtpc_markers_h <- DoHeatmap(h1975_krt17_dtpc, features = top10$gene, raster = FALSE) + 
  NoLegend() +
  theme(axis.text.y = element_text(color = "black", size = 4))


h1975_krt17_dtpc_markers_h #Save also as PDF for gene name clarity






