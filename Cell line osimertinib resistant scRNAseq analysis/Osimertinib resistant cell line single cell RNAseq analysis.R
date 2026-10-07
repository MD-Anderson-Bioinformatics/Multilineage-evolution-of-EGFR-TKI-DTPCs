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
library(devtools)
library(forcats)
library(scCustomize)
library(glmGamPoi)
library(patchwork)
library(magrittr)
library(ggrepel)



##### HCC4006 OR2 single cell preprocessing #####

#Read in h5 dataset
hcc4006_or2_h5 <- Read10X_h5("/EGFR_TKI_resistant_scRNAseq/HCC4006/OR2/cellbender_10x_pbmc_filtered.h5", use.names = TRUE, unique.features = TRUE)

#Create SeuratObject
hcc4006_or2 <- CreateSeuratObject(counts = hcc4006_or2_h5,
                                  project = "HCC4006_OR2")



#Calculate percentage of mitochondrial reads
hcc4006_or2[["percent.mt"]] <- PercentageFeatureSet(hcc4006_or2, pattern = "^MT-")

#Check meta data
head(hcc4006_or2@meta.data, 5)

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (before QC filtering)
unfilt_qc_p <- VlnPlot(hcc4006_or2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfilt_qc_p

#Visualize number of RNA counts vs %MT reads
nRNA_mt_p <- FeatureScatter(hcc4006_or2, feature1 = "nCount_RNA", feature2 = "percent.mt")
nRNA_mt_p

#Visualize number of RNA counts vs number of unique features
nRNA_nFeat_p <- FeatureScatter(hcc4006_or2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nRNA_nFeat_p

#Remove cells with < 200 genes detected and/or >20% mitochondrial reads
hcc4006_or2_1 <- subset(hcc4006_or2, subset = nFeature_RNA > 200 & percent.mt < 20) #removed 300 cells

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (after QC filtering)
filt_qc_p <- VlnPlot(hcc4006_or2_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filt_qc_p


###SCTransform

#Run scTransform, regressing out percent mitochondrial reads
hcc4006_or2_1 <- SCTransform(hcc4006_or2_1, 
                            vst.flavor = "v2", #specifying v2 for updated vst
                            vars.to.regress = "percent.mt", 
                            verbose = FALSE) #requires matrixStats < version 1.2


###PCA

#Run PCA
hcc4006_or2_1 <- RunPCA(hcc4006_or2_1, verbose = FALSE)

#Identify number of PCs that define >80% of variance
pcs = hcc4006_or2_1@reductions$pca@cell.embeddings
pca_var = hcc4006_or2_1@reductions$pca@stdev ** 2
pca_var_cum = sapply(1:length(pca_var), function(x) sum(pca_var[1:x])/sum(pca_var))
npca = which(pca_var_cum>0.80)[1]
npca 


###Find Neighbors

#Find nearest neighbors
hcc4006_or2_1 <- FindNeighbors(hcc4006_or2_1, dims = 1:npca, verbose = FALSE)


##Find Clusters

#Generate clusters
hcc4006_or2_1 <- FindClusters(hcc4006_or2_1, verbose = FALSE)


###UMAP clustering

#Run UMAP
hcc4006_or2_1 <- RunUMAP(hcc4006_or2_1, dims = 1:npca)

#Draw UMAP
DimPlot(hcc4006_or2_1, label = TRUE, reduction = "umap")



###Doublet removal with DoubletFinder

#PK identification
sweep.res <- paramSweep(hcc4006_or2_1, PCs = 1:npca, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #Plot pk (x) vs BCmvn (y)

#Note: use 19 for doublet identification (value in [] in code line below);
#value changes based on maximum of BCmvn plot above

#Use PK determined from above plot
pK.set <- unique(sweep.stats$pK)[19]

#Doublet proportion estimation
nExp_poi <- round(0.08*nrow(hcc4006_or2_1@meta.data))

#Run DoubletFinder using above metrics
hcc4006_or2_1 <- doubletFinder(hcc4006_or2_1, PCs = 1:npca, pN = 0.25, 
                              pK = as.numeric(as.character(pK.set)), 
                              nExp = nExp_poi, reuse.pANN = NULL, 
                              sct = TRUE)

#Subset data for only cells likely as singlets
hcc4006_or2_1 <- subset(hcc4006_or2_1,  
                        DF.classifications_0.25_0.18_313 == "Singlet")


###Leaves 3596 cells (removed  cells)

#Save SeuratRDS (post doublet removal)
SaveSeuratRds(hcc4006_or2_1, "HCC4006_OR2_scTransform_v2_singlets.Rds")


FeaturePlot(hcc4006_or2_1, features = "KRT17") #Note: ~40% cells express KRT17



##### HCC4006 OR7 single cell preprocessing #####

#Read in h5 dataset
hcc4006_or7_h5 <- Read10X_h5("/EGFR_TKI_resistant_scRNAseq/HCC4006/OR7/cellbender_10x_pbmc_filtered.h5", use.names = TRUE, unique.features = TRUE)

#Create SeuratObject
hcc4006_or7 <- CreateSeuratObject(counts = hcc4006_or7_h5,
                                  project = "HCC4006_OR7")


#Calculate percentage of mitochondrial reads
hcc4006_or7[["percent.mt"]] <- PercentageFeatureSet(hcc4006_or7, pattern = "^MT-")

#Check meta data
head(hcc4006_or7@meta.data, 5)

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (before QC filtering)
unfilt_qc_p <- VlnPlot(hcc4006_or7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfilt_qc_p

#Visualize number of RNA counts vs %MT reads
nRNA_mt_p <- FeatureScatter(hcc4006_or7, feature1 = "nCount_RNA", feature2 = "percent.mt")
nRNA_mt_p

#Visualize number of RNA counts vs number of unique features
nRNA_nFeat_p <- FeatureScatter(hcc4006_or7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nRNA_nFeat_p

#Remove cells with < 200 genes detected and/or >20% mitochondrial reads
hcc4006_or7_1 <- subset(hcc4006_or7, subset = nFeature_RNA > 200 & percent.mt < 20) #removed 162 cells

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (after QC filtering)
filt_qc_p <- VlnPlot(hcc4006_or7_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filt_qc_p


###SCTransform

#Run scTransform, regressing out percent mitochondrial reads
hcc4006_or7_1 <- SCTransform(hcc4006_or7_1, 
                             vst.flavor = "v2", #specifying v2 for updated vst
                             vars.to.regress = "percent.mt", 
                             verbose = FALSE) #requires matrixStats < version 1.2


###PCA

#Run PCA
hcc4006_or7_1 <- RunPCA(hcc4006_or7_1, verbose = FALSE)

#Identify number of PCs that define >80% of variance
pcs = hcc4006_or7_1@reductions$pca@cell.embeddings
pca_var = hcc4006_or7_1@reductions$pca@stdev ** 2
pca_var_cum = sapply(1:length(pca_var), function(x) sum(pca_var[1:x])/sum(pca_var))
npca = which(pca_var_cum>0.80)[1]
npca 


###Find Neighbors

#Find nearest neighbors
hcc4006_or7_1 <- FindNeighbors(hcc4006_or7_1, dims = 1:npca, verbose = FALSE)


##Find Clusters

#Generate clusters
hcc4006_or7_1 <- FindClusters(hcc4006_or7_1, verbose = FALSE)


###UMAP clustering

#Run UMAP
hcc4006_or7_1 <- RunUMAP(hcc4006_or7_1, dims = 1:npca)

#Draw UMAP
DimPlot(hcc4006_or7_1, label = TRUE, reduction = "umap")


###Doublet removal with DoubletFinder

#PK identification
sweep.res <- paramSweep(hcc4006_or7_1, PCs = 1:npca, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #Plot pk (x) vs BCmvn (y)

#Note: use 11 for doublet identification (value in [] in code line below);
#value changes based on maximum of BCmvn plot above

#Use PK determined from above plot
pK.set <- unique(sweep.stats$pK)[11]

#Doublet proportion estimation
nExp_poi <- round(0.08*nrow(hcc4006_or7_1@meta.data))

#Run DoubletFinder using above metrics
hcc4006_or7_1 <- doubletFinder(hcc4006_or7_1, PCs = 1:npca, pN = 0.25, 
                               pK = as.numeric(as.character(pK.set)), 
                               nExp = nExp_poi, reuse.pANN = NULL, 
                               sct = TRUE)

#Subset data for only cells likely as singlets
hcc4006_or7_1 <- subset(hcc4006_or7_1,  
                        DF.classifications_0.25_0.1_249 == "Singlet")


###Leaves 2859 cells (removed 411 cells)

#Save SeuratRDS (post doublet removal)
SaveSeuratRds(hcc4006_or7_1, "HCC4006_OR7_scTransform_v2_singlets.Rds")


FeaturePlot(hcc4006_or7_1, features = "KRT17") #Note: subset of cells express KRT17




##### HCC4006 Parental single cell preprocessing #####

#Read in h5 dataset
hcc4006_parental_h5 <- Read10X_h5("/EGFR_TKI_resistant_scRNAseq/HCC4006/Parental/cellbender_10x_pbmc_filtered.h5", use.names = TRUE, unique.features = TRUE)

#Create SeuratObject
hcc4006_parental <- CreateSeuratObject(counts = hcc4006_parental_h5,
                                  project = "HCC4006_Parental")



#Calculate percentage of mitochondrial reads
hcc4006_parental[["percent.mt"]] <- PercentageFeatureSet(hcc4006_parental, pattern = "^MT-")

#Check meta data
head(hcc4006_parental@meta.data, 5)

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (before QC filtering)
unfilt_qc_p <- VlnPlot(hcc4006_parental, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfilt_qc_p

#Visualize number of RNA counts vs %MT reads
nRNA_mt_p <- FeatureScatter(hcc4006_parental, feature1 = "nCount_RNA", feature2 = "percent.mt")
nRNA_mt_p

#Visualize number of RNA counts vs number of unique features
nRNA_nFeat_p <- FeatureScatter(hcc4006_parental, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nRNA_nFeat_p

#Remove cells with < 200 genes detected and/or >20% mitochondrial reads
hcc4006_parental_1 <- subset(hcc4006_parental, subset = nFeature_RNA > 200 & percent.mt < 20) #removed 315 cells

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (after QC filtering)
filt_qc_p <- VlnPlot(hcc4006_parental_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filt_qc_p


###SCTransform

#Run scTransform, regressing out percent mitochondrial reads
hcc4006_parental_1 <- SCTransform(hcc4006_parental_1, 
                             vst.flavor = "v2", #specifying v2 for updated vst
                             vars.to.regress = "percent.mt", 
                             verbose = FALSE) #requires matrixStats < version 1.2


###PCA

#Run PCA
hcc4006_parental_1 <- RunPCA(hcc4006_parental_1, verbose = FALSE)

#Identify number of PCs that define >80% of variance
pcs = hcc4006_parental_1@reductions$pca@cell.embeddings
pca_var = hcc4006_parental_1@reductions$pca@stdev ** 2
pca_var_cum = sapply(1:length(pca_var), function(x) sum(pca_var[1:x])/sum(pca_var))
npca = which(pca_var_cum>0.80)[1]
npca 


###Find Neighbors

#Find nearest neighbors
hcc4006_parental_1 <- FindNeighbors(hcc4006_parental_1, dims = 1:npca, verbose = FALSE)


##Find Clusters

#Generate clusters
hcc4006_parental_1 <- FindClusters(hcc4006_parental_1, verbose = FALSE)


###UMAP clustering

#Run UMAP
hcc4006_parental_1 <- RunUMAP(hcc4006_parental_1, dims = 1:npca)

#Draw UMAP
DimPlot(hcc4006_parental_1, label = TRUE, reduction = "umap")


###Doublet removal with DoubletFinder

#PK identification
sweep.res <- paramSweep(hcc4006_parental_1, PCs = 1:npca, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #Plot pk (x) vs BCmvn (y)

#Note: use 9 for doublet identification (value in [] in code line below);
#value changes based on maximum of BCmvn plot above

#Use PK determined from above plot
pK.set <- unique(sweep.stats$pK)[9]

#Doublet proportion estimation
nExp_poi <- round(0.08*nrow(hcc4006_parental_1@meta.data))

#Run DoubletFinder using above metrics
hcc4006_parental_1 <- doubletFinder(hcc4006_parental_1, PCs = 1:npca, pN = 0.25, 
                               pK = as.numeric(as.character(pK.set)), 
                               nExp = nExp_poi, reuse.pANN = NULL, 
                               sct = TRUE)

#Subset data for only cells likely as singlets
hcc4006_parental_1 <- subset(hcc4006_parental_1,  
                             DF.classifications_0.25_0.08_231  == "Singlet")


###Leaves 2660 cells (removed 546 cells)

#Save SeuratRDS (post doublet removal)
SaveSeuratRds(hcc4006_parental_1, "HCC4006_Parental_scTransform_v2_singlets.Rds")


FeaturePlot(hcc4006_parental_1, features = "KRT17") #Note: 2 cells express KRT17



##### HCC4006 Parental, OR2, and OR7 integrated analysis #####

#Load data
hcc4006_parental_1 <- LoadSeuratRds("HCC4006_Parental_scTransform_v2_singlets.Rds")
hcc4006_or2_1 <- LoadSeuratRds("HCC4006_OR2_scTransform_v2_singlets.Rds")
hcc4006_or7_1 <- LoadSeuratRds("HCC4006_OR7_scTransform_v2_singlets.Rds")

#Merge two Seurat objects
hcc4006_parental_or2_combined <- merge(x = hcc4006_parental_1, y = hcc4006_or2_1, add.cell.ids = c("HCC4006_Parental", "HCC4006_OR2"), project = "HCC4006", merge.data = TRUE)
hcc4006_parental_or2_combined

#Merge Seurat objects
hcc4006_parental_or2_or7_combined <- merge(x = hcc4006_parental_or2_combined, y = hcc4006_or7_1, add.cell.id2 = c("HCC4006_OR7"), project = "HCC4006", merge.data = TRUE)
hcc4006_parental_or2_or7_combined


#Set SCT transformed/normalized values as VariableFeatures, required for PCA
VariableFeatures(hcc4006_parental_or2_or7_combined[["SCT"]]) <- rownames(hcc4006_parental_or2_or7_combined[["SCT"]]@scale.data)


#Note: setting SCT data as variable features is required for SelectIntegrationFeatures & PrepSCTIntegration

#Split object
hcc4006_cond_list <- SplitObject(hcc4006_parental_or2_or7_combined, split.by = "orig.ident")

parental <- hcc4006_cond_list[["HCC4006_Parental"]]
or2 <- hcc4006_cond_list[["HCC4006_OR2"]]
or7 <- hcc4006_cond_list[["HCC4006_OR7"]]

#Select integration features
features <- SelectIntegrationFeatures(object.list = hcc4006_cond_list, nfeatures = 3000)
hcc4006_cond_list <- PrepSCTIntegration(object.list = hcc4006_cond_list, anchor.features = features)


#Find integration anchors
hcc4006_or_anchors <- FindIntegrationAnchors(object.list = hcc4006_cond_list, normalization.method = "SCT",
                                         anchor.features = features)


#Integrate data
hcc4006_or_combined_sct <- IntegrateData(anchorset = hcc4006_or_anchors, normalization.method = "SCT")


#Run PCA
hcc4006_or_combined_sct <- RunPCA(hcc4006_or_combined_sct, verbose = FALSE)

#Run UMAP
hcc4006_or_combined_sct <- RunUMAP(hcc4006_or_combined_sct, reduction = "pca", dims = 1:30, verbose = FALSE)

#Find neighbors (SNN)
hcc4006_or_combined_sct <- FindNeighbors(hcc4006_or_combined_sct, reduction = "pca", dims = 1:30)

#Find clusters
hcc4006_or_combined_sct <- FindClusters(hcc4006_or_combined_sct, resolution = 1)


#Plot UMAP, false coloring by sample condition
DimPlot(hcc4006_or_combined_sct, reduction = "umap", group.by = "orig.ident")

#Plot UMAP, false coloring by seurat cluster
seurat_cluster_umap <- DimPlot(hcc4006_or_combined_sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                               repel = TRUE)

seurat_cluster_umap


#Save SeurateRDS
SaveSeuratRds(hcc4006_or_combined_sct, "HCC4006_Parental_OR2_OR7_scTransformed_v2_singlet_CCAintegrated_annotated.Rds")



###DE marker gene analysis (Parental vs OR2)
Idents(hcc4006_or_combined_sct) <- "orig.ident"

hcc4006_or_combined_sct %<>% PrepSCTFindMarkers(assay = "SCT")

hcc4006_parental_vs_OR2 <- FindMarkers(hcc4006_or_combined_sct, assay = "SCT",
                                           recorrect_umi = FALSE,
                                           ident.1 = "HCC4006_Parental", 
                                           ident.2 = "HCC4006_OR2",
                                           verbose = FALSE)

#Write data to file
write.table(hcc4006_parental_vs_OR2, "HCC4006_parental_vs_OR2_DE_genes.txt", sep = "\t", quote = FALSE)


###DE marker gene analysis (Parental vs OR7)
Idents(hcc4006_or_combined_sct) <- "orig.ident"

hcc4006_or_combined_sct %<>% PrepSCTFindMarkers(assay = "SCT")

hcc4006_parental_vs_OR7 <- FindMarkers(hcc4006_or_combined_sct, assay = "SCT",
                                       recorrect_umi = FALSE,
                                       ident.1 = "HCC4006_Parental", 
                                       ident.2 = "HCC4006_OR7",
                                       verbose = FALSE)


#Write data to file
write.table(hcc4006_parental_vs_OR7, "HCC4006_parental_vs_OR7_DE_genes.txt", sep = "\t", quote = FALSE)



###Percentage of cells expressing KRT17 (HCC4006 Parental, OR2, OR7)

#Create dataframe
hcc4006_p_or_krt17 <- data.frame(Sample_ID = c("Parental", "OR2", "OR7"),
                                 Percent_KRT17 = c(0.5, 24, 5.4))

#Set levels
hcc4006_p_or_krt17$Sample_ID <- factor(hcc4006_p_or_krt17$Sample_ID, levels = c("Parental", "OR2", "OR7"))

#Draw barchart
hcc4006_p_or_krt17_bp <- ggplot(hcc4006_p_or_krt17, aes(x = Sample_ID, y = Percent_KRT17, fill = Sample_ID)) +
                         geom_bar(stat = "identity") +
                         geom_text(aes(label = Percent_KRT17), vjust = -0.5) +
                         theme_bw() +
                         theme(panel.grid = element_blank(),
                               axis.text = element_text(size = 12, color = "black"),
                               axis.title.x = element_blank(),
                               legend.position = "none") +
                         scale_y_continuous("Percent KRT17+ cells", expand = c(0, 0), limits = c(0, 30)) +
                         scale_fill_manual(values = c("Parental" = "lightgreen",
                                                      "OR2" = "springgreen4",
                                                      "OR7" = "springgreen4")) +
                         labs(title = "HCC4006 Parental vs Osi Resistant \nKRT17 expression")


hcc4006_p_or_krt17_bp #Plot: 350 x 500


##### HCC4006 Parental vs OR2 DE gene volcano plot #####

#Read in data
hcc4006_p_or2_de <- HCC4006_parental_vs_OR2_DE_genes_volcano_plot

###KRT17 volcano plot

#Identify genes to label
genes_to_label <- c("KRT17")

#Set labels if in DE list
hcc4006_p_or2_de$label <- ifelse(hcc4006_p_or2_de$Gene %in% genes_to_label, hcc4006_p_or2_de$Gene, NA)

#Draw labeled volcano plot for ETV5, HOPX, and SOX4
hcc4006_p_or2_de_vp <- ggplot(hcc4006_p_or2_de, aes(x = avg_log2FC_OR2_v_Parental, y = Minus_log10_padj)) +
  geom_point(aes(color = avg_log2FC_OR2_v_Parental > 0)) +
  scale_color_manual(values = c("FALSE" = "blue",
                                "TRUE" = "red")) +
  geom_point(data = subset(hcc4006_p_or2_de, label != 0),
             shape = 21,
             fill = NA,
             stroke = 1.2,
             color = "black") +
  geom_hline(yintercept = 1.30102999566398, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  geom_vline(xintercept = -1, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  geom_label_repel(aes(label = label),
                   fill = "white",
                   color = "black",
                   box.padding = 1,
                   point.padding = 0,
                   max.overlaps = Inf,
                   min.segment.length = unit(0, "lines"),
                   segment.color = "black") +
  scale_x_continuous("Average Log2FC", breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  scale_y_continuous("-Log10(padj)",
                     limits = c(0,307),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        legend.position = "none") +
  labs(title = "HCC4006 OR2 vs Parental")

hcc4006_p_or2_de_vp #450 x 500 (WXH)



##### HCC4006 Parental vs OR7 DE gene volcano plot #####

#Read in data
hcc4006_p_or7_de <- HCC4006_parental_vs_OR7_DE_genes_volcano_plot

###KRT17 volcano plot

#Identify genes to label
genes_to_label <- c("KRT17")

#Set labels if in DE list
hcc4006_p_or7_de$label <- ifelse(hcc4006_p_or7_de$Gene %in% genes_to_label, hcc4006_p_or7_de$Gene, NA)

#Draw labeled volcano plot for ETV5, HOPX, and SOX4
hcc4006_p_or7_de_vp <- ggplot(hcc4006_p_or7_de, aes(x = avg_log2FC_OR7_v_Parental, y = Minus_log10_padj)) +
  geom_point(aes(color = avg_log2FC_OR7_v_Parental > 0)) +
  scale_color_manual(values = c("FALSE" = "blue",
                                "TRUE" = "red")) +
  geom_point(data = subset(hcc4006_p_or7_de, label != 0),
             shape = 21,
             fill = NA,
             stroke = 1.2,
             color = "black") +
  geom_hline(yintercept = 1.30102999566398, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  geom_vline(xintercept = -1, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  geom_label_repel(aes(label = label),
                   fill = "white",
                   color = "black",
                   box.padding = 1,
                   point.padding = 0,
                   max.overlaps = Inf,
                   min.segment.length = unit(0, "lines"),
                   segment.color = "black") +
  scale_x_continuous("Average Log2FC", breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  scale_y_continuous("-Log10(padj)",
                     limits = c(0,307),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        legend.position = "none") +
  labs(title = "HCC4006 OR7 vs Parental")

hcc4006_p_or7_de_vp #450 x 500 (WXH)







