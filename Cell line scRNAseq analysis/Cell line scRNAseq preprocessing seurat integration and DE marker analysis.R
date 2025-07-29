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



##### HCC827 DMSO Analysis #####

###Data prep & QC filtering

#Define path to single cell data
data_dir <- 'A:/Osi_DTPC_sc/BM_Thoracic-Head-Neck_10XGenSC861/HCC827_DMSO_analysis/outs/filtered_feature_bc_matrix'

#Check that directory contains required files
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

#Read in data using Read10X
expression_matrix <- Read10X(data.dir = data_dir)

#Create & name Seurat Object
hcc827_dmso <- CreateSeuratObject(counts = expression_matrix, project = "HCC827_DMSO")

#Calculate percentage of mitochondrial reads
hcc827_dmso[["percent.mt"]] <- PercentageFeatureSet(hcc827_dmso, pattern = "^MT-")

#Check meta data
head(hcc827_dmso@meta.data, 5)

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (before QC filtering)
unfilt_qc_p <- VlnPlot(hcc827_dmso, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfilt_qc_p

#Visualize number of RNA counts vs %MT reads
nRNA_mt_p <- FeatureScatter(hcc827_dmso, feature1 = "nCount_RNA", feature2 = "percent.mt")
nRNA_mt_p

#Visualize number of RNA counts vs number of unique features
nRNA_nFeat_p <- FeatureScatter(hcc827_dmso, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nRNA_nFeat_p

#Remove cells with < 200 genes detected and/or >20% mitochondrial reads
hcc827_dmso_1 <- subset(hcc827_dmso, subset = nFeature_RNA > 200 & percent.mt < 20) #removed 171 cells, leaving 4,163

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (after QC filtering)
filt_qc_p <- VlnPlot(hcc827_dmso_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filt_qc_p


###SCTransform

#Run scTransform, regressing out percent mitochondrial reads
hcc827_dmso_1 <- SCTransform(hcc827_dmso_1, 
                             vst.flavor = "v2", #specifying v2 for updated vst
                             vars.to.regress = "percent.mt", 
                             verbose = FALSE) #requires matrixStats < version 1.2

###PCA

#Run PCA
hcc827_dmso_1 <- RunPCA(hcc827_dmso_1, verbose = FALSE)

#Identify number of PCs that define >80% of variance
pcs = hcc827_dmso_1@reductions$pca@cell.embeddings
pca_var = hcc827_dmso_1@reductions$pca@stdev ** 2
pca_var_cum = sapply(1:length(pca_var), function(x) sum(pca_var[1:x])/sum(pca_var))
npca = which(pca_var_cum>0.80)[1]
npca 


###Find Neighbors

#Find nearest neighbors
hcc827_dmso_1 <- FindNeighbors(hcc827_dmso_1, dims = 1:npca, verbose = FALSE)


##Find Clusters

#Generate clusters
hcc827_dmso_1 <- FindClusters(hcc827_dmso_1, verbose = FALSE)


###UMAP clustering

#Run UMAP
hcc827_dmso_1 <- RunUMAP(hcc827_dmso_1, dims = 1:npca)

#Draw UMAP
DimPlot(hcc827_dmso_1, label = TRUE, reduction = "umap")


###Doublet removal with DoubletFinder

#PK identification
sweep.res <- paramSweep(hcc827_dmso_1, PCs = 1:npca, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #Plot pk (x) vs BCmvn (y)

#Note: use 6 for doublet identification (value in [] in code line below);
#value changes based on maximum of BCmvn plot above

#Use PK determined from above plot
pK.set <- unique(sweep.stats$pK)[6]

#Doublet proportion estimation
nExp_poi <- round(0.08*nrow(hcc827_dmso_1@meta.data))

#Run DoubletFinder using above metrics
hcc827_dmso_1 <- doubletFinder(hcc827_dmso_1, PCs = 1:npca, pN = 0.25, 
                               pK = as.numeric(as.character(pK.set)), 
                               nExp = nExp_poi, reuse.pANN = FALSE, 
                               sct = TRUE)

#Subset data for only cells likely as singlets
hcc827_dmso_1 <- subset(hcc827_dmso_1,  
                        DF.classifications_0.25_0.05_333 == "Singlet")

#Leaves 3830 cells (removed 503 cells from start, empty drops, dying cells, & doublets)


#Save SeuratRDS (post doublet removal)
SaveSeuratRds(hcc827_dmso_1, "HCC827_DMSO_1_scTransform_v2_singlets.Rds")




##### HCC827 Osi DTPC Analysis #####

#Define path to single cell data
data_dir <- 'A:/Osi_DTPC_sc/BM_Thoracic-Head-Neck_10XGenSC861/HCC827_Osi_DTPC_analysis/outs/filtered_feature_bc_matrix'

#Check that directory contains required files
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

#Read in data using Read10X
expression_matrix <- Read10X(data.dir = data_dir)

#Create & name Seurat Object
hcc827_osi_dtpc <- CreateSeuratObject(counts = expression_matrix, project = "HCC827_Osi_DTPC")

#Calculate percentage of mitochondrial reads
hcc827_osi_dtpc[["percent.mt"]] <- PercentageFeatureSet(hcc827_osi_dtpc, pattern = "^MT-")

#Check meta data
head(hcc827_osi_dtpc@meta.data, 5)

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (before QC filtering)
unfilt_qc_p <- VlnPlot(hcc827_osi_dtpc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfilt_qc_p

#Visualize number of RNA counts vs %MT reads
nRNA_mt_p <- FeatureScatter(hcc827_osi_dtpc, feature1 = "nCount_RNA", feature2 = "percent.mt")
nRNA_mt_p

#Visualize number of RNA counts vs number of unique features
nRNA_nFeat_p <- FeatureScatter(hcc827_osi_dtpc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nRNA_nFeat_p

#Remove cells with < 200 genes detected and/or >20% mitochondrial reads
hcc827_osi_dtpc_1 <- subset(hcc827_osi_dtpc, subset = nFeature_RNA > 200 & percent.mt < 20) #removed 430 cells, leaving 4466

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (after QC filtering)
filt_qc_p <- VlnPlot(hcc827_osi_dtpc_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filt_qc_p


###SCTransform

#Run scTransform, regressing out percent mitochondrial reads
hcc827_osi_dtpc_1 <- SCTransform(hcc827_osi_dtpc_1, 
                                 vst.flavor = "v2", #specifying v2 for updated vst
                                 vars.to.regress = "percent.mt", 
                                 verbose = FALSE) #requires matrixStats < version 1.2

###PCA

#Run PCA
hcc827_osi_dtpc_1 <- RunPCA(hcc827_osi_dtpc_1, verbose = FALSE)

#Identify number of PCs that define >80% of variance
pcs = hcc827_osi_dtpc_1@reductions$pca@cell.embeddings
pca_var = hcc827_osi_dtpc_1@reductions$pca@stdev ** 2
pca_var_cum = sapply(1:length(pca_var), function(x) sum(pca_var[1:x])/sum(pca_var))
npca = which(pca_var_cum>0.80)[1]
npca 


###Find Neighbors

#Find nearest neighbors
hcc827_osi_dtpc_1 <- FindNeighbors(hcc827_osi_dtpc_1, dims = 1:npca, verbose = FALSE)


##Find Clusters

#Generate clusters
hcc827_osi_dtpc_1 <- FindClusters(hcc827_osi_dtpc_1, verbose = FALSE)


###UMAP clustering

#Run UMAP
hcc827_osi_dtpc_1 <- RunUMAP(hcc827_osi_dtpc_1, dims = 1:npca)

#Draw UMAP
DimPlot(hcc827_osi_dtpc_1, label = TRUE, reduction = "umap")


###Doublet removal with DoubletFinder

#PK identification
sweep.res <- paramSweep(hcc827_osi_dtpc_1, PCs = 1:npca, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #Plot pk (x) vs BCmvn (y)

#Note: use 5 for doublet identification (value in [] in code line below);
#value changes based on maximum of BCmvn plot above

#Use PK determined from above plot
pK.set <- unique(sweep.stats$pK)[5]

#Doublet proportion estimation
nExp_poi <- round(0.08*nrow(hcc827_osi_dtpc_1@meta.data))

#Run DoubletFinder using above metrics
hcc827_osi_dtpc_1 <- doubletFinder(hcc827_osi_dtpc_1, PCs = 1:npca, pN = 0.25, 
                                   pK = as.numeric(as.character(pK.set)), 
                                   nExp = nExp_poi, reuse.pANN = FALSE, 
                                   sct = TRUE)

#Subset data for only cells likely as singlets
hcc827_osi_dtpc_1 <- subset(hcc827_osi_dtpc_1,  
                            DF.classifications_0.25_0.04_357 == "Singlet")

#Leaves 4109 cells (removed 357 cells from start, empty drops, dying cells, & doublets)


#Save SeuratRDS (post doublet removal)
SaveSeuratRds(hcc827_osi_dtpc_1, "HCC827_Osi_DTPC_1_scTransform_v2_singlets.Rds")



##### HCC827 DMSO and HCC827 Osi DTPC integrated analysis  #####

#Load data
hcc827_dmso_1 <- LoadSeuratRds("C:/Users/bbmorris1/Desktop/HCC827_DMSO_1_scTransform_v2_singlets.Rds")
hcc827_osi_dtpc_1 <- LoadSeuratRds("C:/Users/bbmorris1/Desktop/HCC827_Osi_DTPC_1_scTransform_v2_singlets.Rds")

#Merge two Seurat objects
hcc827_combined <- merge(x = hcc827_dmso_1, y = hcc827_osi_dtpc_1, add.cell.ids = c("DMSO", "Osi_DTPC"), project = "HCC827", merge.data = TRUE)
hcc827_combined

#Set SCT transformed/normalized values as VariableFeatures, required for PCA
VariableFeatures(hcc827_combined[["SCT"]]) <- rownames(hcc827_combined[["SCT"]]@scale.data)

#Note: setting SCT data as variable features is required for SelectIntegrationFeatures & PrepSCTIntegration

#Split object
hcc827_cond_list <- SplitObject(hcc827_combined, split.by = "orig.ident")

dmso <- hcc827_cond_list[["HCC827_DMSO"]]
osi_dtpc <- hcc827_cond_list[["HCC827_Osi_DTPC"]]

#Select integration features
features <- SelectIntegrationFeatures(object.list = hcc827_cond_list, nfeatures = 3000)
hcc827_cond_list <- PrepSCTIntegration(object.list = hcc827_cond_list, anchor.features = features)

#Find integration anchors
hcc827_anchors <- FindIntegrationAnchors(object.list = hcc827_cond_list, normalization.method = "SCT",
                                         anchor.features = features)

#Integrate data
hcc827_combined_sct <- IntegrateData(anchorset = hcc827_anchors, normalization.method = "SCT")


#Run PCA
hcc827_combined_sct <- RunPCA(hcc827_combined_sct, verbose = FALSE)

#Run UMAP
hcc827_combined_sct <- RunUMAP(hcc827_combined_sct, reduction = "pca", dims = 1:30, verbose = FALSE)

#Find neighbors (SNN)
hcc827_combined_sct <- FindNeighbors(hcc827_combined_sct, reduction = "pca", dims = 1:30)

#Find clusters
hcc827_combined_sct <- FindClusters(hcc827_combined_sct, resolution = 1)


#Plot UMAP, false coloring by sample condition
DimPlot(hcc827_combined_sct, reduction = "umap", group.by = "orig.ident")

#Plot UMAP, false coloring by seurat cluster
seurat_cluster_umap <- DimPlot(hcc827_combined_sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                               repel = TRUE)

seurat_cluster_umap



#Save SeurateRDS
SaveSeuratRds(hcc827_combined_sct, "HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_annotated.Rds")



###DE marker gene analysis (DMSO vs Osi DTPC)
Idents(hcc827_combined_sct) <- "orig.ident"

hcc827_combined_sct %<>% PrepSCTFindMarkers(assay = "SCT")

hcc827_dmso_vs_osi_response <- FindMarkers(hcc827_combined_sct, assay = "SCT",
                                           recorrect_umi = FALSE,
                                           ident.1 = "HCC827_Osi_DTPC", 
                                           ident.2 = "HCC827_DMSO",
                                           verbose = FALSE)

#Check beginning of marker gene table
head(hcc827_dmso_vs_osi_response)

#Write marker gene table to file
write.table(hcc827_dmso_vs_osi_response, "C:/Users/bbmorris1/Desktop/HCC827 Osi DTPC vs HCC827 DMSO marker genes.txt", sep = "\t", quote = FALSE)



###DE cluster markers gene analysis
Idents(hcc827_combined_sct) <- "seurat_clusters"

hcc827_combined_sct %<>% PrepSCTFindMarkers(assay = "SCT")

#FindAllMarkers (requires normalizing library sizes first)
cluster_markers <- FindAllMarkers(hcc827_combined_sct, assay = "SCT", recorrect_umi = FALSE)

#Check beginning of marker gene table
head(cluster_markers)

#Write marker gene table to file
write.table(cluster_markers, "C:/Users/bbmorris1/Desktop/HCC827 DMSO Osi DTPC seurat cluster marker genes.txt", sep = "\t", quote = FALSE)


#Draw cluster heatamp
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

hcc827_combined_sc_mkr_h <- DoHeatmap(hcc827_combined_sct, features = top10$gene, raster = FALSE) + 
  NoLegend() +
  theme(axis.text.y = element_text(color = "black", size = 4))


hcc827_combined_sc_mkr_h #export as PDF (SVGs raster really slow)



##### HCC4006 DMSO Analysis #####

###Data prep & QC filtering

#Define path to single cell data
data_dir <- 'A:/Osi_DTPC_sc/BM_Thoracic-Head-Neck_10XGenSC861/HCC4006_DMSO_analysis/outs/filtered_feature_bc_matrix'

#Check that directory contains required files
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

#Read in data using Read10X
expression_matrix <- Read10X(data.dir = data_dir)

#Create & name Seurat Object
hcc4006_dmso <- CreateSeuratObject(counts = expression_matrix, project = "HCC4006_DMSO")

#Calculate percentage of mitochondrial reads
hcc4006_dmso[["percent.mt"]] <- PercentageFeatureSet(hcc4006_dmso, pattern = "^MT-")

#Check meta data
head(hcc4006_dmso@meta.data, 5)

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (before QC filtering)
unfilt_qc_p <- VlnPlot(hcc4006_dmso, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfilt_qc_p

#Visualize number of RNA counts vs %MT reads
nRNA_mt_p <- FeatureScatter(hcc4006_dmso, feature1 = "nCount_RNA", feature2 = "percent.mt")
nRNA_mt_p

#Visualize number of RNA counts vs number of unique features
nRNA_nFeat_p <- FeatureScatter(hcc4006_dmso, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nRNA_nFeat_p

#Remove cells with < 200 genes detected and/or >20% mitochondrial reads
hcc4006_dmso_1 <- subset(hcc4006_dmso, subset = nFeature_RNA > 200 & percent.mt < 20) #removed 87 cells, leaving 3155

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (after QC filtering)
filt_qc_p <- VlnPlot(hcc4006_dmso_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filt_qc_p


###SCTransform

#Run scTransform, regressing out percent mitochondrial reads
hcc4006_dmso_1 <- SCTransform(hcc4006_dmso_1, 
                              vst.flavor = "v2", #specifying v2 for updated vst
                              vars.to.regress = "percent.mt", 
                              verbose = FALSE) #requires matrixStats < version 1.2

###PCA

#Run PCA
hcc4006_dmso_1 <- RunPCA(hcc4006_dmso_1, verbose = FALSE)

#Identify number of PCs that define >80% of variance
pcs = hcc4006_dmso_1@reductions$pca@cell.embeddings
pca_var = hcc4006_dmso_1@reductions$pca@stdev ** 2
pca_var_cum = sapply(1:length(pca_var), function(x) sum(pca_var[1:x])/sum(pca_var))
npca = which(pca_var_cum>0.80)[1]
npca 


###Find Neighbors

#Find nearest neighbors
hcc4006_dmso_1 <- FindNeighbors(hcc4006_dmso_1, dims = 1:npca, verbose = FALSE)


##Find Clusters

#Generate clusters
hcc4006_dmso_1 <- FindClusters(hcc4006_dmso_1, verbose = FALSE)


###UMAP clustering

#Run UMAP
hcc4006_dmso_1 <- RunUMAP(hcc4006_dmso_1, dims = 1:npca)

#Draw UMAP
DimPlot(hcc4006_dmso_1, label = TRUE, reduction = "umap")


###Doublet removal with DoubletFinder

#PK identification
sweep.res <- paramSweep(hcc4006_dmso_1, PCs = 1:npca, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #Plot pk (x) vs BCmvn (y)

#Note: use 16 for doublet identification (value in [] in code line below);
#value changes based on maximum of BCmvn plot above

#Use PK determined from above plot
pK.set <- unique(sweep.stats$pK)[16]

#Doublet proportion estimation
nExp_poi <- round(0.08*nrow(hcc4006_dmso_1@meta.data))

#Run DoubletFinder using above metrics
hcc4006_dmso_1 <- doubletFinder(hcc4006_dmso_1, PCs = 1:npca, pN = 0.25, 
                                pK = as.numeric(as.character(pK.set)), 
                                nExp = nExp_poi, reuse.pANN = FALSE, 
                                sct = TRUE)

#Subset data for only cells likely as singlets
hcc4006_dmso_1 <- subset(hcc4006_dmso_1,  
                         DF.classifications_0.25_0.15_252 == "Singlet")

#Leaves 2903 cells (removed 339 cells from start, empty drops, dying cells, & doublets)


#Save SeuratRDS (post doublet removal)
SaveSeuratRds(hcc4006_dmso_1, "HCC4006_DMSO_1_scTransform_v2_singlets.Rds")


##### HCC4006 Osi DTPC Analysis #####

#Define path to single cell data
data_dir <- 'A:/Osi_DTPC_sc/BM_Thoracic-Head-Neck_10XGenSC861/HCC4006_Osi_DTPC_analysis/outs/filtered_feature_bc_matrix'

#Check that directory contains required files
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

#Read in data using Read10X
expression_matrix <- Read10X(data.dir = data_dir)

#Create & name Seurat Object
hcc4006_osi_dtpc <- CreateSeuratObject(counts = expression_matrix, project = "HCC4006_Osi_DTPC")

#Calculate percentage of mitochondrial reads
hcc4006_osi_dtpc[["percent.mt"]] <- PercentageFeatureSet(hcc4006_osi_dtpc, pattern = "^MT-")

#Check meta data
head(hcc4006_osi_dtpc@meta.data, 5)

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (before QC filtering)
unfilt_qc_p <- VlnPlot(hcc4006_osi_dtpc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfilt_qc_p

#Visualize number of RNA counts vs %MT reads
nRNA_mt_p <- FeatureScatter(hcc4006_osi_dtpc, feature1 = "nCount_RNA", feature2 = "percent.mt")
nRNA_mt_p

#Visualize number of RNA counts vs number of unique features
nRNA_nFeat_p <- FeatureScatter(hcc4006_osi_dtpc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nRNA_nFeat_p

#Remove cells with < 200 genes detected and/or >20% mitochondrial reads
hcc4006_osi_dtpc_1 <- subset(hcc4006_osi_dtpc, subset = nFeature_RNA > 200 & percent.mt < 20) #removed  cells, leaving 

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (after QC filtering)
filt_qc_p <- VlnPlot(hcc4006_osi_dtpc_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filt_qc_p


###SCTransform

#Run scTransform, regressing out percent mitochondrial reads
hcc4006_osi_dtpc_1 <- SCTransform(hcc4006_osi_dtpc_1, 
                                  vst.flavor = "v2", #specifying v2 for updated vst
                                  vars.to.regress = "percent.mt", 
                                  verbose = FALSE) #requires matrixStats < version 1.2

###PCA

#Run PCA
hcc4006_osi_dtpc_1 <- RunPCA(hcc4006_osi_dtpc_1, verbose = FALSE)

#Identify number of PCs that define >80% of variance
pcs = hcc4006_osi_dtpc_1@reductions$pca@cell.embeddings
pca_var = hcc4006_osi_dtpc_1@reductions$pca@stdev ** 2
pca_var_cum = sapply(1:length(pca_var), function(x) sum(pca_var[1:x])/sum(pca_var))
npca = which(pca_var_cum>0.80)[1]
npca 


###Find Neighbors

#Find nearest neighbors
hcc4006_osi_dtpc_1 <- FindNeighbors(hcc4006_osi_dtpc_1, dims = 1:npca, verbose = FALSE)


##Find Clusters

#Generate clusters
hcc4006_osi_dtpc_1 <- FindClusters(hcc4006_osi_dtpc_1, verbose = FALSE)


###UMAP clustering

#Run UMAP
hcc4006_osi_dtpc_1 <- RunUMAP(hcc4006_osi_dtpc_1, dims = 1:npca)

#Draw UMAP
DimPlot(hcc4006_osi_dtpc_1, label = TRUE, reduction = "umap")


###Doublet removal with DoubletFinder

#PK identification
sweep.res <- paramSweep(hcc4006_osi_dtpc_1, PCs = 1:npca, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #Plot pk (x) vs BCmvn (y)

#Note: use 5 for doublet identification (value in [] in code line below);
#value changes based on maximum of BCmvn plot above

#Use PK determined from above plot
pK.set <- unique(sweep.stats$pK)[8]

#Doublet proportion estimation
nExp_poi <- round(0.08*nrow(hcc4006_osi_dtpc_1@meta.data))

#Run DoubletFinder using above metrics
hcc4006_osi_dtpc_1 <- doubletFinder(hcc4006_osi_dtpc_1, PCs = 1:npca, pN = 0.25, 
                                    pK = as.numeric(as.character(pK.set)), 
                                    nExp = nExp_poi, reuse.pANN = FALSE, 
                                    sct = TRUE)

#Subset data for only cells likely as singlets
hcc4006_osi_dtpc_1 <- subset(hcc4006_osi_dtpc_1,  
                             DF.classifications_0.25_0.07_458 == "Singlet")

#Leaves 5267 cells (removed 573 cells from start, empty drops, dying cells, & doublets)


#Save SeuratRDS (post doublet removal)
SaveSeuratRds(hcc4006_osi_dtpc_1, "HCC4006_Osi_DTPC_1_scTransform_v2_singlets.Rds")



##### HCC4006 DMSO and HCC4006 Osi DTPC integrated analysis  #####

#Load data
hcc4006_dmso_1 <- LoadSeuratRds("C:/Users/bbmorris1/Desktop/HCC4006_DMSO_1_scTransform_v2_singlets.Rds")
hcc4006_osi_dtpc_1 <- LoadSeuratRds("C:/Users/bbmorris1/Desktop/HCC4006_Osi_DTPC_1_scTransform_v2_singlets.Rds")

#Merge two Seurat objects
hcc4006_combined <- merge(x = hcc4006_dmso_1, y = hcc4006_osi_dtpc_1, add.cell.ids = c("DMSO", "Osi_DTPC"), project = "HCC4006", merge.data = TRUE)
hcc4006_combined

#Set SCT transformed/normalized values as VariableFeatures, required for PCA
VariableFeatures(hcc4006_combined[["SCT"]]) <- rownames(hcc4006_combined[["SCT"]]@scale.data)

#Note: setting SCT data as variable features is required for SelectIntegrationFeatures & PrepSCTIntegration

#Split object
hcc4006_cond_list <- SplitObject(hcc4006_combined, split.by = "orig.ident")

dmso <- hcc4006_cond_list[["HCC4006_DMSO"]]
osi_dtpc <- hcc4006_cond_list[["HCC4006_Osi_DTPC"]]

#Select integration features
features <- SelectIntegrationFeatures(object.list = hcc4006_cond_list, nfeatures = 3000)
hcc4006_cond_list <- PrepSCTIntegration(object.list = hcc4006_cond_list, anchor.features = features)

#Find integration anchors
hcc4006_anchors <- FindIntegrationAnchors(object.list = hcc4006_cond_list, normalization.method = "SCT",
                                          anchor.features = features)

#Integrate data
hcc4006_combined_sct <- IntegrateData(anchorset = hcc4006_anchors, normalization.method = "SCT")


#Run PCA
hcc4006_combined_sct <- RunPCA(hcc4006_combined_sct, verbose = FALSE)

#Run UMAP
hcc4006_combined_sct <- RunUMAP(hcc4006_combined_sct, reduction = "pca", dims = 1:30, verbose = FALSE)

#Find neighbors (SNN)
hcc4006_combined_sct <- FindNeighbors(hcc4006_combined_sct, reduction = "pca", dims = 1:30)

#Find clusters
hcc4006_combined_sct <- FindClusters(hcc4006_combined_sct, resolution = 1)


#Plot UMAP, false coloring by sample condition
DimPlot(hcc4006_combined_sct, reduction = "umap", group.by = "orig.ident")

#Plot UMAP, false coloring by seurat cluster
seurat_cluster_umap <- DimPlot(hcc4006_combined_sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                               repel = TRUE)

seurat_cluster_umap



#Save SeurateRDS
SaveSeuratRds(hcc4006_combined_sct, "HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_annotated.Rds")



###DE marker gene analysis (DMSO vs Osi DTPC)
Idents(hcc4006_combined_sct) <- "orig.ident"

hcc4006_combined_sct %<>% PrepSCTFindMarkers(assay = "SCT")

hcc4006_dmso_vs_osi_response <- FindMarkers(hcc4006_combined_sct, assay = "SCT",
                                            recorrect_umi = FALSE,
                                            ident.1 = "HCC4006_Osi_DTPC", 
                                            ident.2 = "HCC4006_DMSO",
                                            verbose = FALSE)

#Check beginning of marker gene table
head(hcc4006_dmso_vs_osi_response)

#Write marker gene table to file
write.table(hcc4006_dmso_vs_osi_response, "C:/Users/bbmorris1/Desktop/HCC4006 Osi DTPC vs HCC4006 DMSO marker genes.txt", sep = "\t", quote = FALSE)



###DE cluster markers gene analysis
Idents(hcc4006_combined_sct) <- "seurat_clusters"

hcc4006_combined_sct %<>% PrepSCTFindMarkers(assay = "SCT")

#FindAllMarkers (requires normalizing library sizes first)
cluster_markers <- FindAllMarkers(hcc4006_combined_sct, assay = "SCT", recorrect_umi = FALSE)

#Check beginning of marker gene table
head(cluster_markers)

#Write marker gene table to file
write.table(cluster_markers, "C:/Users/bbmorris1/Desktop/HCC4006 DMSO Osi DTPC seurat cluster marker genes.txt", sep = "\t", quote = FALSE)


#Draw cluster heatamp
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

hcc4006_combined_sc_mkr_h <- DoHeatmap(hcc4006_combined_sct, features = top10$gene, raster = FALSE) + 
  NoLegend() +
  theme(axis.text.y = element_text(color = "black", size = 4))


hcc4006_combined_sc_mkr_h #export as PDF (SVGs raster really slow)



##### H1975 DMSO Analysis #####

###Data prep & QC filtering

#Define path to single cell data
data_dir <- 'A:/Osi_DTPC_sc/BM_Thoracic-Head-Neck_10XGenSC861/H1975_DMSO_analysis/outs/filtered_feature_bc_matrix'

#Check that directory contains required files
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

#Read in data using Read10X
expression_matrix <- Read10X(data.dir = data_dir)

#Create & name Seurat Object
h1975_dmso <- CreateSeuratObject(counts = expression_matrix, project = "H1975_DMSO")

#Calculate percentage of mitochondrial reads
h1975_dmso[["percent.mt"]] <- PercentageFeatureSet(h1975_dmso, pattern = "^MT-")

#Check meta data
head(h1975_dmso@meta.data, 5)

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (before QC filtering)
unfilt_qc_p <- VlnPlot(h1975_dmso, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfilt_qc_p

#Visualize number of RNA counts vs %MT reads
nRNA_mt_p <- FeatureScatter(h1975_dmso, feature1 = "nCount_RNA", feature2 = "percent.mt")
nRNA_mt_p

#Visualize number of RNA counts vs number of unique features
nRNA_nFeat_p <- FeatureScatter(h1975_dmso, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nRNA_nFeat_p

#Remove cells with < 200 genes detected and/or >20% mitochondrial reads
h1975_dmso_1 <- subset(h1975_dmso, subset = nFeature_RNA > 200 & percent.mt < 20) #removed 47 cells, leaving 5937 cells

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (after QC filtering)
filt_qc_p <- VlnPlot(h1975_dmso_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filt_qc_p


###SCTransform

#Run scTransform, regressing out percent mitochondrial reads
h1975_dmso_1 <- SCTransform(h1975_dmso_1, 
                            vst.flavor = "v2", #specifying v2 for updated vst
                            vars.to.regress = "percent.mt", 
                            verbose = FALSE) #requires matrixStats < version 1.2

###PCA

#Run PCA
h1975_dmso_1 <- RunPCA(h1975_dmso_1, verbose = FALSE)

#Identify number of PCs that define >80% of variance
pcs = h1975_dmso_1@reductions$pca@cell.embeddings
pca_var = h1975_dmso_1@reductions$pca@stdev ** 2
pca_var_cum = sapply(1:length(pca_var), function(x) sum(pca_var[1:x])/sum(pca_var))
npca = which(pca_var_cum>0.80)[1]
npca 


###Find Neighbors

#Find nearest neighbors
h1975_dmso_1 <- FindNeighbors(h1975_dmso_1, dims = 1:npca, verbose = FALSE)


##Find Clusters

#Generate clusters
h1975_dmso_1 <- FindClusters(h1975_dmso_1, verbose = FALSE)


###UMAP clustering

#Run UMAP
h1975_dmso_1 <- RunUMAP(h1975_dmso_1, dims = 1:npca)

#Draw UMAP
DimPlot(h1975_dmso_1, label = TRUE, reduction = "umap")


###Doublet removal with DoubletFinder

#PK identification
sweep.res <- paramSweep(h1975_dmso_1, PCs = 1:npca, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #Plot pk (x) vs BCmvn (y)

#Note: use 6 for doublet identification (value in [] in code line below);
#value changes based on maximum of BCmvn plot above

#Use PK determined from above plot
pK.set <- unique(sweep.stats$pK)[14]

#Doublet proportion estimation
nExp_poi <- round(0.08*nrow(h1975_dmso_1@meta.data))

#Run DoubletFinder using above metrics
h1975_dmso_1 <- doubletFinder(h1975_dmso_1, PCs = 1:npca, pN = 0.25, 
                              pK = as.numeric(as.character(pK.set)), 
                              nExp = nExp_poi, reuse.pANN = FALSE, 
                              sct = TRUE)

#Subset data for only cells likely as singlets
h1975_dmso_1 <- subset(h1975_dmso_1,  
                       DF.classifications_0.25_0.13_475 == "Singlet")

#Leaves 5462 cells (removed  cells from start, empty drops, dying cells, & doublets)


#Save SeuratRDS (post doublet removal)
SaveSeuratRds(h1975_dmso_1, "H1975_DMSO_1_scTransform_v2_singlets.Rds")



##### H1975 Osi DTPC Analysis #####

#Define path to single cell data
data_dir <- 'A:/Osi_DTPC_sc/BM_Thoracic-Head-Neck_10XGenSC861/H1975_Osi_DTPC_analysis/outs/filtered_feature_bc_matrix'

#Check that directory contains required files
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

#Read in data using Read10X
expression_matrix <- Read10X(data.dir = data_dir)

#Create & name Seurat Object
h1975_osi_dtpc <- CreateSeuratObject(counts = expression_matrix, project = "H1975_Osi_DTPC")

#Calculate percentage of mitochondrial reads
h1975_osi_dtpc[["percent.mt"]] <- PercentageFeatureSet(h1975_osi_dtpc, pattern = "^MT-")

#Check meta data
head(h1975_osi_dtpc@meta.data, 5)

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (before QC filtering)
unfilt_qc_p <- VlnPlot(h1975_osi_dtpc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfilt_qc_p

#Visualize number of RNA counts vs %MT reads
nRNA_mt_p <- FeatureScatter(h1975_osi_dtpc, feature1 = "nCount_RNA", feature2 = "percent.mt")
nRNA_mt_p

#Visualize number of RNA counts vs number of unique features
nRNA_nFeat_p <- FeatureScatter(h1975_osi_dtpc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nRNA_nFeat_p

#Remove cells with < 200 genes detected and/or >20% mitochondrial reads
h1975_osi_dtpc_1 <- subset(h1975_osi_dtpc, subset = nFeature_RNA > 200 & percent.mt < 20) #removed 430 cells, leaving 4466

#Draw QC violin plots for #features, #RNAs, and %MT reads per cell in this sample (after QC filtering)
filt_qc_p <- VlnPlot(h1975_osi_dtpc_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filt_qc_p


###SCTransform

#Run scTransform, regressing out percent mitochondrial reads
h1975_osi_dtpc_1 <- SCTransform(h1975_osi_dtpc_1, 
                                vst.flavor = "v2", #specifying v2 for updated vst
                                vars.to.regress = "percent.mt", 
                                verbose = FALSE) #requires matrixStats < version 1.2

###PCA

#Run PCA
h1975_osi_dtpc_1 <- RunPCA(h1975_osi_dtpc_1, verbose = FALSE)

#Identify number of PCs that define >80% of variance
pcs = h1975_osi_dtpc_1@reductions$pca@cell.embeddings
pca_var = h1975_osi_dtpc_1@reductions$pca@stdev ** 2
pca_var_cum = sapply(1:length(pca_var), function(x) sum(pca_var[1:x])/sum(pca_var))
npca = which(pca_var_cum>0.80)[1]
npca 


###Find Neighbors

#Find nearest neighbors
h1975_osi_dtpc_1 <- FindNeighbors(h1975_osi_dtpc_1, dims = 1:npca, verbose = FALSE)


##Find Clusters

#Generate clusters
h1975_osi_dtpc_1 <- FindClusters(h1975_osi_dtpc_1, verbose = FALSE)


###UMAP clustering

#Run UMAP
h1975_osi_dtpc_1 <- RunUMAP(h1975_osi_dtpc_1, dims = 1:npca)

#Draw UMAP
DimPlot(h1975_osi_dtpc_1, label = TRUE, reduction = "umap")


###Doublet removal with DoubletFinder

#PK identification
sweep.res <- paramSweep(h1975_osi_dtpc_1, PCs = 1:npca, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #Plot pk (x) vs BCmvn (y)

#Note: use 25 for doublet identification (value in [] in code line below);
#value changes based on maximum of BCmvn plot above

#Use PK determined from above plot
pK.set <- unique(sweep.stats$pK)[25]

#Doublet proportion estimation
nExp_poi <- round(0.08*nrow(h1975_osi_dtpc_1@meta.data))

#Run DoubletFinder using above metrics
h1975_osi_dtpc_1 <- doubletFinder(h1975_osi_dtpc_1, PCs = 1:npca, pN = 0.25, 
                                  pK = as.numeric(as.character(pK.set)), 
                                  nExp = nExp_poi, reuse.pANN = FALSE, 
                                  sct = TRUE)

#Subset data for only cells likely as singlets
h1975_osi_dtpc_1 <- subset(h1975_osi_dtpc_1,  
                           DF.classifications_0.25_0.24_502 == "Singlet")

#Leaves 5767 cells (removed  cells from start, empty drops, dying cells, & doublets)


#Save SeuratRDS (post doublet removal)
SaveSeuratRds(h1975_osi_dtpc_1, "H1975_Osi_DTPC_1_scTransform_v2_singlets.Rds")



##### H1975 DMSO and H1975 Osi DTPC integrated analysis  #####

#Load data
h1975_dmso_1 <- LoadSeuratRds("C:/Users/bbmorris1/Desktop/H1975_DMSO_1_scTransform_v2_singlets.Rds")
h1975_osi_dtpc_1 <- LoadSeuratRds("C:/Users/bbmorris1/Desktop/H1975_Osi_DTPC_1_scTransform_v2_singlets.Rds")

#Merge two Seurat objects
h1975_combined <- merge(x = h1975_dmso_1, y = h1975_osi_dtpc_1, add.cell.ids = c("DMSO", "Osi_DTPC"), project = "H1975", merge.data = TRUE)
h1975_combined

#Set SCT transformed/normalized values as VariableFeatures, required for PCA
VariableFeatures(h1975_combined[["SCT"]]) <- rownames(h1975_combined[["SCT"]]@scale.data)

#Note: setting SCT data as variable features is required for SelectIntegrationFeatures & PrepSCTIntegration

#Split object
h1975_cond_list <- SplitObject(h1975_combined, split.by = "orig.ident")

dmso <- h1975_cond_list[["H1975_DMSO"]]
osi_dtpc <- h1975_cond_list[["H1975_Osi_DTPC"]]

#Select integration features
features <- SelectIntegrationFeatures(object.list = h1975_cond_list, nfeatures = 3000)
h1975_cond_list <- PrepSCTIntegration(object.list = h1975_cond_list, anchor.features = features)

#Find integration anchors
h1975_anchors <- FindIntegrationAnchors(object.list = h1975_cond_list, normalization.method = "SCT",
                                        anchor.features = features)

#Integrate data
h1975_combined_sct <- IntegrateData(anchorset = h1975_anchors, normalization.method = "SCT")


#Run PCA
h1975_combined_sct <- RunPCA(h1975_combined_sct, verbose = FALSE)

#Run UMAP
h1975_combined_sct <- RunUMAP(h1975_combined_sct, reduction = "pca", dims = 1:30, verbose = FALSE)

#Find neighbors (SNN)
h1975_combined_sct <- FindNeighbors(h1975_combined_sct, reduction = "pca", dims = 1:30)

#Find clusters
h1975_combined_sct <- FindClusters(h1975_combined_sct, resolution = 1)


#Plot UMAP, false coloring by sample condition
DimPlot(h1975_combined_sct, reduction = "umap", group.by = "orig.ident")

#Plot UMAP, false coloring by seurat cluster
seurat_cluster_umap <- DimPlot(h1975_combined_sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                               repel = TRUE)

seurat_cluster_umap



#Save SeurateRDS
SaveSeuratRds(h1975_combined_sct, "H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_annotated.Rds")



###DE marker gene analysis (DMSO vs Osi DTPC)
Idents(h1975_combined_sct) <- "orig.ident"

h1975_combined_sct %<>% PrepSCTFindMarkers(assay = "SCT")

h1975_dmso_vs_osi_response <- FindMarkers(h1975_combined_sct, assay = "SCT",
                                          recorrect_umi = FALSE,
                                          ident.1 = "H1975_Osi_DTPC", 
                                          ident.2 = "H1975_DMSO",
                                          verbose = FALSE)

#Check beginning of marker gene table
head(h1975_dmso_vs_osi_response)

#Write marker gene table to file
write.table(h1975_dmso_vs_osi_response, "C:/Users/bbmorris1/Desktop/H1975 Osi DTPC vs H1975 DMSO marker genes.txt", sep = "\t", quote = FALSE)



###DE cluster markers gene analysis
Idents(h1975_combined_sct) <- "seurat_clusters"

h1975_combined_sct %<>% PrepSCTFindMarkers(assay = "SCT")

#FindAllMarkers (requires normalizing library sizes first)
cluster_markers <- FindAllMarkers(h1975_combined_sct, assay = "SCT", recorrect_umi = FALSE)

#Check beginning of marker gene table
head(cluster_markers)

#Write marker gene table to file
write.table(cluster_markers, "C:/Users/bbmorris1/Desktop/H1975 DMSO Osi DTPC seurat cluster marker genes.txt", sep = "\t", quote = FALSE)


#Draw cluster heatamp
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

h1975_combined_sc_mkr_h <- DoHeatmap(h1975_combined_sct, features = top10$gene, raster = FALSE) + 
  NoLegend() +
  theme(axis.text.y = element_text(color = "black", size = 4))


h1975_combined_sc_mkr_h #export as PDF (SVGs raster really slow)





