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
library(scales)
library(patchwork)
library(ggrepel)


##### Expanded clinical cohort object construction #####

library(Matrix)
library(SingleCellExperiment)
library(Seurat)
library(data.table)

###Raw counts matrix annotation

#Read in raw counts
raw_counts <- readMM("counts.mtx")
raw_counts

#Transpose counts to correct format (genes x cells)
raw_counts_t <- t(raw_counts)
raw_counts_t

#Read in cell metadata
cell_meta <- fread("counts_cellMeta.csv")

#Read in gene metadata
gene_meta <- fread("counts_geneMeta.csv")

#Set rownames for raw counts t
rownames(raw_counts_t) <- gene_meta$GeneName

#Set colanmes for raw counts t
colnames(raw_counts_t) <- cell_meta$Barcode


###Log2 normalized counts matrix annotation

#Read in log normalized counts
log2norm_counts <- readMM("log2norm.mtx")
log2norm_counts

#Transpose counts to correct format (genes x cells)
log2norm_counts_t <- t(log2norm_counts)
log2norm_counts_t

#Set rownames for log2norm counts t
rownames(log2norm_counts_t) <- gene_meta$GeneName

#Set colanmes for log2norm counts t
colnames(log2norm_counts_t) <- cell_meta$Barcode


#Check column names (cells) exactly match the order of the cell metadata
unique(colnames(log2norm_counts_t) == cell_meta$Barcode) #Returns: TRUE, meaning everything matches
unique(colnames(raw_counts_t) == cell_meta$Barcode) #Returns: TRUE, meaning everything matches
unique(colnames(log2norm_counts_t) == colnames(raw_counts_t)) #Returns: TRUE, meaning everything matches

#Check rownames (genes) exactly match the order of the gene metadata file
unique(rownames(log2norm_counts_t) == gene_meta$GeneName) #Returns: TRUE, meaning everything matches
unique(rownames(raw_counts_t) == gene_meta$GeneName) #Returns: TRUE, meaning everything matches
unique(rownames(log2norm_counts_t) == rownames(raw_counts_t)) #Returns: TRUE, meaning everything matches


#Create single cell experiment object using counts & log2norm data, including with cell level metadata
egfr_sce <- SingleCellExperiment(assays = list(counts = raw_counts_t, log2norm = log2norm_counts_t),
                                 colData = cell_meta)

#Check rownames (gene names) of single cell experiment object
rownames(egfr_sce)

#Check column names (cell barcodes) of single cell experiment object
colnames(egfr_sce)


#Extract count data matrix
egfr_counts_mtx <- assay(egfr_sce, "counts")

#Extract log2norm count data matrix
egfr_log2norm_mtx <- assay(egfr_sce, "log2norm")


#Create Seurat object
egfr_seurat <- CreateSeuratObject(counts = egfr_counts_mtx,
                                  meta.data = as.data.frame(colData(egfr_sce)))


#Add log2norm count data back 
egfr_seurat[["RNA"]] <- SetAssayData(egfr_seurat[["RNA"]],
                                     layer = "data",
                                     new.data = egfr_log2norm_mtx)



#Save combined Seurat object
SaveSeuratRds(egfr_seurat, "EGFR_combined_batches.Rds")



##### Expanded cohort treatment naive sample analysis #####

#Read in data
egfr_seurat <- LoadSeuratRds("EGFR_combined_batches.Rds")

#Filter for treatment naive tumor samples
specimen_ids_to_subset <- c("Biopsy1", "Lung-tumor-10", "JH064",
                            "JH067", "JH095", "JH128", "JH139", "JH338",
                            "JH340", "JH348", "JH384", "JH395", "JH400")

#Subset for treatment naive samples
egfr_epi_tn <- subset(egfr_seurat, subset = batch %in% specimen_ids_to_subset)

#Filter for epithelial cells
egfr_epi_tn <- subset(egfr_epi_tn, subset = General_Annot == "EPITHELIAL")

#Filter for InferCnv annotated tumor cells
egfr_epi_t_tn <- subset(egfr_epi_tn, subset = InfercnvR == "tumor")


#Summarize the number of epithelial tumor cells by Specimen_ID
egfr_epi_stats <- egfr_epi_t_tn@meta.data %>%
  group_by(batch) %>%
  summarise(Epithelial_tumor_cells = n())

write.table(egfr_epi_stats, "MDACC EGFR mutant scRNAseq expanded cohort treatment naive sample tumor cell count by batch ID.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Filter for treatment naive tumor samples with >20 tumor cells
specimen_ids_to_subset <- c("Biopsy1", "Lung-tumor-10", "JH064",
                            "JH067", "JH095", "JH128", "JH139", 
                            "JH348", "JH384", "JH400") #Excluding JH340 as it is LUSC at diagnosis

#Subset for epithelial TN > 20 cells only
egfr_epi_tn <- subset(egfr_seurat, subset = batch %in% specimen_ids_to_subset)

#Filter for epithelial cells
egfr_epi_tn <- subset(egfr_epi_tn, subset = General_Annot == "EPITHELIAL")

#Filter for InferCnv annotated tumor cells
egfr_epi_t_tn <- subset(egfr_epi_tn, subset = InfercnvR == "tumor")

#Change active identity to batch (Specimen ID)
Idents(egfr_epi_t_tn) <- "batch"

#Draw KRT17 expression violin plot across treatment naive samples
krt17_tn_vp <- VlnPlot(egfr_epi_t_tn, features = "KRT17", group.by = "batch", assay = "RNA", layer = "data") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("Biopsy-1" = "#F8766D",
                               "JH064" = "#E9842C",
                               "JH067" = '#D69100',
                               "JH095" = '#BC9D00',
                               "JH128" = '#9CA700',
                               "JH139" = '#6FB000',
                               "Lung-Tumor-10" = '#00BD61',
                               "JH348" = "#00BA38",
                               "JH384" = "#00BE67",
                               "JH400" = "#49B500")) +
  labs(title = "Treatment Naive KRT17 expression")

krt17_tn_vp #Plot: 700 x 400



###Annotate KRT17+ vs KRT17- cells by batch ID

#Label KRT17+ cells
egfr_epi_t_tn$KRT17_exp <- ifelse(FetchData(egfr_epi_t_tn, vars = "KRT17") > 0, 
                                  "KRT17+", "KRT17-")

head(egfr_epi_t_tn)

#Summarize the number of KRT17+ and KRT17- cells by batch ID
egfr_epi_t_tn_krt17_stats <- table(egfr_epi_t_tn@meta.data$batch, egfr_epi_t_tn@meta.data$KRT17_exp)
egfr_epi_t_tn_krt17_stats

#Write data to table
write.table(egfr_epi_t_tn_krt17_stats, "MDACC EGFR mutant scRNAseq expanded cohort treatment naive KRT17 expression counts.txt", sep = "\t", row.names = TRUE, quote = FALSE)


###Annotate KRT17+KRT5+ cells by batch ID

#Label KRT5+ cells
egfr_epi_t_tn$KRT5_exp <- ifelse(FetchData(egfr_epi_t_tn, vars = "KRT5") > 0, 
                                  "KRT5+", "KRT5-")

head(egfr_epi_t_tn)

#Make combined KRT expression metadata entry
egfr_epi_t_tn$KRT_combination <- paste0(egfr_epi_t_tn$KRT17_exp, egfr_epi_t_tn$KRT5_exp)

#Summarize KRT17 & KRT5 co-expression by batch ID
egfr_epi_t_tn_krt17_krt5_stats <- table(egfr_epi_t_tn@meta.data$batch, egfr_epi_t_tn@meta.data$KRT_combination)
egfr_epi_t_tn_krt17_krt5_stats

#Write data to table
write.table(egfr_epi_t_tn_krt17_krt5_stats, "MDACC EGFR mutant scRNAseq expanded cohort treatment naive KRT17 KRT5 coexpression counts.txt", sep = "\t", row.names = TRUE, quote = FALSE)




##### Expanded cohort osimertinib MRD sample analysis #####

#Read in data
egfr_seurat <- LoadSeuratRds("EGFR_combined_batches.Rds")

#Filter for MRD tumor samples
specimen_ids_to_subset <- c("JH297", "JH143", "Lung-Tumor-1", "Lung-Tumor-7", "Lung-tumor-8",
                            "NSTAR1-TumorA", "JH382", "JH380")

#Subset for Osi MRD only
egfr_epi_omrd <- subset(egfr_seurat, subset = batch %in% specimen_ids_to_subset)

#Subset for epithelial only
egfr_epi_omrd <- subset(egfr_epi_omrd, subset = General_Annot == "EPITHELIAL")

#Subset for InferCnv tumor cells
egfr_epi_t_omrd <- subset(egfr_epi_omrd, subset = InfercnvR == "tumor")


#Summarize the number of epithelial tumor cells by Specimen_ID
egfr_epi_stats <- egfr_epi_t_omrd@meta.data %>%
  group_by(batch) %>%
  summarise(Epithelial_tumor_cells = n())

write.table(egfr_epi_stats, "MDACC EGFR mutant scRNAseq expanded cohort osimertinib MRD sample tumor cell count by batch ID.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Filter for Osi MRD tumor samples with > 20 tumor cells
specimen_ids_to_subset <- c("JH297", "JH143", "Lung-Tumor-1", "Lung-Tumor-7", "Lung-tumor-8",
                            "NSTAR1-TumorA", "JH380")

#Subset for Osi MRD only
egfr_epi_omrd <- subset(egfr_seurat, subset = batch %in% specimen_ids_to_subset)

#Subset for epithelial only
egfr_epi_omrd <- subset(egfr_epi_omrd, subset = General_Annot == "EPITHELIAL")

#Subset for InferCnv tumor cells
egfr_epi_t_omrd <- subset(egfr_epi_omrd, subset = InfercnvR == "tumor")


#Change active identity to batch
Idents(egfr_epi_t_omrd) <- "batch"

#Draw KRT17 expression violin plot across Osi MRD samples
krt17_epi_t_mrd_vp <- VlnPlot(egfr_epi_t_omrd, features = "KRT17", group.by = "batch", assay = "RNA", layer = "data") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH297" = "#00B813",
                               "JH143" = "#00C08E",
                               "Lung-Tumor-1" = "#00C0B4",
                               "Lung-Tumor-7" = "#00BDD4",
                               "Lung-tumor-8" = "#00A7FF",
                               "NSTAR1-TumorA" = "#7F96FF",
                               "JH380" = "#619CFF")) +
  labs(title = "Osi MRD KRT17 expression")

krt17_epi_t_mrd_vp #Plot: 700 x 400


###Annotate KRT17+ vs KRT17- cells by Specimen ID

#Label KRT17+ cells
egfr_epi_t_omrd$KRT17_exp <- ifelse(FetchData(egfr_epi_t_omrd, vars = "KRT17") > 0, 
                                    "KRT17+", "KRT17-")

head(egfr_epi_t_omrd)

#Summarize the number of KRT17+ and KRT17- cells by batch ID
egfr_epi_ormd_krt17_stats <- table(egfr_epi_t_omrd@meta.data$batch, egfr_epi_t_omrd@meta.data$KRT17_exp)
egfr_epi_ormd_krt17_stats

#Write data to file
write.table(egfr_epi_ormd_krt17_stats, "MDACC EGFR mutant scRNAseq expanded cohort osimertinib MRD KRT17 expression stats.txt", sep = "\t", row.names = TRUE, quote = FALSE)


###Annotate KRT17+KRT5+ cells by batch ID

#Label KRT5+ cells
egfr_epi_t_omrd$KRT5_exp <- ifelse(FetchData(egfr_epi_t_omrd, vars = "KRT5") > 0, 
                                 "KRT5+", "KRT5-")

head(egfr_epi_t_omrd)

#Make combined KRT expression metadata entry
egfr_epi_t_omrd$KRT_combination <- paste0(egfr_epi_t_omrd$KRT17_exp, egfr_epi_t_omrd$KRT5_exp)

#Summarize KRT17 & KRT5 co-expression by batch ID
egfr_epi_ormd_krt17_krt5_stats <- table(egfr_epi_t_omrd@meta.data$batch, egfr_epi_t_omrd@meta.data$KRT_combination)
egfr_epi_ormd_krt17_krt5_stats

#Write data to table
write.table(egfr_epi_ormd_krt17_krt5_stats, "MDACC EGFR mutant scRNAseq expanded cohort osi MRD KRT17 KRT5 coexpression counts.txt", sep = "\t", row.names = TRUE, quote = FALSE)



##### Expanded cohort osimertinib progression sample analysis #####

#Read in data
egfr_seurat <- LoadSeuratRds("EGFR_combined_batches.Rds")

#Filter for Osi progression tumor samples 
specimen_ids_to_subset <- c("JH033", "JH038", "JH104", "JH304", 
                            "JH305", "JH339", "JH386")

#Subset for osi progression samples only
egfr_epi_opd <- subset(egfr_seurat, subset = batch %in% specimen_ids_to_subset)

#Subset for epithelial cells only
egfr_epi_opd <- subset(egfr_epi_opd, subset = General_Annot == "EPITHELIAL")

#Subset for InferCNV tumor cells only
egfr_epi_t_opd <- subset(egfr_epi_opd, subset = InfercnvR == "tumor")


#Summarize the number of epithelial tumor cells by Specimen_ID
egfr_epi_stats <- egfr_epi_t_opd@meta.data %>%
  group_by(batch) %>%
  summarise(Epithelial_tumor_cells = n())

write.table(egfr_epi_stats, "MDACC EGFR mutant scRNAseq expanded cohort osimertinib progression sample tumor cell count by batch ID.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Filter for Osi progression tumor samples with > 20 tumor cells
specimen_ids_to_subset <- c("JH033", "JH038", "JH104", 
                            "JH304", "JH305", "JH386")

#Subset for osi progression samples only
egfr_epi_opd <- subset(egfr_seurat, subset = batch %in% specimen_ids_to_subset)

#Subset for epithelial cells only
egfr_epi_opd <- subset(egfr_epi_opd, subset = General_Annot == "EPITHELIAL")

#Subset for InferCNV tumor cells only
egfr_epi_t_opd <- subset(egfr_epi_opd, subset = InfercnvR == "tumor")


#Change active identity to batch
Idents(egfr_epi_t_opd) <- "batch"

#Draw KRT17 expression violin plot across osi progression samples
krt17_epi_t_opd <- VlnPlot(egfr_epi_t_opd, features = "KRT17", group.by = "batch", assay = "RNA", layer = "data") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH304" = "#FF62BF", 
                               "JH305" = "#FF6A9A",
                               "JH386" = "magenta")) +
  labs(title = "Osi PD KRT17 expression")

krt17_epi_t_opd


###Annotate KRT17+ vs KRT17- cells by batch

#Label KRT17+ cells
egfr_epi_t_opd$KRT17_exp <- ifelse(FetchData(egfr_epi_t_opd, vars = "KRT17") > 0, 
                                   "KRT17+", "KRT17-")

head(egfr_epi_t_opd)

#Summarize the number of KRT17+ and KRT17- cells by batch ID
egfr_epi_opd_krt17_stats <- table(egfr_epi_t_opd@meta.data$batch, egfr_epi_t_opd@meta.data$KRT17_exp)
egfr_epi_opd_krt17_stats

#Write data to file
write.table(egfr_epi_opd_krt17_stats, "MDACC EGFR mutant scRNAseq expanded cohort osimertinib PD KRT17 expression stats.txt", sep = "\t", row.names = TRUE, quote = FALSE)


###Annotate KRT17+KRT5+ cells by batch ID

#Label KRT5+ cells
egfr_epi_t_opd$KRT5_exp <- ifelse(FetchData(egfr_epi_t_opd, vars = "KRT5") > 0, 
                                   "KRT5+", "KRT5-")

head(egfr_epi_t_opd)

#Make combined KRT expression metadata entry
egfr_epi_t_opd$KRT_combination <- paste0(egfr_epi_t_opd$KRT17_exp, egfr_epi_t_opd$KRT5_exp)

#Summarize KRT17 & KRT5 co-expression by batch ID
egfr_epi_opd_krt17_krt5_stats <- table(egfr_epi_t_opd@meta.data$batch, egfr_epi_t_opd@meta.data$KRT_combination)
egfr_epi_opd_krt17_krt5_stats

#Write data to table
write.table(egfr_epi_opd_krt17_krt5_stats, "MDACC EGFR mutant scRNAseq expanded cohort osi PD KRT17 KRT5 coexpression counts.txt", sep = "\t", row.names = TRUE, quote = FALSE)



###Draw violin plot for TP63 expression

#Draw TP63 expression violin plot across osi progression samples
tp63_epi_t_opd <- VlnPlot(egfr_epi_t_opd, features = "TP63", group.by = "batch", assay = "RNA", layer = "data") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH304" = "#FF62BF", 
                               "JH305" = "#FF6A9A",
                               "JH386" = "magenta")) +
  labs(title = "Osi PD TP63 expression")

tp63_epi_t_opd #Plot: 700 x 400



###Draw violin plot for ASCL1 expression

#Draw ASCL1 expression violin plot across osi progression samples
ascl1_epi_t_opd <- VlnPlot(egfr_epi_t_opd, features = "ASCL1", group.by = "batch", assay = "RNA", layer = "data") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH304" = "#FF62BF", 
                               "JH305" = "#FF6A9A",
                               "JH386" = "magenta")) +
  labs(title = "Osi PD ASCL1 expression")

ascl1_epi_t_opd #Plot: 700 x 400


###Draw violin plot for MET expression

#Draw MET expression violin plot across osi progression samples
met_epi_t_opd <- VlnPlot(egfr_epi_t_opd, features = "MET", group.by = "batch", assay = "RNA", layer = "data") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH304" = "#FF62BF", 
                               "JH305" = "#FF6A9A",
                               "JH386" = "magenta")) +
  labs(title = "Osi PD MET expression")

met_epi_t_opd #Plot: 700 x 400



###Draw violin plot for EGFR expression

#Draw EGFR expression violin plot across osi progression samples
egfr_exp_epi_t_opd <- VlnPlot(egfr_epi_t_opd, features = "EGFR", group.by = "batch", assay = "RNA", layer = "data") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  theme(text = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = c("JH033" = "#BC81FF", 
                               "JH038" = "#E26EF7",
                               "JH104" = "#F763DF", 
                               "JH304" = "#FF62BF", 
                               "JH305" = "#FF6A9A",
                               "JH386" = "magenta")) +
  labs(title = "Osi PD EGFR expression")

egfr_exp_epi_t_opd #Plot: 700 x 400




##### Expanded cohort treatment naive, osimertinib MRD, and osimertinib progression %KRT17 expressing boxplot #####

#Read in data
tn_mrd_pd_krt17_perc <- MDACC_EGFR_mutant_scRNAseq_expanded_cohort_TN_osi_MRD_osi_PD_percentage_cells_expressing_KRT17

#Facet order
tn_mrd_pd_krt17_perc$Timepoint <- factor(tn_mrd_pd_krt17_perc$Timepoint, levels = c("Treatment_naive", "Osi_MRD", "Osi_progression"))

#Define groups for statistical comparisons
my_comparisons <- list(c("Treatment_naive", "Osi_MRD"),
                       c("Osi_MRD", "Osi_progression"),
                       c("Treatment_naive", "Osi_progression"))

#Draw boxplot
tn_mrd_pd_krt17_perc_p <- ggplot(tn_mrd_pd_krt17_perc, aes(x = tn_mrd_pd_krt17_perc$Timepoint, y = tn_mrd_pd_krt17_perc$Percentage_KRT17pos, fill = tn_mrd_pd_krt17_perc$Timepoint)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  scale_y_continuous("Percentage KRT17+", limits = c(0, 120)) +
  theme_bw() +
  theme(text = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(title = "Percentage KRT17+")


tn_mrd_pd_krt17_perc_p #Plot: 450 x 600



##### Expanded cohort treatment naive, osimertinib MRD, and osimertinib progression %KRT17 expressing barchart #####

#Read in data
tn_mrd_pd_krt17_perc <- MDACC_EGFR_mutant_scRNAseq_expanded_cohort_TN_osi_MRD_osi_PD_percentage_cells_expressing_KRT17

#Set sample order
tn_mrd_pd_krt17_perc$Sample_ID <- factor(tn_mrd_pd_krt17_perc$Sample_ID, levels = c("JH064",
                                                                                    "JH139",
                                                                                    "JH400",
                                                                                    "JH095",
                                                                                    "Biopsy1",
                                                                                    "Lung-tumor-10",
                                                                                    "JH067",
                                                                                    "JH384",
                                                                                    "JH348",
                                                                                    "JH128",
                                                                                    "JH380",
                                                                                    "NSTAR1-TumorA",
                                                                                    "JH297",
                                                                                    "JH143",
                                                                                    "Lung-Tumor-1",
                                                                                    "Lung-Tumor-7",
                                                                                    "Lung-tumor-8",
                                                                                    "JH033",
                                                                                    "JH104",
                                                                                    "JH038",
                                                                                    "JH386",
                                                                                    "JH305",
                                                                                    "JH304"))

#Draw barchart
tn_mrd_pd_krt17_perc_bar <- ggplot(tn_mrd_pd_krt17_perc, aes(x = Sample_ID, y = Percentage_KRT17pos, fill = Timepoint)) +
                            geom_bar(stat = "identity") + 
                            geom_text(aes(label = round(Percentage_KRT17pos, 1)), vjust = -0.5) +
                            theme_bw() + 
                            theme(panel.grid = element_blank(),
                                  panel.background = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text = element_text(size = 12, color = "black"),
                                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                  legend.position = "none") +
                            scale_y_continuous("Percentage KRT17+", limits = c(0, 100)) +
                            scale_fill_manual(values = c("Treatment_naive" = "#F8766D",
                                                         "Osi_MRD" = "#00BA38",
                                                         "Osi_progression" = "#619CFF")) +
                            labs(title = "MDACC expanded cohort percentage KRT17 expressing")


tn_mrd_pd_krt17_perc_bar



##### Expanded cohort treatment naive KRT17 KRT5 co expression barchart #####

#Read in data
tn_krt17_krt5 <- MDACC_EGFR_cohort_treatment_naive_percentage_cells_KRT17_KRT5_co_expressing

#Set sample order
tn_krt17_krt5$Sample_ID <- factor(tn_krt17_krt5$Sample_ID, levels = c("JH095",
                                                          "Biopsy1",
                                                          "Lung-tumor-10",
                                                          "JH067",
                                                          "JH384",
                                                          "JH348",
                                                          "JH128"))

#Plot co-expression barchart
tn_krt17_krt5_bp <- ggplot(tn_krt17_krt5, aes(x = Sample_ID, y = Percentage, fill = Type)) +
                    geom_bar(stat = "identity") +
                    theme_bw() +
                    theme(panel.grid = element_blank(),
                          axis.text = element_text(color = "black", size = 12),
                          axis.title.x = element_blank(),
                          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                    scale_y_continuous("Percentage cells co-expressing") +
                    labs(title = "Treatment naive KRT17 KRT5 coexpression")

tn_krt17_krt5_bp #Plot: 600 x 400



##### Expanded cohort osi MRD KRT17 KRT5 co expression barchart #####

#Read in data
omrd_krt17_krt5 <- MDACC_EGFR_cohort_osi_MRD_percentage_cells_KRT17_KRT5_co_expressing

#Set sample order
omrd_krt17_krt5$Sample_ID <- factor(omrd_krt17_krt5$Sample_ID, levels = c("JH380",
                                                                      "NSTAR1-TumorA",
                                                                      "JH297",
                                                                      "JH143",
                                                                      "Lung-Tumor-1",
                                                                      "Lung-Tumor-7",
                                                                      "Lung-tumor-8"))

#Plot co-expression barchart
omrd_krt17_krt5_bp <- ggplot(omrd_krt17_krt5, aes(x = Sample_ID, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous("Percentage cells co-expressing") +
  labs(title = "Osi MRD KRT17 KRT5 coexpression")

omrd_krt17_krt5_bp #Plot: 600 x 400



##### Expanded cohort osi PD KRT17 KRT5 co expression barchart #####

#Read in data
opd_krt17_krt5 <- MDACC_EGFR_cohort_osi_progression_percentage_cells_KRT17_KRT5_co_expressing

#Set sample order
opd_krt17_krt5$Sample_ID <- factor(opd_krt17_krt5$Sample_ID, levels = c("JH033",
                                                                        "JH038",
                                                                          "JH104",
                                                                          "JH386",
                                                                          "JH305",
                                                                          "JH304"))

#Plot co-expression barchart
opd_krt17_krt5_bp <- ggplot(opd_krt17_krt5, aes(x = Sample_ID, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous("Percentage cells co-expressing") +
  labs(title = "Osi Progression KRT17 KRT5 coexpression")

opd_krt17_krt5_bp #Plot: 600 x 400











