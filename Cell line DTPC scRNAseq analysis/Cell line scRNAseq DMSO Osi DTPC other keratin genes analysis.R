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
library(paletteer)
library(scico)
library(ggpubr)


##### HCC827 DMSO and HCC827 Osi DTPC kertain analysis #####

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

###KRT6A

#Define genes of interest
genes_of_interest <- c("KRT6A")

#Extract expression data
krt6a_expression_data <- FetchData(hcc827_combined, vars = genes_of_interest)

#Convert to dataframe
krt6a_expression_data <- as.data.frame(krt6a_expression_data)

#Define thresholds for gene expression
threshold <- 0.1
krt6a_expression_data <- krt6a_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT6A > threshold, 1, 0),
    KRT6A_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc827_combined <- AddMetaData(hcc827_combined, krt6a_expression_data$KRT6A_ExpressionCategory, col.name = "KRT6A_ExpressionCategory")

#Build custom color key for KRT6A_ExpressionCategory
existing_column <- "KRT6A_ExpressionCategory"

label_map <- c(
  "0" = "KRT6A-",
  "1" = "KRT6A+")


#Assign KRT6A_ExpressionCategory data to separate metadata variable
metadata <- hcc827_combined@meta.data

#Match labels to ExpressionCategory binary values using key
metadata$KRT6A_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                   label_map[metadata[[existing_column]]], 
                                   "Other")  # Default "Other" if value doesn't match

#Add new KRT6A labels back to metadata
hcc827_combined <- AddMetaData(hcc827_combined, metadata$KRT6A_ExpressionLabel, col.name = "KRT6A_ExpressionLabel")

#Check new metadata entry
head(hcc827_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT6A-" = "gray",
            "KRT6A+" = "red")

#Plot UMAP with KRT6A pseudocoloring
hcc827_krt6a_dim <- DimPlot(hcc827_combined, group.by = "KRT6A_ExpressionLabel", 
                            split.by = "orig.ident",
                            pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT6A expression",
       color = "KRT6A") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

hcc827_krt6a_dim #Plot: 800 x 400

#Generate table for number of cells expressing KRT6A by condition
krt6a_exp_stats <- table(hcc827_combined$KRT6A_ExpressionLabel, hcc827_combined$orig.ident)
krt6a_exp_stats

#Write table to file
write.table(krt6a_exp_stats, "HCC827 DMSO Osi DTPC KRT6A expression statistics.txt", sep = "\t", quote = FALSE)



###KRT14

#Define genes of interest
genes_of_interest <- c("KRT14")

#Extract expression data
KRT14_expression_data <- FetchData(hcc827_combined, vars = genes_of_interest)

#Convert to dataframe
KRT14_expression_data <- as.data.frame(KRT14_expression_data)

#Define thresholds for gene expression
threshold <- 0.1
KRT14_expression_data <- KRT14_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT14 > threshold, 1, 0),
    KRT14_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc827_combined <- AddMetaData(hcc827_combined, KRT14_expression_data$KRT14_ExpressionCategory, col.name = "KRT14_ExpressionCategory")

#Build custom color key for KRT14_ExpressionCategory
existing_column <- "KRT14_ExpressionCategory"

label_map <- c(
  "0" = "KRT14-",
  "1" = "KRT14+")


#Assign KRT14_ExpressionCategory data to separate metadata variable
metadata <- hcc827_combined@meta.data

#Match labels to ExpressionCategory binary values using key
metadata$KRT14_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT14 labels back to metadata
hcc827_combined <- AddMetaData(hcc827_combined, metadata$KRT14_ExpressionLabel, col.name = "KRT14_ExpressionLabel")

#Check new metadata entry
head(hcc827_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT14-" = "gray",
            "KRT14+" = "red")

#Plot UMAP with KRT14 pseudocoloring
hcc827_KRT14_dim <- DimPlot(hcc827_combined, group.by = "KRT14_ExpressionLabel", 
                            split.by = "orig.ident",
                            pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT14 expression",
       color = "KRT14") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

hcc827_KRT14_dim #Plot: 800 x 400

#Generate table for number of cells expressing KRT14 by condition
KRT14_exp_stats <- table(hcc827_combined$KRT14_ExpressionLabel, hcc827_combined$orig.ident)
KRT14_exp_stats

#Write table to file
write.table(KRT14_exp_stats, "HCC827 DMSO Osi DTPC KRT14 expression statistics.txt", sep = "\t", quote = FALSE)


###KRT15

#Define genes of interest
genes_of_interest <- c("KRT15")

#Extract expression data
KRT15_expression_data <- FetchData(hcc827_combined, vars = genes_of_interest)

#Convert to dataframe
KRT15_expression_data <- as.data.frame(KRT15_expression_data)

#Define thresholds for gene expression
threshold <- 0.1
KRT15_expression_data <- KRT15_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT15 > threshold, 1, 0),
    KRT15_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc827_combined <- AddMetaData(hcc827_combined, KRT15_expression_data$KRT15_ExpressionCategory, col.name = "KRT15_ExpressionCategory")

#Build custom color key for KRT15_ExpressionCategory
existing_column <- "KRT15_ExpressionCategory"

label_map <- c(
  "0" = "KRT15-",
  "1" = "KRT15+")


#Assign KRT15_ExpressionCategory data to separate metadata variable
metadata <- hcc827_combined@meta.data

#Match labels to ExpressionCategory binary values using key
metadata$KRT15_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT15 labels back to metadata
hcc827_combined <- AddMetaData(hcc827_combined, metadata$KRT15_ExpressionLabel, col.name = "KRT15_ExpressionLabel")

#Check new metadata entry
head(hcc827_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT15-" = "gray",
            "KRT15+" = "red")

#Plot UMAP with KRT15 pseudocoloring
hcc827_KRT15_dim <- DimPlot(hcc827_combined, group.by = "KRT15_ExpressionLabel", 
                            split.by = "orig.ident",
                            pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT15 expression",
       color = "KRT15") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

hcc827_KRT15_dim #Plot: 800 x 400

#Generate table for number of cells expressing KRT15 by condition
KRT15_exp_stats <- table(hcc827_combined$KRT15_ExpressionLabel, hcc827_combined$orig.ident)
KRT15_exp_stats

#Write table to file
write.table(KRT15_exp_stats, "HCC827 DMSO Osi DTPC KRT15 expression statistics.txt", sep = "\t", quote = FALSE)


###KRT17

#Define genes of interest
genes_of_interest <- c("KRT17")

#Extract expression data
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

#Match labels to ExpressionCategory binary values using key
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


#Generate table for number of cells expressing KRT17 by condition
KRT17_exp_stats <- table(hcc827_combined$KRT17_ExpressionLabel, hcc827_combined$orig.ident)
KRT17_exp_stats

#Write table to file
write.table(KRT17_exp_stats, "HCC827 DMSO Osi DTPC KRT17 expression statistics.txt", sep = "\t", quote = FALSE)



###KRT5

#Define genes of interest
genes_of_interest <- c("KRT5")

#Extract expression data
KRT5_expression_data <- FetchData(hcc827_combined, vars = genes_of_interest) ##Not found in SCT assay

#Check if non-zero KRT5 expression exists
DefaultAssay(hcc827_combined) <- "RNA"

krt5_expr <- FetchData(hcc827_combined, vars = "KRT5", assay = "RNA", layer = "counts")
summary(krt5_expr$KRT5) #All cells lack KRT5 expression


###KRT13

#Define genes of interest
genes_of_interest <- c("KRT13")

DefaultAssay(hcc827_combined) <- "SCT"

#Extract expression data
krt13_expression_data <- FetchData(hcc827_combined, vars = genes_of_interest)

#Convert to dataframe
krt13_expression_data <- as.data.frame(krt13_expression_data)

#Define thresholds for gene expression (adjust as needed)
threshold <- 0.1
krt13_expression_data <- krt13_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT13 > threshold, 1, 0),
    KRT13_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc827_combined <- AddMetaData(hcc827_combined, krt13_expression_data$KRT13_ExpressionCategory, col.name = "KRT13_ExpressionCategory")

#Build custom color key for KRT13_ExpressionCategory
existing_column <- "KRT13_ExpressionCategory"

label_map <- c(
  "0" = "KRT13-",
  "1" = "KRT13+")


#Assign KRT13_ExpressionCategory data to separate metadata variable
metadata <- hcc827_combined@meta.data

#Match labels to ExpressionCategory binary values using key
metadata$KRT13_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT13 labels back to metadata
hcc827_combined <- AddMetaData(hcc827_combined, metadata$KRT13_ExpressionLabel, col.name = "KRT13_ExpressionLabel")

#Check new metadata entry
head(hcc827_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT13-" = "gray",
            "KRT13+" = "red")

#Plot UMAP with KRT13 pseudocoloring
hcc827_krt13_dim <- DimPlot(hcc827_combined, group.by = "KRT13_ExpressionLabel", 
                            split.by = "orig.ident",
                            pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT13 expression",
       color = "KRT13") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

hcc827_krt13_dim #Plot: 800 x 400


#Generate table for number of cells expressing KRT13 by condition
KRT13_exp_stats <- table(hcc827_combined$KRT13_ExpressionLabel, hcc827_combined$orig.ident)
KRT13_exp_stats

#Write table to file
write.table(KRT13_exp_stats, "HCC827 DMSO Osi DTPC KRT13 expression statistics.txt", sep = "\t", quote = FALSE)


###KRT combination expression patterns

#Make combined KRT expression metadata entry
hcc827_combined$KRT_combination <- paste0(hcc827_combined$KRT6A_ExpressionLabel, hcc827_combined$KRT14_ExpressionLabel, hcc827_combined$KRT15_ExpressionLabel, hcc827_combined$KRT17_ExpressionLabel, hcc827_combined$KRT13_ExpressionLabel)

#Create table summarizing KRT expression patterns across conditions
KRT_combo_exp_stats <- table(hcc827_combined$KRT_combination, hcc827_combined$orig.ident)
KRT_combo_exp_stats

#Write table to file
write.table(KRT_combo_exp_stats, "HCC827 DMSO Osi DTPC keratin coexpression statistics.txt", sep = "\t", quote = FALSE)



##### HCC827 DMSO and HCC827 Osi DTPC KRT6A KRT13 KRT14 KRT15 KRT5 percent expressing barchart #####

#Read in data
hcc827_krt_exp <- HCC827_DMSO_Osi_DTPC_KRT6A_KRT13_KRT14_KRT15_KRT5_percent_expressing

#Set category order
hcc827_krt_exp$Category <- factor(hcc827_krt_exp$Category, levels = c("KRT6A-", "KRT6A+", "KRT13-", "KRT13+", "KRT14-", "KRT14+", "KRT15-", "KRT15+", "KRT5-", "KRT5+"))

#Plot barchart
hcc827_krt_exp_bp <- ggplot(hcc827_krt_exp, aes(x = Category, y = Percentage, fill = Sample_ID)) + 
                     geom_bar(stat = "identity", position = "dodge") + 
                     theme_bw() +
                     theme(panel.grid = element_blank(),
                           panel.background = element_blank(),
                           axis.text = element_text(size = 12, color = "black"), 
                           axis.title.x = element_blank(),
                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                     scale_fill_manual(values = c("HCC827_DMSO" = "lightblue",
                                                  "HCC827_Osi_DTPC" = "blue")) +
                     labs(title = "HCC827 DMSO Osi DTPC Keratin expression")

hcc827_krt_exp_bp #Plot: 800 x 400



##### HCC827 DMSO and HCC827 Osi DTPC KRT17 KRT6A KRT13 KRT14 KRT15 KRT5 percent coexpression barchart #####

#Read in data
hcc827_krt_coexp <- HCC827_DMSO_Osi_DTPC_KRT17_KRT6A_KRT13_KRT14_KRT15_KRT5_percent_coexpression

#Plot barchart
hcc827_krt_coexp_bp <- ggplot(hcc827_krt_coexp, aes(x = Category, y = Percentage, fill = Sample_ID)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("HCC827_DMSO" = "lightblue",
                               "HCC827_Osi_DTPC" = "blue")) +
  labs(title = "HCC827 DMSO Osi DTPC Keratin coexpression")

hcc827_krt_coexp_bp #Plot: 800 x 800




##### HCC4006 DMSO and HCC4006 Osi DTPC kertain analysis #####

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

###KRT6A

#Define genes of interest
genes_of_interest <- c("KRT6A")

#Extract expression data
krt6a_expression_data <- FetchData(hcc4006_combined, vars = genes_of_interest) ###No KRT6A expression found in SCT


#Check if non-zero KRT5 expression data exists
DefaultAssay(hcc4006_combined) <- "RNA"

krt6A_expr <- FetchData(hcc4006_combined, vars = "KRT6A", assay = "RNA", layer = "counts")
summary(krt6A_expr$KRT6A) ##No KRT6A expression detected



###KRT14

#Define genes of interest
genes_of_interest <- c("KRT14")

DefaultAssay(hcc4006_combined) <- "SCT"

#Extract expression data
KRT14_expression_data <- FetchData(hcc4006_combined, vars = genes_of_interest)

#Convert to dataframe
KRT14_expression_data <- as.data.frame(KRT14_expression_data)

#Define thresholds for gene expression
threshold <- 0.1
KRT14_expression_data <- KRT14_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT14 > threshold, 1, 0),
    KRT14_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc4006_combined <- AddMetaData(hcc4006_combined, KRT14_expression_data$KRT14_ExpressionCategory, col.name = "KRT14_ExpressionCategory")

#Build custom color key for KRT14_ExpressionCategory
existing_column <- "KRT14_ExpressionCategory"

label_map <- c(
  "0" = "KRT14-",
  "1" = "KRT14+")


#Assign KRT14_ExpressionCategory data to separate metadata variable
metadata <- hcc4006_combined@meta.data

#Match labels to ExpressionCategory binary values using key
metadata$KRT14_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT14 labels back to metadata
hcc4006_combined <- AddMetaData(hcc4006_combined, metadata$KRT14_ExpressionLabel, col.name = "KRT14_ExpressionLabel")

#Check new metadata entry
head(hcc4006_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT14-" = "gray",
            "KRT14+" = "red")

#Plot UMAP with KRT14 pseudocoloring
hcc4006_KRT14_dim <- DimPlot(hcc4006_combined, group.by = "KRT14_ExpressionLabel", 
                            split.by = "orig.ident",
                            pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT14 expression",
       color = "KRT14") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

hcc4006_KRT14_dim #Plot: 800 x 400

#Generate table for number of cells expressing KRT14 by condition
KRT14_exp_stats <- table(hcc4006_combined$KRT14_ExpressionLabel, hcc4006_combined$orig.ident)
KRT14_exp_stats

#Write table to file
write.table(KRT14_exp_stats, "HCC4006 DMSO Osi DTPC KRT14 expression statistics.txt", sep = "\t", quote = FALSE)



###KRT15

#Define genes of interest
genes_of_interest <- c("KRT15")

#Extract expression data
KRT15_expression_data <- FetchData(hcc4006_combined, vars = genes_of_interest)

#Convert to dataframe
KRT15_expression_data <- as.data.frame(KRT15_expression_data)

#Define thresholds for gene expression
threshold <- 0.1
KRT15_expression_data <- KRT15_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT15 > threshold, 1, 0),
    KRT15_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
hcc4006_combined <- AddMetaData(hcc4006_combined, KRT15_expression_data$KRT15_ExpressionCategory, col.name = "KRT15_ExpressionCategory")

#Build custom color key for KRT15_ExpressionCategory
existing_column <- "KRT15_ExpressionCategory"

label_map <- c(
  "0" = "KRT15-",
  "1" = "KRT15+")


#Assign KRT15_ExpressionCategory data to separate metadata variable
metadata <- hcc4006_combined@meta.data

#Match labels to ExpressionCategory binary values using key
metadata$KRT15_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT15 labels back to metadata
hcc4006_combined <- AddMetaData(hcc4006_combined, metadata$KRT15_ExpressionLabel, col.name = "KRT15_ExpressionLabel")

#Check new metadata entry
head(hcc4006_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT15-" = "gray",
            "KRT15+" = "red")

#Plot UMAP with KRT15 pseudocoloring
hcc4006_KRT15_dim <- DimPlot(hcc4006_combined, group.by = "KRT15_ExpressionLabel", 
                            split.by = "orig.ident",
                            pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT15 expression",
       color = "KRT15") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

hcc4006_KRT15_dim #Plot: 800 x 400

#Generate table for number of cells expressing KRT15 by condition
KRT15_exp_stats <- table(hcc4006_combined$KRT15_ExpressionLabel, hcc4006_combined$orig.ident)
KRT15_exp_stats

#Write table to file
write.table(KRT15_exp_stats, "HCC4006 DMSO Osi DTPC KRT15 expression statistics.txt", sep = "\t", quote = FALSE)


###KRT17

#Define genes of interest
genes_of_interest <- c("KRT17")

#Extract expression data
KRT17_expression_data <- FetchData(hcc4006_combined, vars = genes_of_interest)

#Convert to dataframe
KRT17_expression_data <- as.data.frame(KRT17_expression_data)

#Define thresholds for gene expression
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

#Match labels to ExpressionCategory binary values using key
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

#Generate table for number of cells expressing KRT17 by condition
KRT17_exp_stats <- table(hcc4006_combined$KRT17_ExpressionLabel, hcc4006_combined$orig.ident)
KRT17_exp_stats

#Write table to file
write.table(KRT17_exp_stats, "HCC4006 DMSO Osi DTPC KRT17 expression statistics.txt", sep = "\t", quote = FALSE)


###KRT5

#Define genes of interest
genes_of_interest <- c("KRT5")

#Extract expression data
KRT5_expression_data <- FetchData(hcc4006_combined, vars = genes_of_interest) #Not found in SCT

#Check if non-zero KRT5 expression exists
DefaultAssay(hcc4006_combined) <- "RNA"

krt5_expr <- FetchData(hcc4006_combined, vars = "KRT5", assay = "RNA", layer = "counts")
summary(krt5_expr$KRT5) ##Rare KRT5 expression detected


#Define thresholds for gene expression
KRT5_expression_data <- ifelse(krt5_expr$KRT5 > 0, 
                               "KRT5+", "KRT5-")

hcc4006_combined$KRT5_ExpressionCategory <- KRT5_expression_data

#Generate table for number of cells expressing KRT5 by condition
KRT5_exp_stats <- table(hcc4006_combined$KRT5_ExpressionCategory, hcc4006_combined$orig.ident)
KRT5_exp_stats

#Write table to file
write.table(KRT5_exp_stats, "HCC4006 DMSO Osi DTPC KRT5 expression statistics.txt", sep = "\t", quote = FALSE)


###KRT13

#Define genes of interest
genes_of_interest <- c("KRT13")

DefaultAssay(hcc4006_combined) <- "SCT"

#Extract expression data
krt13_expression_data <- FetchData(hcc4006_combined, vars = genes_of_interest) #No expression in SCT

#Check if non-zero KRTt13 expression exists
DefaultAssay(hcc4006_combined) <- "RNA"

krt13_expr <- FetchData(hcc4006_combined, vars = "KRT13", assay = "RNA", layer = "counts")
summary(krt13_expr$KRT13) ##Rare KRT13 expression detected


#Define thresholds for gene expression
krt13_expression_data <- ifelse(krt13_expr$KRT13 > 0, 
                               "KRT13+", "KRTt13-")

hcc4006_combined$KRT13_ExpressionCategory <- krt13_expression_data


#Generate table for number of cells expressing KRT13 by condition
KRT13_exp_stats <- table(hcc4006_combined$KRT13_ExpressionCategory, hcc4006_combined$orig.ident)
KRT13_exp_stats

#Write table to file
write.table(KRT13_exp_stats, "HCC4006 DMSO Osi DTPC KRT13 expression statistics.txt", sep = "\t", quote = FALSE)




###Combination expression pattern

#Make combined KRT expression metadata entry
hcc4006_combined$KRT_combination <- paste0(hcc4006_combined$KRT14_ExpressionLabel, hcc4006_combined$KRT15_ExpressionLabel, hcc4006_combined$KRT17_ExpressionLabel, hcc4006_combined$KRT5_ExpressionCategory, hcc4006_combined$KRT13_ExpressionCategory)

#Create table summarizing KRT expression patterns across conditions
KRT_combo_exp_stats <- table(hcc4006_combined$KRT_combination, hcc4006_combined$orig.ident)
KRT_combo_exp_stats

#Write table to file
write.table(KRT_combo_exp_stats, "HCC4006 DMSO Osi DTPC keratin coexpression statistics.txt", sep = "\t", quote = FALSE)



##### HCC4006 DMSO and HCC4006 Osi DTPC KRT6A KRT13 KRT14 KRT15 KRT5 percent expressing barchart #####

#Read in data
hcc4006_krt_exp <- HCC4006_DMSO_Osi_DTPC_KRT6A_KRT13_KRT14_KRT15_KRT5_percent_expressing

#Set category order
hcc4006_krt_exp$Category <- factor(hcc4006_krt_exp$Category, levels = c("KRT6A-", "KRT6A+", "KRT13-", "KRT13+", "KRT14-", "KRT14+", "KRT15-", "KRT15+", "KRT5-", "KRT5+"))

#Plot barchart
hcc4006_krt_exp_bp <- ggplot(hcc4006_krt_exp, aes(x = Category, y = Percentage, fill = Sample_ID)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("HCC4006_DMSO" = "palegreen",
                               "HCC4006_Osi_DTPC" = "green3")) +
  labs(title = "HCC4006 DMSO Osi DTPC Keratin expression")

hcc4006_krt_exp_bp #Plot: 800 x 400



##### HCC4006 DMSO and HCC4006 Osi DTPC KRT17 KRT6A KRT13 KRT14 KRT15 KRT5 percent coexpression barchart #####

#Read in data
hcc4006_krt_coexp <- HCC4006_DMSO_Osi_DTPC_KRT17_KRT6A_KRT13_KRT14_KRT15_KRT5_percent_coexpressing

#Plot barchart
hcc4006_krt_coexp_bp <- ggplot(hcc4006_krt_coexp, aes(x = Category, y = Percentage, fill = Sample_ID)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("HCC4006_DMSO" = "palegreen",
                               "HCC4006_Osi_DTPC" = "green3")) +
  labs(title = "HCC4006 DMSO Osi DTPC Keratin coexpression")

hcc4006_krt_coexp_bp #Plot: 1000 x 800



##### H1975 DMSO and H1975 Osi DTPC kertain analysis #####

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

###KRT6A

#Define genes of interest
genes_of_interest <- c("KRT6A")

#Extract expression data
krt6a_expression_data <- FetchData(h1975_combined, vars = genes_of_interest) ###No KRT6A expression in SCT

#Check if non-zero KRT5 expression data exists
DefaultAssay(h1975_combined) <- "RNA"

krt6A_expr <- FetchData(h1975_combined, vars = "KRT6A", assay = "RNA", layer = "counts")
summary(krt6A_expr$KRT6A) ##Rare KRT6A expression detected


#Define thresholds for gene expression
KRT6A_expression_data <- ifelse(krt6A_expr$KRT6A > 0, 
                               "KRT6A+", "KRT6A-")

h1975_combined$KRT6A_ExpressionCategory <- KRT6A_expression_data

#Generate table for number of cells expressing KRT6A by condition
KRT6A_exp_stats <- table(h1975_combined$KRT6A_ExpressionCategory, h1975_combined$orig.ident)
KRT6A_exp_stats

#Write table to file
write.table(KRT6A_exp_stats, "H1975 DMSO Osi DTPC KRT6A expression statistics.txt", sep = "\t", quote = FALSE)



###KRT14

#Define genes of interest
genes_of_interest <- c("KRT14")

DefaultAssay(h1975_combined) <- "SCT"

#Extract expression data
KRT14_expression_data <- FetchData(h1975_combined, vars = genes_of_interest)

#Convert to dataframe
KRT14_expression_data <- as.data.frame(KRT14_expression_data)

#Define thresholds for gene expression
threshold <- 0.1
KRT14_expression_data <- KRT14_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT14 > threshold, 1, 0),
    KRT14_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
h1975_combined <- AddMetaData(h1975_combined, KRT14_expression_data$KRT14_ExpressionCategory, col.name = "KRT14_ExpressionCategory")

#Build custom color key for KRT14_ExpressionCategory
existing_column <- "KRT14_ExpressionCategory"

label_map <- c(
  "0" = "KRT14-",
  "1" = "KRT14+")


#Assign KRT14_ExpressionCategory data to separate metadata variable
metadata <- h1975_combined@meta.data

#Match labels to ExpressionCategory binary values using key
metadata$KRT14_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT14 labels back to metadata
h1975_combined <- AddMetaData(h1975_combined, metadata$KRT14_ExpressionLabel, col.name = "KRT14_ExpressionLabel")

#Check new metadata entry
head(h1975_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT14-" = "gray",
            "KRT14+" = "red")

#Plot UMAP with KRT14 pseudocoloring
h1975_KRT14_dim <- DimPlot(h1975_combined, group.by = "KRT14_ExpressionLabel", 
                            split.by = "orig.ident",
                            pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT14 expression",
       color = "KRT14") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

h1975_KRT14_dim #Plot: 800 x 400

#Generate table for number of cells expressing KRT14 by condition
KRT14_exp_stats <- table(h1975_combined$KRT14_ExpressionLabel, h1975_combined$orig.ident)
KRT14_exp_stats

#Write table to file
write.table(KRT14_exp_stats, "H1975 DMSO Osi DTPC KRT14 expression statistics.txt", sep = "\t", quote = FALSE)



###KRT15

#Define genes of interest
genes_of_interest <- c("KRT15")

#Extract expression data
KRT15_expression_data <- FetchData(h1975_combined, vars = genes_of_interest)

#Convert to dataframe
KRT15_expression_data <- as.data.frame(KRT15_expression_data)

#Define thresholds for gene expression
threshold <- 0.1
KRT15_expression_data <- KRT15_expression_data %>%
  mutate(
    Gene1_expr = ifelse(KRT15 > threshold, 1, 0),
    KRT15_ExpressionCategory = paste0(Gene1_expr))

#Add ExpressionCategory to metadata
h1975_combined <- AddMetaData(h1975_combined, KRT15_expression_data$KRT15_ExpressionCategory, col.name = "KRT15_ExpressionCategory")

#Build custom color key for KRT15_ExpressionCategory
existing_column <- "KRT15_ExpressionCategory"

label_map <- c(
  "0" = "KRT15-",
  "1" = "KRT15+")


#Assign KRT15_ExpressionCategory data to separate metadata variable
metadata <- h1975_combined@meta.data

#Match labels to ExpressionCategory binary values using key
metadata$KRT15_ExpressionLabel <- ifelse(metadata[[existing_column]] %in% names(label_map), 
                                         label_map[metadata[[existing_column]]], 
                                         "Other")  # Default "Other" if value doesn't match

#Add new KRT15 labels back to metadata
h1975_combined <- AddMetaData(h1975_combined, metadata$KRT15_ExpressionLabel, col.name = "KRT15_ExpressionLabel")

#Check new metadata entry
head(h1975_combined@meta.data)

#Define custom color mapping (required for proper visualization)
colors <- c("KRT15-" = "gray",
            "KRT15+" = "red")

#Plot UMAP with KRT15 pseudocoloring
h1975_KRT15_dim <- DimPlot(h1975_combined, group.by = "KRT15_ExpressionLabel", 
                            split.by = "orig.ident",
                            pt.size = 1) +
  scale_color_manual(values = colors) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(title = "KRT15 expression",
       color = "KRT15") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        strip.text.x = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

h1975_KRT15_dim #Plot: 800 x 400

#Generate table for number of cells expressing KRT15 by condition
KRT15_exp_stats <- table(h1975_combined$KRT15_ExpressionLabel, h1975_combined$orig.ident)
KRT15_exp_stats

#Write table to file
write.table(KRT15_exp_stats, "H1975 DMSO Osi DTPC KRT15 expression statistics.txt", sep = "\t", quote = FALSE)


###KRT17

#Define genes of interest
genes_of_interest <- c("KRT17")

#Extract expression data
KRT17_expression_data <- FetchData(h1975_combined, vars = genes_of_interest)

#Convert to dataframe
KRT17_expression_data <- as.data.frame(KRT17_expression_data)

#Define thresholds for gene expression
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

#Match labels to ExpressionCategory binary values using key
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

#Generate table for number of cells expressing KRT17 by condition
KRT17_exp_stats <- table(h1975_combined$KRT17_ExpressionLabel, h1975_combined$orig.ident)
KRT17_exp_stats

#Write table to file
write.table(KRT17_exp_stats, "H1975 DMSO Osi DTPC KRT17 expression statistics.txt", sep = "\t", quote = FALSE)


###KRT5

#Define genes of interest
genes_of_interest <- c("KRT5")

#Extract expression data
KRT5_expression_data <- FetchData(h1975_combined, vars = genes_of_interest) #Not found in SCT

#Check if non-zero KRT5 expression data exists
DefaultAssay(h1975_combined) <- "RNA"

krt5_expr <- FetchData(h1975_combined, vars = "KRT5", assay = "RNA", layer = "counts")
summary(krt5_expr$KRT5) ##Rare KRT5 expression detected


#Define thresholds for gene expression (adjust as needed)
KRT5_expression_data <- ifelse(krt5_expr$KRT5 > 0, 
                               "KRT5+", "KRT5-")

h1975_combined$KRT5_ExpressionCategory <- KRT5_expression_data

#Generate table for number of cells expressing KRT5 by condition
KRT5_exp_stats <- table(h1975_combined$KRT5_ExpressionCategory, h1975_combined$orig.ident)
KRT5_exp_stats

#Write table to file
write.table(KRT5_exp_stats, "H1975 DMSO Osi DTPC KRT5 expression statistics.txt", sep = "\t", quote = FALSE)


###KRT13

#Define genes of interest
genes_of_interest <- c("KRT13")

#Extract expression data
KRT13_expression_data <- FetchData(h1975_combined, vars = genes_of_interest) #No SCT expression

#Check if non-zero KRT13 expression data exists
DefaultAssay(h1975_combined) <- "RNA"

krt13_expr <- FetchData(h1975_combined, vars = "KRT13", assay = "RNA", layer = "counts")
summary(krt13_expr$KRT5) ##No KRT13 expression detected


###Combination expression pattern

#Make combined KRT expression metadata entry
h1975_combined$KRT_combination <- paste0(h1975_combined$KRT6A_ExpressionCategory, h1975_combined$KRT14_ExpressionLabel, h1975_combined$KRT15_ExpressionLabel, h1975_combined$KRT17_ExpressionLabel, h1975_combined$KRT5_ExpressionCategory)

#Create table summarizing KRT expression patterns across conditions
KRT_combo_exp_stats <- table(h1975_combined$KRT_combination, h1975_combined$orig.ident)
KRT_combo_exp_stats

#Write table to file
write.table(KRT_combo_exp_stats, "H1975 DMSO Osi DTPC keratin coexpression statistics.txt", sep = "\t", quote = FALSE)



##### H1975 DMSO and H1975 Osi DTPC KRT6A KRT13 KRT14 KRT15 KRT5 percent expressing barchart #####

#Read in data
h1975_krt_exp <- H1975_DMSO_Osi_DTPC_KRT6A_KRT13_KRT14_KRT15_KRT5_percent_expressing

#Set category order
h1975_krt_exp$Category <- factor(h1975_krt_exp$Category, levels = c("KRT6A-", "KRT6A+", "KRT13-", "KRT13+", "KRT14-", "KRT14+", "KRT15-", "KRT15+", "KRT5-", "KRT5+"))

#Plot barchart
h1975_krt_exp_bp <- ggplot(h1975_krt_exp, aes(x = Category, y = Percentage, fill = Sample_ID)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("H1975_DMSO" = "indianred1",
                               "H1975_Osi_DTPC" = "red3")) +
  labs(title = "H1975 DMSO Osi DTPC Keratin expression")

h1975_krt_exp_bp #Plot: 800 x 400



##### H1975 DMSO and H1975 Osi DTPC KRT17 KRT6A KRT13 KRT14 KRT15 KRT5 percent coexpression barchart #####

#Read in data
h1975_krt_coexp <- H1975_DMSO_Osi_DTPC_KRT17_KRT6A_KRT13_KRT14_KRT15_KRT5_percent_coexpressing

#Plot barchart
h1975_krt_coexp_bp <- ggplot(h1975_krt_coexp, aes(x = Category, y = Percentage, fill = Sample_ID)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("H1975_DMSO" = "indianred1",
                               "H1975_Osi_DTPC" = "red")) +
  labs(title = "H1975 DMSO Osi DTPC Keratin coexpression")

h1975_krt_coexp_bp #Plot: 1200 x 800







