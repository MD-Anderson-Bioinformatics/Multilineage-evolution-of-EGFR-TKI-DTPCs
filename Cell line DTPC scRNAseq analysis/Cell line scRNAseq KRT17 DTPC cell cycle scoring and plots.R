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



##### HCC827 Osi DTPC Lito cell cycle by KRT17 status #####

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


#Create table with KRT17+ & KRT17- DTPCs by Lito cell cycle state
krt17_lito_cc <- as.data.frame(table(hcc827_dtpc@meta.data$ExpressionLabel, 
                                     hcc827_dtpc@meta.data$cell_cycle_expanded))

#Add column names
colnames(krt17_lito_cc) <- c("KRT17_status", "Lito_CC", "Number_of_cells")

krt17_lito_cc

#Calculate percentage of cells occupying each cell cycle state by KRT17 status
cell_counts_pct <- krt17_lito_cc %>%
  group_by(KRT17_status) %>%
  mutate(total_per_krt17 = sum(Number_of_cells),
         pct = (Number_of_cells / total_per_krt17) * 100)

#Round percentage to 2 decimals for plotting
cell_counts_pct$pct <- round(cell_counts_pct$pct, 2)

#Draw stacked barchart for KRT17+ vs KRT17- Osi DTPC Lito CC states
krt17_lito_cc_bp <- ggplot(cell_counts_pct, aes(x = Lito_CC, y = pct, fill = Lito_CC)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = pct), vjust = -1) +
  facet_wrap(~ KRT17_status) +
  scale_y_continuous("Percentage of cells", limits = c(0, 100)) +
  scale_fill_manual(values = c("G1S" = "lightgreen", 
                               "S" = "blue", 
                               "G2M" = "purple",
                               "G0" = "red",
                               "M" = "magenta",
                               "MG1" = "gold")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_blank()) +
  labs(title = "HCC827 Osi DTPCs")


krt17_lito_cc_bp #Plot: 600 x 350



##### HCC4006 Osi DTPC Lito cell cycle by KRT17 status #####

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


#Create table with KRT17+ & KRT17- DTPCs by Lito cell cycle state
krt17_lito_cc <- as.data.frame(table(hcc4006_dtpc@meta.data$ExpressionLabel, 
                                     hcc4006_dtpc@meta.data$cell_cycle_expanded))

#Add column names
colnames(krt17_lito_cc) <- c("KRT17_status", "Lito_CC", "Number_of_cells")

krt17_lito_cc

#Calculate percentage of cells occupying each cell cycle state by KRT17 status
cell_counts_pct <- krt17_lito_cc %>%
  group_by(KRT17_status) %>%
  mutate(total_per_krt17 = sum(Number_of_cells),
         pct = (Number_of_cells / total_per_krt17) * 100)

#Round percentage to 2 decimals for plotting
cell_counts_pct$pct <- round(cell_counts_pct$pct, 2)

#Draw stacked barchart for KRT17+ vs KRT17- Osi DTPC Lito CC states
krt17_lito_cc_bp <- ggplot(cell_counts_pct, aes(x = Lito_CC, y = pct, fill = Lito_CC)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = pct), vjust = -1) +
  facet_wrap(~ KRT17_status) +
  scale_y_continuous("Percentage of cells", limits = c(0, 100)) +
  scale_fill_manual(values = c("G1S" = "lightgreen", 
                               "S" = "blue", 
                               "G2M" = "purple",
                               "G0" = "red",
                               "M" = "magenta",
                               "MG1" = "gold")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_blank()) +
  labs(title = "HCC4006 Osi DTPCs")


krt17_lito_cc_bp #Plot: 600 x 350



##### H1975 Osi DTPC Lito cell cycle by KRT17 status #####

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


#Create table with KRT17+ & KRT17- DTPCs by Lito cell cycle state
krt17_lito_cc <- as.data.frame(table(h1975_dtpc@meta.data$ExpressionLabel, 
                                     h1975_dtpc@meta.data$cell_cycle_expanded))

#Add column names
colnames(krt17_lito_cc) <- c("KRT17_status", "Lito_CC", "Number_of_cells")

krt17_lito_cc

#Calculate percentage of cells occupying each cell cycle state by KRT17 status
cell_counts_pct <- krt17_lito_cc %>%
  group_by(KRT17_status) %>%
  mutate(total_per_krt17 = sum(Number_of_cells),
         pct = (Number_of_cells / total_per_krt17) * 100)

#Round percentage to 2 decimals for plotting
cell_counts_pct$pct <- round(cell_counts_pct$pct, 2)

#Draw stacked barchart for KRT17+ vs KRT17- Osi DTPC Lito CC states
krt17_lito_cc_bp <- ggplot(cell_counts_pct, aes(x = Lito_CC, y = pct, fill = Lito_CC)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = pct), vjust = -1) +
  facet_wrap(~ KRT17_status) +
  scale_y_continuous("Percentage of cells", limits = c(0, 100)) +
  scale_fill_manual(values = c("G1S" = "lightgreen", 
                               "S" = "blue", 
                               "G2M" = "purple",
                               "G0" = "red",
                               "M" = "magenta",
                               "MG1" = "gold")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_blank()) +
  labs(title = "H1975 Osi DTPCs")


krt17_lito_cc_bp #Plot: 600 x 350









