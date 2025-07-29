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
library(ggpubr)
library(devtools)
library(forcats)
library(scCustomize)
library(scico)


##### HCC827 DMSO and HCC827 Osi DTPC integrated UCell scoring analysis #####

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###AT2 gene scoring

#Read in AEC2 genes score
aec2_genes <- AEC2_genes

#Set AT2 geen features
aec2_genes_features <- list(c(aec2_genes$Gene))

#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=aec2_genes_features, name="aec2")


###AT2 Faceted UMAP

#Create dataframe
at2_umap_data <- as.data.frame(Embeddings(hcc827_combined, "umap"))
at2_umap_data$at2 <- hcc827_combined@meta.data$signature_1aec2
at2_umap_data$orig.ident <- hcc827_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
at2_umap_f <- ggplot(at2_umap_data, aes(x = umap_1, y = umap_2, colour = at2)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "AT2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

at2_umap_f #Plot: 600 x 300


#Save plot
ggsave(at2_umap_f, file = "2025.6.17 HCC827 DMSO Osi DTPC AEC2 ucell faceted umap.svg", width = 6, height = 3)


###AT2 violin plot
hcc827_at2_vp <- VlnPlot(hcc827_combined, features = "signature_1aec2", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "AEC2")

hcc827_at2_vp #Plot: 350 x 500


###ETV5 feature plot
hcc827_etv5_fp <- FeaturePlot(hcc827_combined, features = "ETV5", split.by = "orig.ident",
                              pt.size = 1, max.cutoff = 4) +
  theme(legend.position = "right")

hcc827_etv5_fp_1 <- hcc827_etv5_fp & scale_color_viridis_c(option = "magma", name = "ETV5") & 
  scale_y_continuous("UMAP2") & 
  scale_x_continuous("UMAP1") & 
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

#Save plot
ggsave(hcc827_etv5_fp_1, file = "2025.6.17 HCC827 DMSO Osi DTPC ETV5 faceted feature plot.svg", width = 6.5, height = 3)



###AT1 gene scoring

#Read in AT1 genes score
aec1_genes <- AEC1_genes

#Set AT1 gene features
aec1_genes_features <- list(c(aec1_genes$Gene))

#Score cells using Ucell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=aec1_genes_features, name="aec1")


###AT1 Faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc827_combined, "umap"))
umap_data$at1 <- hcc827_combined@meta.data$signature_1aec1
umap_data$orig.ident <- hcc827_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
at1_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = at1)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "AT1") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

at1_umap_f

#Save plot
ggsave(at1_umap_f, file = "2025.6.17 HCC827 DMSO Osi DTPC AEC1 ucell faceted umap.svg", width = 6, height = 3)


###AT1 violin plot
hcc827_at1_vp <- VlnPlot(hcc827_combined, features = "signature_1aec1", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "AEC1")

hcc827_at1_vp #Plot: 350 x 500


###HOPX feature plot
hcc827_hopx_fp <- FeaturePlot(hcc827_combined, features = "HOPX", split.by = "orig.ident",
                              pt.size = 1, max.cutoff = 4)  +
  theme(legend.position = "right")

hcc827_hopx_fp_1 <- hcc827_hopx_fp & scale_color_viridis_c(option = "magma", name = "HOPX") & 
  scale_y_continuous("UMAP2") & 
  scale_x_continuous("UMAP1") & 
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

#Save plot
ggsave(hcc827_hopx_fp_1, file = "2025.6.17 HCC827 DMSO Osi DTPC HOPX faceted feature plot.svg", width = 6.5, height = 3)



###SOX4 feature plot
hcc827_sox4_fp <- FeaturePlot(hcc827_combined, features = "SOX4", split.by = "orig.ident",
                              pt.size = 1, max.cutoff = 4)  +
  theme(legend.position = "right")

hcc827_sox4_fp_1 <- hcc827_sox4_fp & scale_color_viridis_c(option = "magma", name = "SOX4") & 
  scale_y_continuous("UMAP2") & 
  scale_x_continuous("UMAP1") & 
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

#Save plot
ggsave(hcc827_sox4_fp_1, file = "2025.6.17 HCC827 DMSO Osi DTPC SOX4 faceted feature plot.svg", width = 6.5, height = 3)



###AT0 gene scoring

#Read in AT0 genes score
at0_genes <- AT0_genes

#Set AT0 gene features
at0_genes_features <- list(c(at0_genes$Gene))

#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=at0_genes_features, name="at0")


###AT0 Faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc827_combined, "umap"))
umap_data$at0 <- hcc827_combined@meta.data$signature_1at0
umap_data$orig.ident <- hcc827_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
at0_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = at0)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "AT0") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

at0_umap_f

#Save plot
ggsave(at0_umap_f, file = "2025.6.17 HCC827 DMSO Osi DTPC AT0 ucell faceted umap.svg", width = 6, height = 3)


###AT0 violin plot
hcc827_at0_vp <- VlnPlot(hcc827_combined, features = "signature_1at0", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "AT0")

hcc827_at0_vp #Plot: 350 x 500



###Aberrant basaloid scoring (Wang)

#Read in aberrant basaloid genes
wang_ab_basaloid_genes <- Wang_aberrant_basaloid_signature

#Set aberrant basaloid features
wang_ab_basaloid_genes_features <- list(c(wang_ab_basaloid_genes$Gene))

#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=wang_ab_basaloid_genes_features, name="ab_basaloid_wang")


###Aberrant basaloid Faceted UMAP (Wang)

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc827_combined, "umap"))
umap_data$ab_basaloid_wang <- hcc827_combined@meta.data$signature_1ab_basaloid_wang
umap_data$orig.ident <- hcc827_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
wang_ab_basaloid_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = ab_basaloid_wang)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "Aberrant \nBasaloid") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

wang_ab_basaloid_umap_f


#Save plot
ggsave(wang_ab_basaloid_umap_f, file = "2025.6.18 HCC827 DMSO Osi DTPC Aberrant basaloid ucell faceted umap.svg", width = 6, height = 3)


###Aberrant basaloid violin plot
hcc827_ab_basaloid_vp <- VlnPlot(hcc827_combined, features = "signature_1ab_basaloid_wang", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "Aberrant Basaloid")

hcc827_ab_basaloid_vp #Plot: 350 x 500




##### HCC4006 DMSO and HCC4006 Osi DTPC integrated UCell scoring analysis #####

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###AT2 gene scoring

#Read in AEC2 genes score
aec2_genes <- AEC2_genes

#Set AT2 geen features
aec2_genes_features <- list(c(aec2_genes$Gene))

#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=aec2_genes_features, name="aec2")


###AT2 Faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc4006_combined, "umap"))
umap_data$at2 <- hcc4006_combined@meta.data$signature_1aec2
umap_data$orig.ident <- hcc4006_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
at2_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = at2)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "AT2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

at2_umap_f #Plot: 600 x 300

#Save plot
ggsave(at2_umap_f, file = "2025.6.17 HCC4006 DMSO Osi DTPC AEC2 ucell faceted umap.svg", width = 6, height = 3)


###AT2 violin plot
hcc4006_at2_vp <- VlnPlot(hcc4006_combined, features = "signature_1aec2", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "AEC2")

hcc4006_at2_vp #Plot: 350 x 500


###ETV5 feature plot
hcc4006_etv5_fp <- FeaturePlot(hcc4006_combined, features = "ETV5", split.by = "orig.ident",
                               pt.size = 1, max.cutoff = 4) +
  theme(legend.position = "right")

hcc4006_etv5_fp_1 <- hcc4006_etv5_fp & scale_color_viridis_c(option = "magma", name = "ETV5") & 
  scale_y_continuous("UMAP2") & 
  scale_x_continuous("UMAP1") & 
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

#Save plot
ggsave(hcc4006_etv5_fp_1, file = "2025.6.17 HCC4006 DMSO Osi DTPC ETV5 faceted feature plot.svg", width = 6.5, height = 3)



###AT1 gene scoring

#Read in AT1 genes score
aec1_genes <- AEC1_genes

#Set AT1 gene features
aec1_genes_features <- list(c(aec1_genes$Gene))

#Score cells using Ucell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=aec1_genes_features, name="aec1")


###AT1 Faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc4006_combined, "umap"))
umap_data$at1 <- hcc4006_combined@meta.data$signature_1aec1
umap_data$orig.ident <- hcc4006_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
at1_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = at1)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "AT1") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

at1_umap_f

#Save plot
ggsave(at1_umap_f, file = "2025.6.17 HCC4006 DMSO Osi DTPC AEC1 ucell faceted umap.svg", width = 6, height = 3)


###AT1 violin plot
hcc4006_at1_vp <- VlnPlot(hcc4006_combined, features = "signature_1aec1", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "AEC1")

hcc4006_at1_vp #Plot: 350 x 500


###HOPX feature plot
hcc4006_hopx_fp <- FeaturePlot(hcc4006_combined, features = "HOPX", split.by = "orig.ident",
                               pt.size = 1, max.cutoff = 4) +
  theme(legend.position = "right")

hcc4006_hopx_fp_1 <- hcc4006_hopx_fp & scale_color_viridis_c(option = "magma", name = "HOPX") & 
  scale_y_continuous("UMAP2") & 
  scale_x_continuous("UMAP1") & 
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

#Save plot
ggsave(hcc4006_hopx_fp_1, file = "2025.6.17 HCC4006 DMSO Osi DTPC HOPX faceted feature plot.svg", width = 6.5, height = 3)


###SOX4 feature plot
hcc4006_sox4_fp <- FeaturePlot(hcc4006_combined, features = "SOX4", split.by = "orig.ident",
                               pt.size = 1, max.cutoff = 4) +
  theme(legend.position = "right")

hcc4006_sox4_fp_1 <- hcc4006_sox4_fp & scale_color_viridis_c(option = "magma", name = "SOX4") & 
  scale_y_continuous("UMAP2") & 
  scale_x_continuous("UMAP1") & 
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

#Save plot
ggsave(hcc4006_sox4_fp_1, file = "2025.6.17 HCC4006 DMSO Osi DTPC SOX4 faceted feature plot.svg", width = 6.5, height = 3)



###AT0 gene scoring

#Read in AT0 genes score
at0_genes <- AT0_genes

#Set AT0 gene features
at0_genes_features <- list(c(at0_genes$Gene))

#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=at0_genes_features, name="at0")


###AT0 Faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc4006_combined, "umap"))
umap_data$at0 <- hcc4006_combined@meta.data$signature_1at0
umap_data$orig.ident <- hcc4006_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
at0_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = at0)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "AT0") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

at0_umap_f

#Save plot
ggsave(at0_umap_f, file = "2025.6.17 HCC4006 DMSO Osi DTPC AT0 ucell faceted umap.svg", width = 6, height = 3)


###AT0 violin plot
hcc4006_at0_vp <- VlnPlot(hcc4006_combined, features = "signature_1at0", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "AT0")

hcc4006_at0_vp #Plot: 350 x 500


###Aberrant basaloid scoring (Wang)

#Read in aberrant basaloid genes
wang_ab_basaloid_genes <- Wang_aberrant_basaloid_signature

#Set aberrant basaloid features
wang_ab_basaloid_genes_features <- list(c(wang_ab_basaloid_genes$Gene))

#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=wang_ab_basaloid_genes_features, name="ab_basaloid_wang")


###Aberrant basaloid Faceted UMAP (Wang)

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc4006_combined, "umap"))
umap_data$ab_basaloid_wang <- hcc4006_combined@meta.data$signature_1ab_basaloid_wang
umap_data$orig.ident <- hcc4006_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
wang_ab_basaloid_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = ab_basaloid_wang)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "Aberrant \nBasaloid") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

wang_ab_basaloid_umap_f

#Save plot
ggsave(wang_ab_basaloid_umap_f, file = "2025.6.18 HCC4006 DMSO Osi DTPC Aberrant basaloid ucell faceted umap.svg", width = 6, height = 3)


###Aberrant basaloid violin plot
hcc4006_ab_basaloid_vp <- VlnPlot(hcc4006_combined, features = "signature_1ab_basaloid_wang", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "Aberrant Basaloid")

hcc4006_ab_basaloid_vp #Plot: 350 x 500



##### H1975 DMSO and H1975 Osi DTPC integrated UCell scoring analysis #####

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

###AT2 gene scoring

#Read in AEC2 genes score
aec2_genes <- AEC2_genes

#Set AT2 geen features
aec2_genes_features <- list(c(aec2_genes$Gene))

#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=aec2_genes_features, name="aec2")


###AT2 Faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(h1975_combined, "umap"))
umap_data$at2 <- h1975_combined@meta.data$signature_1aec2
umap_data$orig.ident <- h1975_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
at2_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = at2)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "AT2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

at2_umap_f #Plot: 600 x 300

#Save plot
ggsave(at2_umap_f, file = "2025.6.17 H1975 DMSO Osi DTPC AEC2 ucell faceted umap.svg", width = 6, height = 3)


###AT2 violin plot
h1975_at2_vp <- VlnPlot(h1975_combined, features = "signature_1aec2", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "AEC2")

h1975_at2_vp #Plot: 350 x 500


###ETV5 feature plot
h1975_etv5_fp <- FeaturePlot(h1975_combined, features = "ETV5", split.by = "orig.ident",
                             pt.size = 1, max.cutoff = 4) +
  theme(legend.position = "right")

h1975_etv5_fp_1 <- h1975_etv5_fp & scale_color_viridis_c(option = "magma", name = "ETV5") & 
  scale_y_continuous("UMAP2") & 
  scale_x_continuous("UMAP1") & 
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))


#Save plot
ggsave(h1975_etv5_fp_1, file = "2025.6.17 H1975 DMSO Osi DTPC ETV5 faceted feature plot.svg", width = 6.5, height = 3)



###AT1 gene scoring

#Read in AT1 genes score
aec1_genes <- AEC1_genes

#Set AT1 gene features
aec1_genes_features <- list(c(aec1_genes$Gene))

#Score cells using Ucell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=aec1_genes_features, name="aec1")


###AT1 Faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(h1975_combined, "umap"))
umap_data$at1 <- h1975_combined@meta.data$signature_1aec1
umap_data$orig.ident <- h1975_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
at1_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = at1)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "AT1") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

at1_umap_f

#Save plot
ggsave(at1_umap_f, file = "2025.6.17 H1975 DMSO Osi DTPC AEC1 ucell faceted umap.svg", width = 6, height = 3)


###AT1 violin plot
h1975_at1_vp <- VlnPlot(h1975_combined, features = "signature_1aec1", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "AEC1")

h1975_at1_vp #Plot: 350 x 500


###HOPX feature plot
h1975_hopx_fp <- FeaturePlot(h1975_combined, features = "HOPX", split.by = "orig.ident",
                             pt.size = 1, max.cutoff = 4) +
  theme(legend.position = "right")

h1975_hopx_fp_1 <- h1975_hopx_fp & scale_color_viridis_c(option = "magma", name = "HOPX") & 
  scale_y_continuous("UMAP2") & 
  scale_x_continuous("UMAP1") & 
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))


#Save plot
ggsave(h1975_hopx_fp_1, file = "2025.6.17 H1975 DMSO Osi DTPC HOPX faceted feature plot.svg", width = 6.5, height = 3)



###SOX4 feature plot
h1975_sox4_fp <- FeaturePlot(h1975_combined, features = "SOX4", split.by = "orig.ident",
                             pt.size = 1, max.cutoff = 4) +
  theme(legend.position = "right")

h1975_sox4_fp_1 <- h1975_sox4_fp & scale_color_viridis_c(option = "magma", name = "SOX4") & 
  scale_y_continuous("UMAP2") & 
  scale_x_continuous("UMAP1") & 
  theme(axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

#Save plot
ggsave(h1975_sox4_fp_1, file = "2025.6.17 H1975 DMSO Osi DTPC SOX4 faceted feature plot.svg", width = 6.5, height = 3)


###AT0 gene scoring

#Read in AT0 genes score
at0_genes <- AT0_genes

#Set AT0 gene features
at0_genes_features <- list(c(at0_genes$Gene))

#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=at0_genes_features, name="at0")


###AT0 Faceted UMAP

#Create dataframe
umap_data <- as.data.frame(Embeddings(h1975_combined, "umap"))
umap_data$at0 <- h1975_combined@meta.data$signature_1at0
umap_data$orig.ident <- h1975_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
at0_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = at0)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "AT0") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

at0_umap_f

#Save plot
ggsave(at0_umap_f, file = "2025.6.17 H1975 DMSO Osi DTPC AT0 ucell faceted umap.svg", width = 6, height = 3)


###AT0 violin plot
h1975_at0_vp <- VlnPlot(h1975_combined, features = "signature_1at0", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "AT0")

h1975_at0_vp #Plot: 350 x 500



###Aberrant basaloid scoring (Wang)

#Read in aberrant basaloid genes
wang_ab_basaloid_genes <- Wang_aberrant_basaloid_signature

#Set aberrant basaloid features
wang_ab_basaloid_genes_features <- list(c(wang_ab_basaloid_genes$Gene))

#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=wang_ab_basaloid_genes_features, name="ab_basaloid_wang")


###Aberrant basaloid Faceted UMAP (Wang)

#Create dataframe
umap_data <- as.data.frame(Embeddings(h1975_combined, "umap"))
umap_data$ab_basaloid_wang <- h1975_combined@meta.data$signature_1ab_basaloid_wang
umap_data$orig.ident <- h1975_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
wang_ab_basaloid_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = ab_basaloid_wang)) +
  geom_point() +
  scale_color_scico(palette = "vik") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "Aberrant \nBasaloid") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

wang_ab_basaloid_umap_f

#Save plot
ggsave(wang_ab_basaloid_umap_f, file = "2025.6.18 H1975 DMSO Osi DTPC Aberrant basaloid ucell faceted umap.svg", width = 6, height = 3)


###Aberrant basaloid violin plot
h1975_ab_basaloid_vp <- VlnPlot(h1975_combined, features = "signature_1ab_basaloid_wang", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "Aberrant Basaloid")

h1975_ab_basaloid_vp #Plot: 350 x 500






