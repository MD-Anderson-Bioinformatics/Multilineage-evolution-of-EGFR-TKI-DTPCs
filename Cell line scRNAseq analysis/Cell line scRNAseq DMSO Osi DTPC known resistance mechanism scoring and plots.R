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


##### HCC827 DMSO and HCC827 Osi DTPC integrated UCell scoring analysis #####

#Load in data
hcc827_combined <- LoadSeuratRds("HCC827_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")

###EMT scoring

#Read in EMT score
emt_genes <- Hallmark_EMT_genes

#Define emt features
emt_genes_features <- list(c(emt_genes$Gene))

#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=emt_genes_features, name="Hallmark_EMT")


###Draw Hallmark EMT Faceted UMAP (DMSO vs. Osi DTPC)

#Create dataframe
emt_umap_data <- as.data.frame(Embeddings(hcc827_combined, "umap"))
emt_umap_data$emt <- hcc827_combined@meta.data$signature_1Hallmark_EMT
emt_umap_data$orig.ident <- hcc827_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
emt_umap_f <- ggplot(emt_umap_data, aes(x = umap_1, y = umap_2, colour = emt)) +
  geom_point() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "Hallmark \nEMT") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

emt_umap_f #Plot: 600 x 300

#Save plot
ggsave(emt_umap_f, file = "2025.6.17 HCC827 DMSO Osi DTPC Hallmark EMT ucell faceted umap.svg", width = 6, height = 3)


###Hallmark EMT violin plot 
hcc827_emt_vp <- VlnPlot(hcc827_combined, features = "signature_1Hallmark_EMT", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "Hallmark EMT")

hcc827_emt_vp #Plot: 350 x 500



###IL6 pathway

#Read in IL6 signaling genes
il6_genes <- Hallmark_IL6_JAK_STAT_signaling_genes

#Set IL6 signaling features
il6_genes_features <- list(c(il6_genes$Gene))

#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=il6_genes_features, name="il6")


###IL6 signaling Faceted UMAP

#Create dataframe
il6_sig_umap_data <- as.data.frame(Embeddings(hcc827_combined, "umap"))
il6_sig_umap_data$il6_sig <- hcc827_combined@meta.data$signature_1il6
il6_sig_umap_data$orig.ident <- hcc827_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
il6_sig_umap_f <- ggplot(il6_sig_umap_data, aes(x = umap_1, y = umap_2, colour = il6_sig)) +
  geom_point() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "IL6 \nSignaling") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

il6_sig_umap_f #Plot: 600 x 300

#Save plot
ggsave(il6_sig_umap_f, file = "2025.6.17 HCC827 DMSO Osi DTPC Hallmark IL6 signaling ucell faceted umap.svg", width = 6, height = 3)


###IL6 signaling violin plot (6/17/2025 update; no points)
hcc827_il6_vp <- VlnPlot(hcc827_combined, features = "signature_1il6", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "IL6 Signaling")

hcc827_il6_vp #350 x 500



###RTK signaling scoring

#Read in Reactome RTK signaling genes
rtk_sig_genes <- Reactome_signaling_by_RTKs_genes

#Set Reactome RTK signaling features
rtk_sig_genes_features <- list(c(rtk_sig_genes$Gene))

#Score cells using UCell
hcc827_combined <- AddModuleScore_UCell(hcc827_combined, 
                                        features=rtk_sig_genes_features, name="RTK_sig")


###RTK signaling Faceted UMAP

#Create dataframe
rtk_sig_umap_data <- as.data.frame(Embeddings(hcc827_combined, "umap"))
rtk_sig_umap_data$rtk_sig <- hcc827_combined@meta.data$signature_1RTK_sig
rtk_sig_umap_data$orig.ident <- hcc827_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
rtk_sig_umap_f <- ggplot(rtk_sig_umap_data, aes(x = umap_1, y = umap_2, colour = rtk_sig)) +
  geom_point() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "RTK \nSignaling") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

rtk_sig_umap_f #Plot: 600 x 300

#Save plot
ggsave(rtk_sig_umap_f, file = "2025.6.17 HCC827 DMSO Osi DTPC Reactome RTK signaling ucell faceted umap.svg", width = 6, height = 3)


###RTK signaling signaling violin plot
hcc827_rtk_sig_vp <- VlnPlot(hcc827_combined, features = "signature_1RTK_sig", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "RTK Signaling")

hcc827_rtk_sig_vp #350 x 500




##### HCC4006 DMSO and HCC4006 Osi DTPC integrated UCell scoring analysis #####

#Load in data
hcc4006_combined <- LoadSeuratRds("HCC4006_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###EMT scoring

#Read in EMT score
emt_genes <- Hallmark_EMT_genes

#Define emt features
emt_genes_features <- list(c(emt_genes$Gene))

#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=emt_genes_features, name="Hallmark_EMT")


###Draw Hallmark EMT Faceted UMAP (DMSO vs. Osi DTPC)

#Create dataframe
emt_umap_data <- as.data.frame(Embeddings(hcc4006_combined, "umap"))
emt_umap_data$emt <- hcc4006_combined@meta.data$signature_1Hallmark_EMT
emt_umap_data$orig.ident <- hcc4006_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
emt_umap_f <- ggplot(emt_umap_data, aes(x = umap_1, y = umap_2, colour = emt)) +
  geom_point() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "Hallmark \nEMT") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

emt_umap_f #Plot: 600 x 300

#Save plot
ggsave(emt_umap_f, file = "2025.6.17 HCC4006 DMSO Osi DTPC Hallmark EMT ucell faceted umap.svg", width = 6, height = 3)


###Hallmark EMT violin plot
hcc4006_emt_vp <- VlnPlot(hcc4006_combined, features = "signature_1Hallmark_EMT", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "Hallmark EMT")

hcc4006_emt_vp #Plot: 350 x 500



###IL6 pathway

#Read in IL6 signaling genes
il6_genes <- Hallmark_IL6_JAK_STAT_signaling_genes

#Set IL6 signaling features
il6_genes_features <- list(c(il6_genes$Gene))

#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=il6_genes_features, name="il6")


###IL6 signaling Faceted UMAP

#Create dataframe
il6_sig_umap_data <- as.data.frame(Embeddings(hcc4006_combined, "umap"))
il6_sig_umap_data$il6_sig <- hcc4006_combined@meta.data$signature_1il6
il6_sig_umap_data$orig.ident <- hcc4006_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
il6_sig_umap_f <- ggplot(il6_sig_umap_data, aes(x = umap_1, y = umap_2, colour = il6_sig)) +
  geom_point() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "IL6 \nSignaling") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

il6_sig_umap_f #Plot: 600 x 300

#Save plot
ggsave(il6_sig_umap_f, file = "2025.6.17 HCC4006 DMSO Osi DTPC Hallmark IL6 signaling ucell faceted umap.svg", width = 6, height = 3)


###IL6 signaling violin plot
hcc4006_il6_vp <- VlnPlot(hcc4006_combined, features = "signature_1il6", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "IL6 Signaling")

hcc4006_il6_vp #350 x 500



###RTK signaling scoring

#Read in Reactome RTK signaling genes
rtk_sig_genes <- Reactome_signaling_by_RTKs_genes

#Set Reactome RTK signaling features
rtk_sig_genes_features <- list(c(rtk_sig_genes$Gene))

#Score cells using UCell
hcc4006_combined <- AddModuleScore_UCell(hcc4006_combined, 
                                         features=rtk_sig_genes_features, name="RTK_sig")


###RTK signaling Faceted UMAP

#Create dataframe
rtk_sig_umap_data <- as.data.frame(Embeddings(hcc4006_combined, "umap"))
rtk_sig_umap_data$rtk_sig <- hcc4006_combined@meta.data$signature_1RTK_sig
rtk_sig_umap_data$orig.ident <- hcc4006_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
rtk_sig_umap_f <- ggplot(rtk_sig_umap_data, aes(x = umap_1, y = umap_2, colour = rtk_sig)) +
  geom_point() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "RTK \nSignaling") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

rtk_sig_umap_f #Plot: 600 x 300

#Save plot
ggsave(rtk_sig_umap_f, file = "2025.6.17 HCC4006 DMSO Osi DTPC Reactome RTK signaling ucell faceted umap.svg", width = 6, height = 3)


###RTK signaling signaling violin plot
hcc4006_rtk_sig_vp <- VlnPlot(hcc4006_combined, features = "signature_1RTK_sig", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "RTK Signaling")

hcc4006_rtk_sig_vp #350 x 500



##### H1975 DMSO and H1975 Osi DTPC integrated UCell scoring analysis #####

#Load in data
h1975_combined <- LoadSeuratRds("H1975_DMSO_Osi_DTPC_scTransformed_v2_singlet_CCAintegrated_LitoCC_annotated.Rds")


###EMT scoring

#Read in EMT score
emt_genes <- Hallmark_EMT_genes

#Define emt features
emt_genes_features <- list(c(emt_genes$Gene))

#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=emt_genes_features, name="Hallmark_EMT")


###Draw Hallmark EMT Faceted UMAP (DMSO vs. Osi DTPC)

#Create dataframe
emt_umap_data <- as.data.frame(Embeddings(h1975_combined, "umap"))
emt_umap_data$emt <- h1975_combined@meta.data$signature_1Hallmark_EMT
emt_umap_data$orig.ident <- h1975_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
emt_umap_f <- ggplot(emt_umap_data, aes(x = umap_1, y = umap_2, colour = emt)) +
  geom_point() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "Hallmark \nEMT") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

emt_umap_f #Plot: 600 x 300


#Save plot
ggsave(emt_umap_f, file = "2025.6.17 H1975 DMSO Osi DTPC Hallmark EMT ucell faceted umap.svg", width = 6, height = 3)


###Hallmark EMT violin plot
h1975_emt_vp <- VlnPlot(h1975_combined, features = "signature_1Hallmark_EMT", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "Hallmark EMT")

h1975_emt_vp #Plot: 350 x 500



###IL6 pathway

#Read in IL6 signaling genes
il6_genes <- Hallmark_IL6_JAK_STAT_signaling_genes

#Set IL6 signaling features
il6_genes_features <- list(c(il6_genes$Gene))

#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=il6_genes_features, name="il6")


###IL6 signaling Faceted UMAP

#Create dataframe
il6_sig_umap_data <- as.data.frame(Embeddings(h1975_combined, "umap"))
il6_sig_umap_data$il6_sig <- h1975_combined@meta.data$signature_1il6
il6_sig_umap_data$orig.ident <- h1975_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
il6_sig_umap_f <- ggplot(il6_sig_umap_data, aes(x = umap_1, y = umap_2, colour = il6_sig)) +
  geom_point() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "IL6 \nSignaling") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

il6_sig_umap_f #Plot: 600 x 300

#Save plot
ggsave(il6_sig_umap_f, file = "2025.6.17 H1975 DMSO Osi DTPC Hallmark IL6 signaling ucell faceted umap.svg", width = 6, height = 3)


###IL6 signaling violin plot
h1975_il6_vp <- VlnPlot(h1975_combined, features = "signature_1il6", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "IL6 Signaling")

h1975_il6_vp #350 x 500



###RTK signaling scoring

#Read in Reactome RTK signaling genes
rtk_sig_genes <- Reactome_signaling_by_RTKs_genes

#Set Reactome RTK signaling features
rtk_sig_genes_features <- list(c(rtk_sig_genes$Gene))

#Score cells using UCell
h1975_combined <- AddModuleScore_UCell(h1975_combined, 
                                       features=rtk_sig_genes_features, name="RTK_sig")



###RTK signaling Faceted UMAP

#Create dataframe
rtk_sig_umap_data <- as.data.frame(Embeddings(h1975_combined, "umap"))
rtk_sig_umap_data$rtk_sig <- h1975_combined@meta.data$signature_1RTK_sig
rtk_sig_umap_data$orig.ident <- h1975_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
rtk_sig_umap_f <- ggplot(rtk_sig_umap_data, aes(x = umap_1, y = umap_2, colour = rtk_sig)) +
  geom_point() +
  scale_color_viridis_c(option = "turbo") +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "RTK \nSignaling") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

rtk_sig_umap_f #Plot: 600 x 300

#Save plot
ggsave(rtk_sig_umap_f, file = "2025.6.17 H1975 DMSO Osi DTPC Reactome RTK signaling ucell faceted umap.svg", width = 6, height = 3)


###RTK signaling signaling violin plot (6/17/2025 update; no points)
h1975_rtk_sig_vp <- VlnPlot(h1975_combined, features = "signature_1RTK_sig", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "RTK Signaling")

h1975_rtk_sig_vp #350 x 500



