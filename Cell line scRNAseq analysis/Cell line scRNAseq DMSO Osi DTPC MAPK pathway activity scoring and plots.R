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



##### HCC827 DMSO and HCC827 Osi DTPC integrated UCell scoring analysis #####

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


###Draw MAPK Pathway Activity Score Faceted UMAP (DMSO vs. Osi DTPC)

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc827_combined, "umap"))
umap_data$mapk_act <- hcc827_combined@meta.data$signature_1MAPK_act
umap_data$orig.ident <- hcc827_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
mapk_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = mapk_act)) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "MAPK \nPathway\nActivity \nScore") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

mapk_umap_f #Plot: 600 x 300

#Save plot 
ggsave(mapk_umap_f, file = "2025.6.17 HCC827 DMSO Osi DTPC MAPK Pathway activity score ucell faceted umap.svg", width = 6, height = 3)


###MAPK Pathway Activity Score violin plot
hcc827_mapk_vp <- VlnPlot(hcc827_combined, features = "signature_1MAPK_act", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "MAPK Pathway \nActivity Score")

hcc827_mapk_vp #Plot: 350 x 500


##### HCC4006 DMSO and HCC4006 Osi DTPC integrated UCell scoring analysis #####

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


###Draw MAPK Pathway Activity Score Faceted UMAP (DMSO vs. Osi DTPC)

#Create dataframe
umap_data <- as.data.frame(Embeddings(hcc4006_combined, "umap"))
umap_data$mapk_act <- hcc4006_combined@meta.data$signature_1MAPK_act
umap_data$orig.ident <- hcc4006_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
mapk_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = mapk_act)) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "MAPK \nPathway\nActivity \nScore") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

mapk_umap_f #Plot: 600 x 300

#Save plot
ggsave(mapk_umap_f, file = "2025.6.17 HCC4006 DMSO Osi DTPC MAPK Pathway activity score ucell faceted umap.svg", width = 6, height = 3)



###MAPK Pathway Activity Score violin plot
hcc4006_mapk_vp <- VlnPlot(hcc4006_combined, features = "signature_1MAPK_act", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "MAPK Pathway \nActivity Score")

hcc4006_mapk_vp #Plot: 350 x 500



##### H1975 DMSO and H1975 Osi DTPC integrated UCell scoring analysis #####

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


###Draw MAPK Pathway Activity Score Faceted UMAP (DMSO vs. Osi DTPC)

#Create dataframe
umap_data <- as.data.frame(Embeddings(h1975_combined, "umap"))
umap_data$mapk_act <- h1975_combined@meta.data$signature_1MAPK_act
umap_data$orig.ident <- h1975_combined@meta.data$orig.ident

#Draw faceted UMAP by orig.ident (DMSO vs DTPC)
mapk_umap_f <- ggplot(umap_data, aes(x = umap_1, y = umap_2, colour = mapk_act)) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~orig.ident) +
  scale_x_continuous("UMAP1") +
  scale_y_continuous("UMAP2") +
  labs(color = "MAPK \nPathway\nActivity \nScore") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"))

mapk_umap_f #Plot: 600 x 300

#Save plot
ggsave(mapk_umap_f, file = "2025.6.17 H1975 DMSO Osi DTPC MAPK Pathway activity score ucell faceted umap.svg", width = 6, height = 3)


###MAPK Pathway Activity Score violin plot 
h1975_mapk_vp <- VlnPlot(h1975_combined, features = "signature_1MAPK_act", group.by = "orig.ident", pt.size = FALSE) +
  stat_compare_means(label.x = 1.3) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Activity") +
  labs(title = "MAPK Pathway \nActivity Score")

h1975_mapk_vp #Plot: 350 x 500

