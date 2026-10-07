library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(ggrepel)



##### HCC827 Osi DTPC vs DMSO DE genes volcano plots #####

#Read in data
hcc827_de <- HCC827_Osi_DTPC_v_DMSO_DE_marker_genes_volcano_plot

###ETV5, HOPX, SOX4 volcano plot

#Note: for genes with padj = 0, using -Minuslog10(padj) matching dataset limit

#Identify genes to label
genes_to_label <- c("ETV5", "HOPX", "SOX4")

#Set labels if in DE list
hcc827_de$label <- ifelse(hcc827_de$Gene %in% genes_to_label, hcc827_de$Gene, NA)

#Draw labeled volcano plot for ETV5, HOPX, and SOX4
hcc827_de_vp <- ggplot(hcc827_de, aes(x = avg_log2FC_Osi_DTPC_v_DMSO, y = Minus_log10_padj)) +
  geom_point(aes(color = avg_log2FC_Osi_DTPC_v_DMSO > 0)) +
  scale_color_manual(values = c("FALSE" = "blue",
                                "TRUE" = "red")) +
  geom_point(data = subset(hcc827_de, label != 0),
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
  scale_x_continuous("Average Log2FC") +
  scale_y_continuous("-Log10(padj)",
                     limits = c(0,307),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        legend.position = "none") +
  labs(title = "HCC827")

hcc827_de_vp #400 x 500 (WXH)


###KRT17 volcano plot

#Identify genes to label
genes_to_label <- c("KRT17")

#Set labels if in DE list
hcc827_de$label <- ifelse(hcc827_de$Gene %in% genes_to_label, hcc827_de$Gene, NA)

#Draw labeled volcano plot for ETV5, HOPX, and SOX4
hcc827_de_krt17_vp <- ggplot(hcc827_de, aes(x = avg_log2FC_Osi_DTPC_v_DMSO, y = Minus_log10_padj)) +
  geom_point(aes(color = avg_log2FC_Osi_DTPC_v_DMSO > 0)) +
  scale_color_manual(values = c("FALSE" = "blue",
                                "TRUE" = "red")) +
  geom_point(data = subset(hcc827_de, label != 0),
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
  scale_x_continuous("Average Log2FC") +
  scale_y_continuous("-Log10(padj)",
                     limits = c(0,307),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        legend.position = "none") +
  labs(title = "HCC827")

hcc827_de_krt17_vp #400 x 500 (WXH)



##### HCC4006 Osi DTPC vs DMSO DE genes volcano plots #####

#Read in data
hcc4006_de <- HCC4006_Osi_DTPC_v_DMSO_DE_marker_genes_volcano_plot

###ETV5, HOPX, SOX4 volcano plot

#Note: for genes with padj = 0, using -Minuslog10(padj) matching dataset limit

#Identify genes to label
genes_to_label <- c("ETV5", "HOPX", "SOX4")

#Set labels if in DE list
hcc4006_de$label <- ifelse(hcc4006_de$Gene %in% genes_to_label, hcc4006_de$Gene, NA)

#Draw labeled volcano plot for ETV5, HOPX, and SOX4
hcc4006_de_vp <- ggplot(hcc4006_de, aes(x = avg_log2FC_Osi_DTPC_v_DMSO, y = Minus_log10_padj)) +
  geom_point(aes(color = avg_log2FC_Osi_DTPC_v_DMSO > 0)) +
  scale_color_manual(values = c("FALSE" = "blue",
                                "TRUE" = "red")) +
  geom_point(data = subset(hcc4006_de, label != 0),
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
  scale_x_continuous("Average Log2FC") +
  scale_y_continuous("-Log10(padj)",
                     limits = c(0,307),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        legend.position = "none") +
  labs(title = "HCC4006")

hcc4006_de_vp #400 x 500 (WXH)


###KRT17 volcano plot

#Identify genes to label
genes_to_label <- c("KRT17")

#Set labels if in DE list
hcc4006_de$label <- ifelse(hcc4006_de$Gene %in% genes_to_label, hcc4006_de$Gene, NA)

#Draw labeled volcano plot for ETV5, HOPX, and SOX4
hcc4006_de_krt17_vp <- ggplot(hcc4006_de, aes(x = avg_log2FC_Osi_DTPC_v_DMSO, y = Minus_log10_padj)) +
  geom_point(aes(color = avg_log2FC_Osi_DTPC_v_DMSO > 0)) +
  scale_color_manual(values = c("FALSE" = "blue",
                                "TRUE" = "red")) +
  geom_point(data = subset(hcc4006_de, label != 0),
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
  scale_x_continuous("Average Log2FC") +
  scale_y_continuous("-Log10(padj)",
                     limits = c(0,307),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        legend.position = "none") +
  labs(title = "HCC4006")

hcc4006_de_krt17_vp #400 x 500 (WXH)




##### H1975 Osi DTPC vs DMSO DE genes volcano plots #####

#Read in data
h1975_de <- H1975_Osi_DTPC_v_DMSO_DE_marker_genes_volcano_plot

###ETV5, HOPX, SOX4 volcano plot

#Note: for genes with padj = 0, using -Minuslog10(padj) matching dataset limit

#Identify genes to label
genes_to_label <- c("ETV5", "HOPX", "SOX4")

#Set labels if in DE list
h1975_de$label <- ifelse(h1975_de$Gene %in% genes_to_label, h1975_de$Gene, NA)

#Draw labeled volcano plot for ETV5, HOPX, and SOX4
h1975_de_vp <- ggplot(h1975_de, aes(x = avg_log2FC_Osi_DTPC_v_DMSO, y = Minus_log10_padj)) +
  geom_point(aes(color = avg_log2FC_Osi_DTPC_v_DMSO > 0)) +
  scale_color_manual(values = c("FALSE" = "blue",
                                "TRUE" = "red")) +
  geom_point(data = subset(h1975_de, label != 0),
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
  scale_x_continuous("Average Log2FC") +
  scale_y_continuous("-Log10(padj)",
                     limits = c(0,307),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        legend.position = "none") +
  labs(title = "H1975")

h1975_de_vp #400 x 500 (WXH)


###KRT17 volcano plot

#Identify genes to label
genes_to_label <- c("KRT17")

#Set labels if in DE list
h1975_de$label <- ifelse(h1975_de$Gene %in% genes_to_label, h1975_de$Gene, NA)

#Draw labeled volcano plot for ETV5, HOPX, and SOX4
h1975_de_krt17_vp <- ggplot(h1975_de, aes(x = avg_log2FC_Osi_DTPC_v_DMSO, y = Minus_log10_padj)) +
  geom_point(aes(color = avg_log2FC_Osi_DTPC_v_DMSO > 0)) +
  scale_color_manual(values = c("FALSE" = "blue",
                                "TRUE" = "red")) +
  geom_point(data = subset(h1975_de, label != 0),
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
  scale_x_continuous("Average Log2FC") +
  scale_y_continuous("-Log10(padj)",
                     limits = c(0,307),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        legend.position = "none") +
  labs(title = "H1975")

h1975_de_krt17_vp #400 x 500 (WXH)






