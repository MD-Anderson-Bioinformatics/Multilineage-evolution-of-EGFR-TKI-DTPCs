library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(ggpubr)


##### EGFR TKI DTPC shared upregulated & downregulated HOMER known motif analysis plots #####

###Upregulated genes

#Read in data
osi_dtpc_u_km <- HCC827_HCC4006_H1975_DTPC_shared_upregulated_HOMER_known_motif_plot

#Plot sigmoid plot of shared upregulated motif rank vs. -log10(motif pvalue)
osi_dtpc_u_km_p <- ggplot(osi_dtpc_km, aes(x = osi_dtpc_km$Rank, y = osi_dtpc_km$`q-value (Benjamini)`)) +
  geom_point(alpha = 1, aes(color = ifelse(Label != 0, 'red', 'black'))) +
  geom_label_repel(aes(label = osi_dtpc_km$Label),
                   box.padding = .1,
                   ylim = c(-0.05, -0.01)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.position = "none") +
  scale_x_reverse("Rank", breaks = c(30, 20, 10, 1)) +
  scale_y_reverse("BH q-value", limits = c(0.05, 0)) 

osi_dtpc_u_km_p


###Downregulated genes

#Read in data
osi_dtpc_d_km <- HCC827_HCC4006_H1975_DTPC_shared_downregulated_HOMER_known_motif_plot

#Plot sigmoid plot of shared upregulated motif rank vs. -log10(motif pvalue)
osi_dtpc_d_km_p <- ggplot(osi_dtpc_d_km, aes(x = osi_dtpc_d_km$Rank, y = osi_dtpc_d_km$`q-value (Benjamini)`)) +
  geom_point(alpha = 1, aes(color = ifelse(Label != 0, 'red', 'black'))) +
  geom_label_repel(aes(label = osi_dtpc_d_km$Label),
                   box.padding = .1,
                   ylim = c(-0.05, -0.01)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.position = "none") +
  scale_x_reverse("Rank", breaks = c(30, 20, 10, 1)) +
  scale_y_reverse("BH q-value", limits = c(0.05, 0)) 

osi_dtpc_d_km_p





