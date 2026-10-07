library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)


##### HCC827 HCC4006 H1975 Osi DTPC vs DMSO Bulk RNAseq ORA plot #####

#Read in data
bulk_ora <- HCC827_HCC4006_H1975_Osi_DTPC_vs_DMSO_bulk_RNAseq_ORA_combined_results

bulk_ora$color <- ifelse(bulk_ora$enrichmentRatio > 0, "up", "down")

#Draw ORA barchart for consensus upregulated & downregulated genes
bulk_ora_p <- ggplot(bulk_ora, aes(x = reorder(bulk_ora$description, bulk_ora$enrichmentRatio), y = bulk_ora$enrichmentRatio, fill = bulk_ora$color)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  scale_y_continuous("Enrichment Ratio") +
  scale_fill_manual(values = c("up" = "red",
                               "down" = "blue")) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        legend.position = "none")

bulk_ora_p #1000 x 1700









