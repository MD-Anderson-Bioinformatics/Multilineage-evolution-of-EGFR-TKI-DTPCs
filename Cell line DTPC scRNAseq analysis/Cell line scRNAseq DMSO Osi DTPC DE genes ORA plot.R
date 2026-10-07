library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)



##### HCC827 HCC4006 H1975 Osi DTPC vs DMSO scRNAseq ORA plot #####


#Read in data
sc_ora <- HCC827_HCC4006_H1975_Osi_DTPC_vs_DMSO_scRNAseq_ORA_combined_results

sc_ora$color <- ifelse(sc_ora$enrichmentRatio > 0, "up", "down")

#Create a unique index for each row to use on the x-axis
sc_ora <- sc_ora %>%
  arrange(enrichmentRatio) %>% 
  mutate(bar_id = row_number())  

#Draw ORA barchart for consensus upregulated & downregulated genes
sc_ora_p <- ggplot(sc_ora, aes(x = factor(bar_id), y = enrichmentRatio, fill = color)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  scale_y_continuous("Enrichment Ratio") +
  scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
  scale_x_discrete(labels = sc_ora$description) + 
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        legend.position = "none")

sc_ora_p #1000 x 1700


