library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(org.Hs.eg.db)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)



##### EGFR TKI DTPC shared upregulated genes (with all FDR < 0.05) #####

### HCC827, HCC4006, H1975 all shared ###

#Merge HCC827 osi DTPC upregulated genes with HCC4006 osi DTPC upregulated genes
#by Entrez ID
file1 <- HCC827_osimertinib_DTPC_upregulated_genes_entrez_id
file2 <- HCC4006_osi_DTPC_upregulated_genes_entrez_id
matched <- merge(file1, file2, by = "gene_id")
matched

#Merge HCC827 + HCC4006 shared upregulated genes with H1975 osi DTPC upregulated genes
#by Entrez ID

file3 <- H1975_osi_DTPC_sig_upregulated_genes_entrez_ids
matched_f <- merge(matched, file3, by = "gene_id")
matched_f

write.table(matched_f, "C:/Users/bbmorris1/Desktop/HCC827 HCC4006 H1975 osi DTPC shared upregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Convert entrez IDs back to gene symbols
library(annotate)

#Read in gene ids
dtpc_up <- as.character(HCC827_HCC4006_H1975_osi_DTPC_shared_upregulated_genes$gene_id)

dtpc_up_ann <- select(org.Hs.eg.db, 
                      keys = dtpc_up, 
                      column = c('SYMBOL', 'ENTREZID'), 
                      keytype = 'ENTREZID')

write.table(dtpc_up_ann, "C:/Users/bbmorris1/Desktop/HCC827 HCC4006 H1975 osi DTPC shared upregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)



### HCC827 HCC4006 H1975 shared upregulated genes venn diagram and euler plots ###

library(ggvenn)

#Create list of each model's osi DTPC significantly upregulated genes
up_d <- list('HCC827' = HCC827_osimertinib_DTPC_upregulated_genes_entrez_id$gene_id,
             'HCC4006' = HCC4006_osi_DTPC_upregulated_genes_entrez_id$gene_id,
             'H1975' = H1975_osi_DTPC_sig_upregulated_genes_entrez_ids$gene_id)

#Draw venn diagram
up_vd <- ggVennDiagram(up_d, label = "count", set_size = 3.2, edge_size = 0.5)

#Annotate with custom color palette
up_vd + scale_fill_distiller(palette = "RdBu", direction = -1)


#Draw euler plot for visualizing upregulated genes
library(eulerr)

euler_u <- euler(c("HCC827" = 1512, "HCC4006" = 1756, "H1975" = 2106, 
                   "HCC827&HCC4006" = 779, "HCC4006&H1975" = 1815, "HCC827&H1975" = 816,
                   "HCC827&HCC4006&H1975" = 1792),
                 shape = "ellipse")

up_e_p <- plot(euler_u, quantities = TRUE)
up_e_p



##### EGFR TKI DTPC shared downregulated genes (with all FDR < 0.05) #####

### HCC827, HCC4006, H1975 all shared ###

#Merge HCC827 osi DTPC downregulated genes with HCC4006 osi DTPC downregulated genes
#by Entrez ID
file1 <- HCC827_osi_DTPC_significantly_downregulated_genes_entrez_id
file2 <- HCC4006_osi_DTPC_downregulated_genes_entrez_ids
matched <- merge(file1, file2, by = "gene_id")
matched


#Merge HCC827 + HCC4006 shared downregulated genes with H1975 osi DTPC downregulated genes
#by Entrez ID

file3 <- H1975_osi_DTPC_sig_downregulated_genes_entrez_ids
matched_f <- merge(matched, file3, by = "gene_id")
matched_f

write.table(matched_f, "C:/Users/bbmorris1/Desktop/HCC827 HCC4006 H1975 osi DTPC shared downregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Convert entrez IDs back to gene symbols
library(annotate)

#Read in gene ids
dtpc_down <- as.character(HCC827_HCC4006_H1975_osi_DTPC_shared_downregulated_genes$gene_id)

dtpc_down_ann <- select(org.Hs.eg.db, 
                        keys = dtpc_down, 
                        column = c('SYMBOL', 'ENTREZID'), 
                        keytype = 'ENTREZID')

write.table(dtpc_down_ann, "C:/Users/bbmorris1/Desktop/HCC827 HCC4006 H1975 osi DTPC shared downregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)


### HCC827 HCC4006 H1975 shared downregulated genes venn diagram and euler plots ###

library(ggvenn)

#Create list of each model's osi DTPC significantly downregulated genes
dwn_d <- list('HCC827' = HCC827_osi_DTPC_significantly_downregulated_genes_entrez_id$gene_id,
              'HCC4006' = HCC4006_osi_DTPC_downregulated_genes_entrez_ids$gene_id,
              'H1975' = H1975_osi_DTPC_sig_downregulated_genes_entrez_ids$gene_id)

#Draw venn diagram
dwn_vd <- ggVennDiagram(dwn_d, label = "count", set_size = 3.2, edge_size = 0.5)

#Annotate with custom color palette
dwn_vd + scale_fill_distiller(palette = "RdBu", direction = -1)


#Draw euler plot for visualizing downregulated genes
library(eulerr)

euler_dwn <- euler(c("HCC827" = 1477, "HCC4006" = 1640, "H1975" = 1909, 
                     "HCC827&HCC4006" = 520, "HCC4006&H1975" = 1583, "HCC827&H1975" = 1072,
                     "HCC827&HCC4006&H1975" = 1968),
                   shape = "ellipse")

dwn_e_p <- plot(euler_dwn, quantities = TRUE)
dwn_e_p

