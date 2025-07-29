library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(SeuratExtend)
library(SeuratDisk)
library(sctransform)
library(scran)
library(ggpubr)
library(harmony)
library(patchwork)
library(magrittr)


##### HCC827, HCC4006, H1975 consensus DTPC signature analysis #####

###Upregulated genes

#Merge HCC827, HCC4006, and H1975 DTPC upregulated genes to identify conserved upregulated genes @ DTPC
file1 <- HCC827_DTPC_padj_sig_upregulated_genes
file2 <- HCC4006_DTPC_padj_sig_upregulated_genes
matched <- merge(file1, file2, by = "Gene")
matched

file3 <- matched
file4 <- H1975_DTPC_padj_sig_upregulated_genes
matched_2 <- merge(file3, file4, by = "Gene")
matched_2

write.table(matched_2, "C:/Users/bbmorris1/Desktop/HCC827 HCC4006 H1975 scRNAseq consensus DTPC upregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Downregulated genes

#Merge HCC827, HCC4006, and H1975 DTPC downregulated genes to identify conserved downregulated genes @ DTPC
file1 <- HCC827_DTPC_padj_sig_downregulated_genes
file2 <- HCC4006_DTPC_padj_sig_downregulated_genes
matched <- merge(file1, file2, by = "Gene")
matched

file3 <- matched
file4 <- H1975_DTPC_padj_sig_downregulated_genes
matched_2 <- merge(file3, file4, by = "Gene")
matched_2

write.table(matched_2, "C:/Users/bbmorris1/Desktop/HCC827 HCC4006 H1975 scRNAseq conserved DTPC downregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Upregulated genes Euler plot (and individual pair wise overlaps)

#Read in data
hcc827_u <- HCC827_DTPC_padj_sig_upregulated_genes
hcc4006_u <- HCC4006_DTPC_padj_sig_upregulated_genes
h1975_u <- H1975_DTPC_padj_sig_upregulated_genes


#Merge HCC827 & HCC4006 Osi DTPC upregulated gene lists
hcc827_hcc4006_u <- merge(hcc827_u, hcc4006_u, by = "Gene")
hcc827_hcc4006_u
write.table(hcc827_hcc4006_u, "C:/Users/bbmorris1/Desktop/HCC827 HCC4006 Osi DTPC shared upregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Merge HCC827 & H1975 Osi DTPC upregulated gene lists
hcc827_h1975_u <- merge(hcc827_u, h1975_u, by = "Gene")
hcc827_h1975_u

write.table(hcc827_h1975_u, "C:/Users/bbmorris1/Desktop/HCC827 H1975 Osi DTPC shared upregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Merge HCC4406 & H1975 Osi DTPC upregulated gene lists
hcc4006_h1975_u <- merge(hcc4006_u, h1975_u, by = "Gene")
hcc4006_h1975_u

write.table(hcc4006_h1975_u, "C:/Users/bbmorris1/Desktop/HCC4006 H1975 Osi DTPC shared upregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Plot shared upregulated gene Euler plot
library(eulerr)

euler_u <- euler(c("HCC827" = 2176, "HCC4006" = 6261, "H1975" = 7743, 
                   "HCC827&HCC4006" = 1079, "HCC4006&H1975" = 1660, "HCC827&H1975" = 1281,
                   "HCC827&HCC4006&H1975" = 694),
                 shape = "ellipse")

up_e_p <- plot(euler_u, quantities = TRUE)
up_e_p



###Downregulated genes Euler plot (and individual pair wise overlaps)

#Read in data
hcc827_d <- HCC827_DTPC_padj_sig_downregulated_genes
hcc4006_d <- HCC4006_DTPC_padj_sig_downregulated_genes
h1975_d <- H1975_DTPC_padj_sig_downregulated_genes

#Merge HCC827 & HCC4006 Osi DTPC downregulated gene lists
hcc827_h4006_d <- merge(hcc827_d, hcc4006_d, by = "Gene")
hcc827_h4006_d

write.table(hcc827_h4006_d, "C:/Users/bbmorris1/Desktop/HCC827 HCC4006 Osi DTPC shared downregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)  


#Merge HCC827 & H1975 Osi DTPC downregulated gene lists
hcc827_h1975_d <- merge(hcc827_d, h1975_d, by = "Gene")
hcc827_h1975_d

write.table(hcc827_h1975_d, "C:/Users/bbmorris1/Desktop/HCC827 H1975 Osi DTPC shared downregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Merge HCC4006 & H1975 Osi DTPC downregulated gene lists
hcc4006_h1975_d <- merge(hcc4006_d, h1975_d, by = "Gene")
hcc4006_h1975_d

write.table(hcc4006_h1975_d, "C:/Users/bbmorris1/Desktop/HCC4006 H1975 Osi DTPC shared downregulated genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Plot shared downregulated gene Euler plot
euler_d <- euler(c("HCC827" = 11407, "HCC4006" = 11201, "H1975" = 5028, 
                   "HCC827&HCC4006" = 8730, "HCC4006&H1975" = 2868, "HCC827&H1975" = 3131,
                   "HCC827&HCC4006&H1975" = 2337),
                 shape = "ellipse")

down_e_p <- plot(euler_d, quantities = TRUE)
down_e_p



### Entrez ID mapping ###
library(org.Hs.eg.db)

#Read in gene ids
dtpc_consensus_up <- HCC827_HCC4006_H1975_Osi_DTPC_shared_upregulated_gene_names$Gene

dtpc_consensus_up_ann <- mapIds(org.Hs.eg.db, 
                                keys = dtpc_consensus_up, 
                                column = c('ENTREZID'), 
                                keytype = 'SYMBOL')
dtpc_consensus_up_ann

write.table(dtpc_consensus_up_ann, "C:/Users/bbmorris1/Desktop/HCC827 HCC4006 H1975 Osi DTPC shared upregulated genes annotated.txt", sep = "\t", quote = FALSE)

#Note: needed to manually match 6 genes to updated Entrez IDs, 79 entries lack entrez ids (are unknown transcripts/lncRNAs)



