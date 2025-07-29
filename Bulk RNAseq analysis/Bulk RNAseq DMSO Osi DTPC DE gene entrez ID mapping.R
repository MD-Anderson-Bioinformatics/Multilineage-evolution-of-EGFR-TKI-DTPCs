library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(org.Hs.eg.db)


##### HCC827 Osimertinib DTPC upregulated genes analysis #####

#Load in homo sapien genome database annotations
hs_ann <- org.Hs.eg.db

#Load in list of genes significantly upgregulated in HCC827 Osimertinib DTPCs

dtpc_up <- HCC827_osimertinib_DTPC_upregulated_significant$Gene

#Annotate each gene with Entrez ID required for HOMER analysis
dtpc_up_ann <- select(hs_ann, 
                      keys = dtpc_up, 
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
dtpc_up_ann

#Change second column name to 'gene_id'
colnames(dtpc_up_ann) <- c("Symbol", "gene_id")

dtpc_up_ann

#Write data to file
write.table(dtpc_up_ann, "C:/Users/bbmorris1/Desktop/HCC827 osimertinib DTPC upregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Note: 441 genes don't have matching Entrez ID, likely due to alias mismatch,
###will manually curate aliases and try again


#Re-mapping genes to Entrez IDs after manual curation

#Load in homo sapien genome database annotations
hs_ann <- org.Hs.eg.db

#Load in list of genes significantly upgregulated in HCC827 Osimertinib DTPCs

dtpc_up <-  HCC827_osimertinib_DTPC_upregulated_significant$Gene

#Annotate each gene with Entrez ID required for HOMER analysis
dtpc_up_ann <- select(hs_ann, 
                      keys = dtpc_up, 
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
dtpc_up_ann

#Change second column name to 'gene_id'
colnames(dtpc_up_ann) <- c("Symbol", "gene_id")

dtpc_up_ann

#Write data to file
write.table(dtpc_up_ann, "C:/Users/bbmorris1/Desktop/HCC827 osimertinib DTPC upregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)



##### HCC827 Osimertinib DTPC downregulated genes analysis #####

#Load in homo sapien genome database annotations
hs_ann <- org.Hs.eg.db

#Load in list of genes significantly downregulated in HCC827 Osimertinib DTPCs
dtpc_down <-  HCC827_osi_DTPC_significantly_downregulated_genes$Gene

#Annotate each gene with Entrez ID required for HOMER analysis
dtpc_down_ann <- select(hs_ann, 
                        keys = dtpc_down, 
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "SYMBOL")
dtpc_down_ann

#Change second column name to 'gene_id'
colnames(dtpc_down_ann) <- c("Symbol", "gene_id")

dtpc_down_ann

#Write data to file
write.table(dtpc_down_ann, "C:/Users/bbmorris1/Desktop/HCC827 osimertinib DTPC downregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)

###Note: 159 genes don't match Entrez IDs, likely due to alias mismatch
###Re-mapping these genes by hand


###Note: Only 3/159 genes don't have Entrez IDs after manual annotation




##### HCC4006 Osi DTPC RNAseq analysis #####

#Load in homo sapien genome database annotations
hs_ann <- org.Hs.eg.db

#Load in list of genes significantly upgregulated in HCC4006 osimertinib DTPCs

dtpc_up <-  HCC4006_osi_DTPC_upregulated_genes$Gene

#Annotate each gene with Entrez ID required for HOMER analysis
dtpc_up_ann <- select(hs_ann, 
                      keys = dtpc_up, 
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
dtpc_up_ann

#Change second column name to 'gene_id'
colnames(dtpc_up_ann) <- c("Symbol", "gene_id")

dtpc_up_ann

#Write data to file
write.table(dtpc_up_ann, "C:/Users/bbmorris1/Desktop/HCC4006 osi DTPC upregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Note: 77 genes don't have matching Entrez ID, likely due to alias mismatch,
###will manually curate aliases and try again


##### updated analysis (updating to upregulated FDR < 0.05) #####

#Merge list of HCC4006 osi DTPC upregulated genes (q < 0.05) with previous list of
#upregulated genes (q < 0.001) with annotated aliases
file1 <- HCC4006_osi_DTPC_upregulated_sig_0_05
file2 <- HCC4006_osi_DTPC_upregulated_sig_0_001_with_alias
matched <- join(file1, file2, type = "left", match = "all")
matched

write.table(matched, "C:/Users/bbmorris1/Desktop/HCC4006 osi DTPC upregulated sig 0.05 with alias.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Load in homo sapien genome database annotations
hs_ann <- org.Hs.eg.db

#Load in list of genes significantly upgregulated in HCC4006 osimertinib DTPCs
dtpc_up <-  HCC4006_osi_DTPC_upregulated_genes$Gene

#Annotate each gene with Entrez ID required for HOMER analysis
dtpc_up_ann <- select(hs_ann, 
                      keys = dtpc_up, 
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
dtpc_up_ann

#Change second column name to 'gene_id'
colnames(dtpc_up_ann) <- c("Symbol", "gene_id")

dtpc_up_ann

#Write data to file
write.table(dtpc_up_ann, "C:/Users/bbmorris1/Desktop/HCC4006 osi DTPC upregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Note: 199/6148 genes don't have matching Entrez ID, likely due to alias mismatch,
###will manually curate aliases and try again


###Note: only 2/199 genes (both lncRNAs) did not have Entrez IDs after matching aliases



##### HCC4006 osi DTPC analysis (downregulated FDR < 0.05) #####

#Load in homo sapien genome database annotations
hs_ann <- org.Hs.eg.db

#Load in list of genes significantly upgregulated in HCC4006 osimertinib DTPCs
dtpc_down <-  HCC4006_osi_DTPC_sig_downregulated_genes$Gene

#Annotate each gene with Entrez ID required for HOMER analysis
dtpc_down_ann <- select(hs_ann, 
                        keys = dtpc_down, 
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "SYMBOL")
dtpc_down_ann

#Change second column name to 'gene_id'
colnames(dtpc_down_ann) <- c("Symbol", "gene_id")

dtpc_down_ann

#Write data to file
write.table(dtpc_down_ann, "C:/Users/bbmorris1/Desktop/HCC4006 osi DTPC downregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Note: 173/5727 genes don't have matching Entrez ID, likely due to alias mismatch,
###will manually curate aliases and try again


###Note: 6/173 genes don't have matching Entrez ID after manual curation



##### H1975 osi DTPC upregulaed RNAseq analysis (FDR < 0.05) #####

#Load in homo sapien genome database annotations
hs_ann <- org.Hs.eg.db

#Load in list of genes significantly upgregulated in H1975 Osimertinib DTPCs

dtpc_up <- H1975_osi_DTPC_sig_upregulated_genes$Gene

#Annotate each gene with Entrez ID required for HOMER analysis
dtpc_up_ann <- select(hs_ann, 
                      keys = dtpc_up, 
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL")
dtpc_up_ann

#Change second column name to 'gene_id'
colnames(dtpc_up_ann) <- c("Symbol", "gene_id")

dtpc_up_ann

#Write data to file
write.table(dtpc_up_ann, "C:/Users/bbmorris1/Desktop/H1975 osi DTPC sig upregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Note: 239/6549 genes don't have matching Entrez ID, likely due to alias mismatch,
###will manually curate aliases and try again


###Note: 5 genes didn't have matching Entrez IDs after manual curation



##### H1975 osi DTPC downregulated RNAseq analysis (FDR < 0.05) #####

#Load in homo sapien genome database annotations
hs_ann <- org.Hs.eg.db

#Load in list of genes significantly upgregulated in H1975 Osimertinib DTPCs
dtpc_down <- H1975_osi_DTPC_sig_downregulated_genes$Gene

#Annotate each gene with Entrez ID required for HOMER analysis
dtpc_down_ann <- select(hs_ann, 
                        keys = dtpc_down, 
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "SYMBOL")
dtpc_down_ann

#Change second column name to 'gene_id'
colnames(dtpc_down_ann) <- c("Symbol", "gene_id")

dtpc_down_ann

#Write data to file
write.table(dtpc_down_ann, "C:/Users/bbmorris1/Desktop/H1975 osi DTPC sig downregulated genes annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###Note: 175/6547 genes don't have matching Entrez ID, likely due to alias mismatch,
###will manually curate aliases and try again


###Note: 6/175 genes didn't have matching Entrez IDs after manual curation










