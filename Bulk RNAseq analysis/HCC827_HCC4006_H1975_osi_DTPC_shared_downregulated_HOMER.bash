module load perl/5.36.1
module load homer/4.11

configureHomer.pl -list

findMotifs.pl HCC827_HCC4006_H1975_osi_DTPC_shared_downregulated_genes_entrez_ids.txt human /Downregulated_genes_out_dir/ -len 8,10,12