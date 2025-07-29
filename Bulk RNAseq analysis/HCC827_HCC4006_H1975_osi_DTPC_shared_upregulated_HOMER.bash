#!/bin/bash

#BSUB -W 3:00
#BSUB -o /rsrch4/home/thor_hn_med_onc_rsch/bbmorris1/HOMER/EGFR_DTPC/DTPC_shared/Upregulated_genes
#BSUB -e /rsrch4/scratch/thor_hn_med_onc_rsch/bbmorris1/error_logs
#BSUB -cwd /rsrch4/home/thor_hn_med_onc_rsch/bbmorris1/HOMER/EGFR_DTPC/DTPC_shared
#BSUB -q e80short
#BSUB -M 950
#BSUB -R rusage[mem=950]
#BSUB -u bbmorris1@mdanderson.org
#BSUB -J HCC827_HCC4006_H1975_osi_DTPC_shared_upregulated_HOMER

module load perl/5.36.1
module load homer/4.11

configureHomer.pl -list

findMotifs.pl HCC827_HCC4006_H1975_osi_DTPC_shared_upregulated_genes_entrez_ids.txt human /rsrch4/home/thor_hn_med_onc_rsch/bbmorris1/HOMER/EGFR_DTPC/DTPC_shared/Upregulated_genes -len 8,10,12