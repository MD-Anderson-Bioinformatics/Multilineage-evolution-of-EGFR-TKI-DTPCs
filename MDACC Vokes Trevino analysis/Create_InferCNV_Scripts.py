#!/usr/bin/env python
# coding: utf-8

# In[1]:


#patient numbers for analysis
samples = [patient_01, patient_02, patient_03, patient_04, patient_05, patient_06, patient_07, patient_08, patient_09, patient_10, 
           patient_11, patient_12, patient_13, patient_14, patient_15, patient_16, patient_17, patient_18, patient_19, patient_20, 
           patient_21, patient_22, patient_23, patient_24, patient_25, patient_26, patient_27, patient_28, patient_29, patient_30,
           patient_31, patient_32, patient_33, patient_34, patient_35, patient_36, patient_37, patient_38, patient_39, patient_40,
           patient_41, patient_42, patient_43, patient_44, patient_45, patient_46, patient_47, patient_48, patient_49, patient_50,
           patient_51, patient_52, patient_53, patient_54, patient_55]

# In[2]:


#make strings
new = []
for s in samples:
    new.append(str(s))
samples = new


# In[3]:


#add ATAC
add = [patient_08]
for a in add:
    samples.append(str(a)+'_ATAC')


# In[4]:


print(samples)


# In[5]:


print(len(samples))


# In[6]:


#cluster into 4 groups for initial analysis
groups = [4]*len(samples)


# In[7]:


print(groups)


# In[8]:


print(len(groups))


# In[12]:


#make scripts for HPC
for i in range(0,len(samples)):
    text_file = open("RinferCNV_"+ str(i+1) + ".lsf", "w")
    text_file.write("#BSUB -J R-InferCNV-%s\n#BSUB -W 240:00\n" % str(i+1))
    text_file.write("#BSUB -o /rsrch3/home/thor_hn_med_onc_rsch/sgtrevino1/\n")
    text_file.write("#BSUB -e /rsrch3/home/thor_hn_med_onc_rsch/sgtrevino1/\n")
    text_file.write("#BSUB -cwd /rsrch3/home/thor_hn_med_onc_rsch/sgtrevino1/2023_02_08_batch/ScanPy/InferCNVscripts/\n")
    text_file.write("#BSUB -q e80long\n")
    text_file.write("#BSUB -n 10\n")
    text_file.write("#BSUB -M 512\n")
    text_file.write("#BSUB -R rusage[mem=512]\n")
    text_file.write("#BSUB -N\n")
    text_file.write("#BSUB -B\n")
    text_file.write("#BSUB -u sgtrevino1@mdanderson.org\n\n")
    text_file.write("module load R/3.6.0\n")
    
    text_file.write("R CMD BATCH RinferCNV_run%d.R Routput%d.out\n" % (i+1,i+1))
    text_file.close()
    
    text_file = open("RinferCNV_run"+ str(i+1) + ".R", "w")
    text_file.write("library('infercnv')\n\n")
    text_file.write("# create the infercnv object\n")
    text_file.write("infercnv_obj = CreateInfercnvObject(raw_counts_matrix='/rsrch3/home/thor_hn_med_onc_rsch/sgtrevino1/2023_02_08_batch/ScanPy/INFERCNV_R_PATIENT/INFERCNV_input_%s_MYELOID.csv',\n" % samples[i])
    text_file.write("\t\t\tannotations_file='/rsrch3/home/thor_hn_med_onc_rsch/sgtrevino1/2023_02_08_batch/ScanPy/INFERCNV_R_PATIENT/annotation_%s_MYELOID.txt',\n" % samples[i])
    text_file.write('\t\t\tdelim="\\%s",\n' % 't')
    text_file.write("\t\t\tgene_order_file='/rsrch3/home/thor_hn_med_onc_rsch/sgtrevino1/geneorder.txt',\n")
    text_file.write('\t\t\tref_group_names=c("normal"))\n\n')
    text_file.write('out_dir="/rsrch4/home/thor_hn_med_onc/nvokes/Vokes_Le_scRNAseq_treatment_samples/INFERCNV_Le_Vokes_4/InferCNV_%s_out/"' % samples[i])
    text_file.write("# perform infercnv operations to reveal cnv signal\n")
    text_file.write("infercnv_obj = infercnv::run(infercnv_obj,\n\t\t\tcutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics\n\t\t\tout_dir=out_dir,\n \t\t\tcluster_by_groups=FALSE,\n \t\t\tdenoise=TRUE,\n\t\t\tk_obs_groups = %d,\n\t\t\tscale_data = TRUE,\n\t\t\tHMM=FALSE,\n\t\t\tnum_threads=10,\n\t\t\thclust_method = 'ward.D2',\n\t\t\tno_prelim_plot=TRUE,\n\t\t\t)\n" % groups[i])
    text_file.close()
    


# In[ ]:




