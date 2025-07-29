#!/usr/bin/env python
# coding: utf-8


import palantir
import scanpy as sc
import numpy as np
import pandas as pd

# Plotting
import matplotlib.pyplot as plt


#read in subset
adata = sc.read('/Vokes_Le_scRNA_EGFR_5_25/CodeForBen/EGFRSub.h5ad')


metadirectory = "/Le_Vokes_sample_collection/2025_05_EGFR_annotations/"


#get metadata
Meta_df = pd.read_excel(metadirectory + 'Sequenced_Samples_Sample_Database_EGFR.xlsx')


#Add updated Meta data labels
MetaID = Meta_df['specimen_id'].tolist()
SampleID = adata.obs['Specimen_ID'].value_counts().index.tolist()


for ID in SampleID:
    if ID not in MetaID:
        print(ID)


ResMechMeta = Meta_df['EGFR_resistance_mechanism'].tolist()
annot = {}
for i in range(0,len(MetaID)):
    annot[MetaID[i]] = ResMechMeta[i] 
adata.obs['Res_mech'] = adata.obs['Specimen_ID'].map(annot).astype('category')


ResMechMeta = Meta_df['specimen_histology'].tolist()
annot = {}
for i in range(0,len(MetaID)):
    annot[MetaID[i]] = ResMechMeta[i] 
adata.obs['Histology'] = adata.obs['Specimen_ID'].map(annot).astype('category')


ResMechMeta = Meta_df['Sample_organ'].tolist()
annot = {}
for i in range(0,len(MetaID)):
    annot[MetaID[i]] = ResMechMeta[i] 
adata.obs['Location'] = adata.obs['Specimen_ID'].map(annot).astype('category')


ResMechMeta = Meta_df['Context'].tolist()
annot = {}
for i in range(0,len(MetaID)):
    annot[MetaID[i]] = ResMechMeta[i] 
adata.obs['Context'] = adata.obs['Specimen_ID'].map(annot).astype('category')


ResMechMeta = Meta_df['Treatment_class'].tolist()
annot = {}
for i in range(0,len(MetaID)):
    annot[MetaID[i]] = ResMechMeta[i] 
adata.obs['Treatment_class'] = adata.obs['Specimen_ID'].map(annot).astype('category')


adata.obs['Specimen_ID'].value_counts()


# # Make UMAP

adata.X = adata.raw.X


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5)


sc.pp.scale(adata, max_value=10)


# Set seed
intialization = 3120
sc.pp.pca(adata, random_state=intialization, n_comps = 100 ,svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs= 100, log=True, show = True)#, save = "PCA_variance.png")
sc.pl.pca(adata,color=["BatchLegend"])#, save = "PCA.png")


sc.pp.neighbors(adata, random_state=intialization, n_neighbors=15, n_pcs=50,  method = 'umap', metric = 'euclidean')


sc.tl.umap(adata, random_state=intialization,spread = 1,min_dist=0.5)


del adata.uns['Specimen_ID_colors']


sc.pl.umap(adata,color=['Specimen_ID'],show = True,title = "",frameon=False, legend_fontsize=10, legend_fontoutline=2,na_in_legend=False)


# # Diffusion maps

#palantir.utils.run_pca(adata,n_components=300,use_hvg=False)


# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(adata, n_components =10)


#choose # of diffusion components based on eigengap
ms_data = palantir.utils.determine_multiscale_space(adata)


# # Create force directed layout for visualization

#choose inital positions from UMAP
init_pos = adata.obsm['X_umap'].tolist()


sc.tl.draw_graph(adata,layout='fa',obsp='DM_Kernel',maxiter=500,init_pos='X_umap')


sc.pl.draw_graph(adata,color=['Specimen_ID'],
                layout='fa',frameon=False,
                edges=False,
                title="")


#run magic imputation


imputed_X = palantir.utils.run_magic_imputation(adata)


sc.pl.embedding(
    adata,
     basis='X_draw_graph_fa',
    layer="MAGIC_imputed_data",
    color=["KRT17"],
    frameon=False,
    cmap='viridis'
)
plt.show()


#lognormalized for comparison
sc.pl.embedding(
    adata,
    basis='X_draw_graph_fa',
    #layer="MAGIC_imputed_data",
    color=["KRT17"],
    frameon=False,
    ncols=5,
    cmap='viridis',
    #vmax=0.25

)


# # Palintir

#choose a root cell from TN sample JH064
roots = adata[adata.obs['Specimen_ID'].isin(['JH064']),:].obs.index.tolist()[0:1]


#plot root on trajectory map
palantir.plot.highlight_cells_on_umap(adata, roots,embedding_basis='X_draw_graph_fa')
plt.show()


#view endpoints from previous runs
palantir.plot.highlight_cells_on_umap(adata, ['CTAACTTTCACTCCTG-1-JH304','ATCCACCTCTTGGGTA-1-Lung-Tumor-7','CATCGGGAGTCGAGTG-1-JH038','TTCTCCTGTCGCATAT-1-JH033','TCACAAGGTTCCCTTG-1-JH067'],embedding_basis='X_draw_graph_fa')
plt.show()


#rename the terminal state cell labels for plot label clarity
terminal_states = pd.Series(
    ['JH304','Lung-Tumor-7','JH038','JH033','JH067'],
    index=['CTAACTTTCACTCCTG-1-JH304','ATCCACCTCTTGGGTA-1-Lung-Tumor-7','CATCGGGAGTCGAGTG-1-JH038','TTCTCCTGTCGCATAT-1-JH033','TCACAAGGTTCCCTTG-1-JH067'],
)


#run palintir 
pr_res = palantir.core.run_palantir(
    adata, 'ACACCCTGTAGAAGGA-1-JH064', num_waypoints=1500,terminal_states=terminal_states
)


#visualize results
palantir.plot.plot_palantir_results(adata,embedding_basis='X_draw_graph_fa')
plt.show()


#select branch cells
masks = palantir.presults.select_branch_cells(adata, q=1e-4, eps=1e-4)


#visualize branches
fig = palantir.plot.plot_branch_selection(adata,embedding_basis='X_draw_graph_fa')
plt.show()


#compute gene trends on trajectory 
gene_trends = palantir.presults.compute_gene_trends(
    adata,
    expression_key="MAGIC_imputed_data",
)


#plot KRT17 along each trajectory in pseudotime
genes = ["KRT17"]
palantir.plot.plot_gene_trends(adata, genes)
plt.show()


#plot KRT17 along each trajectory as a heatmap
genes = ["KRT17"]
palantir.plot.plot_gene_trend_heatmaps(adata, genes,basefigsize=(7,1.0),scaling='z-score')
plt.show()


#Custom Plots 


#make the same plots as in the plot_branch_selection() function but with KRT17 as cell color
genedf = sc.get.obs_df(adata, keys=["KRT17"],use_raw=True)


KRT17 = genedf['KRT17'].tolist()


branch = ['Lung-Tumor-7','JH038','JH304','JH067','JH033']



for i in range(0,len(masks[0])):
    glist=[]
    for j in range(0,len(masks[:,i])):
        if masks[j,i] == False:
            glist.append(np.nan)
        else:
            glist.append(KRT17[j])
    genedf[branch[i]+'_KRT17'] = glist


del genedf['KRT17']


columns = genedf.columns.tolist()


for i in range(0,len(columns)):
    adata.obs[columns[i]] = genedf[columns[i]]


branch = ['Lung-Tumor-7','JH038','JH304','JH067','JH033']


#trajectory maps
for i in range(0,len(branch)):
    title = "Branch " + branch[i]
    data = branch[i] + "_KRT17"
    sc.pl.embedding(
    adata,
    basis='X_draw_graph_fa',
    color=[data],
    frameon=False,
    cmap='viridis',
    title=title,

)


branch = ['Lung-Tumor-7','JH038','JH304','JH067','JH033']
fatep = adata.obsm['palantir_fate_probabilities']
pt = adata.obs['palantir_pseudotime'].tolist()


#fate probability plots
for i in range(0,len(branch)):
    x = pt
    y = fatep[branch[i]].tolist()
    title = "Branch " + branch[i]
    data = branch[i] + "_KRT17"
    plt.figure(figsize=(10, 3))
    # Get a colormap and set the color for bad values (NaN)
    cmap = plt.cm.viridis
    cmap.set_bad('lightgrey',alpha=0.1) # Set NaNs to red
    scatter=plt.scatter(x,y,s=9.0,c=adata.obs[data].tolist(),cmap=cmap, plotnonfinite=True)
    plt.colorbar(scatter, label='KRT17')
    plt.title(title,fontsize=8)
    plt.ylabel("Fate Probability",fontsize=8)
    plt.xlabel("Pseudotime",fontsize=8)
    plt.show()




