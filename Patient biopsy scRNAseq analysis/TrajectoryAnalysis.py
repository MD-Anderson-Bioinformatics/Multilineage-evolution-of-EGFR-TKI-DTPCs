#!/usr/bin/env python
# coding: utf-8


import palantir
import scanpy as sc
import numpy as np
import pandas as pd

# Plotting
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


#read in subset
adata = sc.read('EGFRSub_wCyto.h5ad')


metadirectory = "/Le_Vokes_sample_collection/"


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

#plot cytotrace pseudotime

ct_TN = adata[adata.obs['Context'].isin(['TN'])].obs['Specimen_ID'].value_counts().index.tolist()
ct_RD=adata[adata.obs['Context'].isin(['RD'])].obs['Specimen_ID'].value_counts().index.tolist()
ct_PD=adata[adata.obs['Context'].isin(['PD'])].obs['Specimen_ID'].value_counts().index.tolist()

markerscore=['ct_pseudotime'] 
for score in markerscore:
    Samples = adata.obs['Specimen_ID'].value_counts().index.tolist()

    tindex = []
    order = []
    for  n in ct_TN:
        order.append(n)
        tindex.append(0)
    for  n in ct_RD:
        order.append(n)
        tindex.append(1)
    for  n in ct_PD:
        order.append(n)
        tindex.append(2)


    nsamples = len(Samples)#24
    PLOT = []
    porder = []
    pindex = []
    for i in range(0,len(order)):
        if order[i] in Samples[0:nsamples]:
            porder.append(order[i])
            pindex.append(tindex[i])
            
            PLOT.append(adata[adata.obs['Specimen_ID'].isin([order[i]]),:].obs[score].tolist())
           

    fig,ax1=plt.subplots(figsize=(22,6))
    ax1.set(
        axisbelow=True,  # Hide the grid behind plot objects
        title= score,
        xlabel='Specimen ID',
        ylabel='Score',
    )
    plt.grid(True)
    
    bp = ax1.boxplot(PLOT)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')

    box_colors = ['tab:blue', 'tab:orange','tab:green','tab:red']
    box_colors = ['#33a02c','#a6cee3','#1f78b4','#b2df8a']
    num_boxes = len(PLOT)
    medians = np.empty(num_boxes)
    for i in range(num_boxes):
        box = bp['boxes'][i]
        box_x = []
        box_y = []
        for j in range(5):
            box_x.append(box.get_xdata()[j])
            box_y.append(box.get_ydata()[j])
        box_coords = np.column_stack([box_x, box_y])
       
        ax1.add_patch(Polygon(box_coords, facecolor=box_colors[pindex[i]]))
        
        med = bp['medians'][i]
        median_x = []
        median_y = []
        for j in range(2):
            median_x.append(med.get_xdata()[j])
            median_y.append(med.get_ydata()[j])
            ax1.plot(median_x, median_y, 'k')
        medians[i] = median_y[0]

    ax1.set_xticklabels(porder,rotation=45, fontsize=8)


    fig.text(0.96, 0.98, f'TN',
             bbox=dict(facecolor=box_colors[0],alpha=0.8,linewidth=0), color='black', weight='demi',
             size='x-small')
    fig.text(0.96, 0.94, '    RD   ',
             backgroundcolor=box_colors[1],
             color='black', weight='demi', size='x-small')
    fig.text(0.96, 0.9, '   PD    ', bbox=dict(facecolor=box_colors[2],alpha=0.8,linewidth=0),color='black', weight='demi',
             size='x-small')




# # Make UMAP

adata.X = adata.raw.X


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5)


sc.pp.scale(adata, max_value=10)


# Set seed
intialization = 3120
sc.pp.pca(adata, random_state=intialization, n_comps = 100 ,svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs= 100, log=True, show = True)#, save = "PCA_variance.png")
sc.pl.pca(adata,color=["BatchLegend"])#, save = "PCA.png")


sc.pp.neighbors(adata, random_state=intialization, n_neighbors=30, n_pcs=50,  method = 'umap', metric = 'euclidean')


sc.tl.umap(adata, random_state=intialization,spread = 1,min_dist=0.5)


del adata.uns['Specimen_ID_colors']


sc.pl.umap(adata,color=['Specimen_ID'],show = True,title = "",frameon=False, legend_fontsize=10, legend_fontoutline=2,na_in_legend=False)


# # Diffusion maps


# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(adata, n_components =10)


#choose # of diffusion components based on eigengap
ms_data = palantir.utils.determine_multiscale_space(adata)


# # Create force directed layout for visualization

#choose inital positions from UMAP
init_pos = adata.obsm['X_umap'].tolist()


sc.tl.draw_graph(adata,layout='fa',obsp='DM_Kernel',maxiter=500,init_pos='X_umap')

#color by specimen id
sc.pl.draw_graph(adata,color=['Specimen_ID'],
                layout='fa',frameon=False,
                edges=False,
                title="")

#color by context
cdict={}
cdict['PD'] = 'tab:blue'
cdict['RD'] = 'tab:orange'
cdict['TN'] = 'tab:green'

sc.pl.draw_graph(adata,color=['Context'],
                layout='fa',frameon=False,
                edges=False, palette=cdict,
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


#find unique endpoints for ensemble of waypoint values:
nw = 1000
endpoints=[]
for i in range(0,26):
    print(nw)
    print()
    pr_res = palantir.core.run_palantir(
    adata, 'CTGCTGTGTTCAGTAC-1-JH064', num_waypoints=nw, knn=50,use_early_cell_as_start =False,n_jobs=20
    )
    nep = adata.obsm['palantir_fate_probabilities'].columns.tolist()
    for ep in nep:
        if ep not in endpoints:
            endpoints.append(ep)
    waypoints = pr_res.waypoints
    palantir.plot.highlight_cells_on_umap(adata,waypoints,embedding_basis='X_draw_graph_fa')
    plt.show()
    palantir.plot.plot_palantir_results(adata,embedding_basis='X_draw_graph_fa')
    plt.show()
    nw = nw+250
    print()


#view endpoints from previous runs
palantir.plot.highlight_cells_on_umap(adata,endpoints,embedding_basis='X_draw_graph_fa')
plt.show()


#rename the terminal state cell labels for plot label clarity
terminal_states = pd.Series(
    ['Lung-Tumor-10','JH386','JH067','JH304','Lung-Tumor-7','JH033'],
    index=['CCGGACATCTGAGCAT-1-Lung-tumor-10', 'CGTTCTGAGTGCGATG-1-JH386', 'TCACAAGGTTCCCTTG-1-JH067', 'TGTTCCGGTACACCGC-1-JH304', 'TGACTAGGTGCAACTT-1-Lung-Tumor-7', 'TCACGAAGTGTGAATA-1-JH033']
)


#run palintir w updated labels
pr_res = palantir.core.run_palantir(
    adata, 'CTGCTGTGTTCAGTAC-1-JH064', num_waypoints=5000,terminal_states=terminal_states, knn=50,use_early_cell_as_start =False
)


#visualize results
palantir.plot.plot_palantir_results(adata,embedding_basis='X_draw_graph_fa')
plt.show()


#select branch cells
masks = palantir.presults.select_branch_cells(adata, q=1e-5, eps=1e-2)


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
branch = fatep.columns.tolist()

#make the same plots as in the plot_branch_selection() function but with KRT17 as cell color
genedf = sc.get.obs_df(adata, keys=["KRT17"],use_raw=True)


KRT17 = genedf['KRT17'].tolist()


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




