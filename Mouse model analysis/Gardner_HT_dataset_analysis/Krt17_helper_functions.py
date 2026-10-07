# Python file for defining helper functions for Krt17 Analysis:
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hc

import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib
from matplotlib import rcParams
from adjustText import adjust_text

from scipy.stats import zscore
from copy import deepcopy
from kneed.knee_locator import KneeLocator
from sklearn.metrics.cluster import adjusted_rand_score
from skimage import filters
from scipy import stats

import scipy
import sys
import os
import glob
import anndata
import itertools
import phenograph
import sklearn
import time
import pickle

from os import path
from pathlib import Path
from scipy.ndimage import convolve
from statsmodels.nonparametric.smoothers_lowess import lowess

output_dir_base = '/workdir/varmus_single_cell/merged_pipeline_out/RD_Krt17_0225/Revisions_0626/'
output_dir = output_dir_base
git_dir = './scripts/'

# Helper function to flatten nested lists:
def flatten(l):
    if all(isinstance(x, list) for x in l):
        return [item for sublist in l for item in sublist]
    elif all(isinstance(x, np.ndarray) for x in l):
        return [item for sublist in l for item in sublist]
    else:
        return l

# Helper function for a quick look at variable distributions:
def stratify(adata, var1, var2, normalize=False):
    for x in sorted(adata.obs[var1].unique()):
        print(x)
        print(adata[adata.obs[var1] == x].obs[var2].value_counts(normalize=normalize))
        print('\n', '#########################', '\n')
    
    
# run bash command in background if stdoutfile does not exist. If it does, output contents of stdoutfile.
# If force=true, overwrite and re-run command.
def run_in_background(command, stdoutfile, stderrfile="", force=False, quiet=False, 
                      wait=False):
    if force or not path.exists(stdoutfile):
        command = 'bash -c \'' + command + '\' > ' + stdoutfile
        if stderrfile:
            command = command + ' 2>' + stderrfile
        else:
            command = command + ' 2>&1'
        if not wait:
            command = command + ' &'
        print('calling ' + command + '\n')
        os.system(command)
    else:
        wait=True
    if wait:
        if not quiet:
            print("Output from stdout file " + stdoutfile)
            os.system(f'cat {stdoutfile}')
        else:
            print("Output from stdout file " + stdoutfile + " is suppressed")
        if stderrfile and path.exists(stderrfile):
            print("Output from stderr file " + stderrfile)
            os.system(f'cat {stdoutfile}')
            

# Helper function for filtering lowly expressed genes
def filter_genes(adata, minCell):
    hasCell = adata.layers['X'] > 0
    numCell = np.array(hasCell.sum(axis=0))
    tooSmall = [i < minCell for i in numCell]
    is_invalid = np.zeros(adata.shape[1], np.bool)
    is_invalid[np.where(tooSmall)[0]] = True
    print(f"Removing {np.sum(is_invalid)} out of {len(is_invalid)} genes with less than {minCell} cells")
    adata.var['n_cells_by_counts'] = numCell
    return is_invalid


# Define functions for vmin and vmax of genes expression:
def get_vmin(values):
    p=0.05
    vmin=float(np.quantile(values, p))
    if vmin<=0:
        vmin=0
    return float(vmin)

def get_vmax(values):
    p=0.95
    vmax=float(np.quantile(values, p))
    if vmax<=0:
        vmax=0.1
    return float(vmax)


# barplot function:
def cellType_frac_barplot(data, group_var, cellType_var, group_order=None, cmap=None, scaled=False, output_dir=output_dir, prefix='', figsize=(10,8)):

    # Calculate group fractions:
    df_frac = data.obs.groupby([group_var, cellType_var], observed=False).size().unstack(fill_value=0)
    df_frac = df_frac.div(df_frac.sum(axis=1), axis=0).T

    # Define group order if not provided:
    if group_order==None:
        group_order=data.obs[group_var].value_counts().index.tolist()
    df_frac = df_frac[group_order]
    
    # Scale groups by the median sample size if requested:
    if scaled:
        sample_size_scale = data.obs[group_var].value_counts().loc[df_frac.columns.tolist()].values / \
                            int(data.obs[group_var].value_counts().median())
        df_frac = df_frac * sample_size_scale

    # Format dataframe:
    arc_sample_frac_df = pd.DataFrame()
    arc_list = []
    sample_list = []
    frac_list = []
    
    for col in df_frac:
        arc_list = arc_list + df_frac.index.tolist()
        sample_list = sample_list + [col for x in range(len(df_frac))]
        frac_list = frac_list + df_frac[col].values.tolist()
        
    arc_sample_frac_df['Arc'] = arc_list
    arc_sample_frac_df['condition'] = sample_list
    arc_sample_frac_df['Frac'] = frac_list
    tmp = arc_sample_frac_df.pivot(index='condition', columns='Arc')
    tmp.index.name = None
    tmp.columns.name = None
    tmp.columns = tmp.columns.droplevel()
    tmp = tmp.loc[group_order,:]

    # Create batplot:
    matplotlib.style.use('default')
    fig, ax = plt.subplots(figsize=figsize)
    sns.despine()

    if cmap!=None:
        tmp.plot.bar(stacked=True, ax=ax, 
                     color=[cmap[x] for x in tmp.columns.tolist()],
                     edgecolor='black',
                     linewidth=1.25,
                     width=0.55)
    else:
        print('No cmap provided, using random colors.')
        tmp.plot.bar(stacked=True, ax=ax, 
                     #color=[cmap[x] for x in tmp.columns.tolist()],
                     edgecolor='black',
                     linewidth=1.25,
                     width=0.55)
    
    ax.legend(bbox_to_anchor=(1.01, 1.01))

    # Save figure to png:
    cellType_var_save = cellType_var.replace(' ', '_').replace('/', '_')
    group_var_save = group_var.replace(' ', '_').replace('/', '_')
    
    if scaled:
        fn = f'/{prefix}{cellType_var_save}_by_{group_var_save}_Barplot_scaled.svg'
    else:
        fn = f'/{prefix}{cellType_var_save}_by_{group_var_save}_Barplot.svg'
        
    plt.tight_layout()
    print(f'Figure saved : {output_dir}/{fn}')
    plt.savefig(f'{output_dir}/{fn}', dpi=400, transparent=True, bbox_inches='tight')



### Functions to select optimal number of nPCs using linear approximation ###
# A linear deviation threshold of 1e-5 is the default, as this has empirically worked well for most cases
# Returns the optimal number of PCs to keep (sensitive to starting number of PCs)
# The optimize_dev_thr function can be used to test the threshold used here.

def optimize_pcs(adata, pca_name='pca', dev_thr=1e-5, save=True, fn='PCA_linear_deviation_optimization.png'):

    plt.style.use('default')
    
    if pca_name not in adata.uns.keys():
        exit(f'ERROR: {pca_name} not in adata.')
    
    x = range(len(adata.uns[pca_name]['variance_ratio']))
    y = np.cumsum(adata.uns[pca_name]['variance_ratio'])
    
    deviation = []
    slopes = []
    
    #dev_thr = 1e-5
    thr_hit=False
    
    opt_n_pcs_linear_dev=0
    
    fig = plt.subplots(figsize=(6,6))
    
    for i in range(len(x)-2):
        
        tmp_slope = (y[i+1] - y[i]) / (x[i+1] - x[i])
        tmp_int = y[i] - tmp_slope*x[i]
        
        plt.scatter(x, y, s=10, color='blue')
        
        y_line = [tmp_slope*x + tmp_int for x in x]
        
        tmp_dev = abs(y_line[i+2] - y[i+2])
        
        deviation.append(tmp_dev)
        
        slopes.append(tmp_slope)
        
        plt.plot(x, y_line, color='grey', linewidth=0.5)
        
        if (tmp_dev < dev_thr) and (thr_hit==False):
            #print(f'Optimal nPCs = {i+1}')
            #print(f'Optimal PCs index = {i}')
            #print(f'Deviation = {tmp_dev}')
            
            opt_n_pcs_linear_dev = i
            
            plt.scatter(x[i], y[i], s=15, color='red')
            plt.axvline(x[i], color='red', linestyle='--', linewidth=1)
            
            thr_hit=True
        
        plt.xlabel('PCs Ranked')
        plt.ylabel('Cumulative Variance')
        plt.title(f'Cumulative Explained Variance, i = {opt_n_pcs_linear_dev} | nPCs = {opt_n_pcs_linear_dev+1}')
        
        plt.ylim([y[0]*0.9, y[-1]*1.1])

    if save:
        fn = f'{output_dir}/{fn}'
        plt.tight_layout()
        print(fn)
        plt.savefig(fn, dpi=400)

    return opt_n_pcs_linear_dev+1
        
# Function to test different linear deviation thresholds for optimizing nPCs

# We want to choose the deviation threshold that is the inflection point in the nPC -by- threshold plot. 
# The number of PCs is stable for a wide range of thresholds (representing another valid choice) and nPCs shoots up after becoming slightly more stringent. 
# We take the point right before the number of PCs begins to shoot up
def optimize_dev_thr(adata, pca_name, save=False):
    
    plt.style.use('default')
    if pca_name not in adata.uns.keys():
        exit(f'ERROR: {pca_name} not in adata.')
    
    x = range(len(adata.uns[pca_name]['variance_ratio']))
    y = np.cumsum(adata.uns[pca_name]['variance_ratio'])
    
    end_line_slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
    end_line_int = y[-1] - end_line_slope*x[-1]
    
    deviation = []
    slopes = []
    dev_thr_list = [0.001, 0.0009, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003, 0.0002,
                    0.0001, 0.00009, 0.00008, 0.00007, 0.00006, 0.00005, 0.00004, 0.00003, 0.00002, 
                    0.00001, 0.000009, 0.000008, 0.000007, 0.000006, 0.000005, 0.000004, 0.000003, 0.000002,
                    0.00000001, 0.000000001, 0.0000000001]
    
    thr_hit=False
    
    opt_n_pcs_linear_dev_list=[]
    
    for dev_thr in dev_thr_list:
        
        print(dev_thr)
        thr_hit=False
        
        for i in range(len(x)-2):
            
            tmp_slope = (y[i+1] - y[i]) / (x[i+1] - x[i])
            tmp_int = y[i] - tmp_slope*x[i]
            
            plt.scatter(x, y, s=10, color='blue')
            
            y_line = [tmp_slope*x + tmp_int for x in x]
            
            tmp_dev = abs(y_line[i+2] - y[i+2])
            
            deviation.append(tmp_dev)
            
            slopes.append(tmp_slope)
            
            plt.plot(x, y_line, color='grey', linewidth=0.5)
            
            if (tmp_dev < dev_thr) and (thr_hit==False):
                print(f'Optimal nPCs = {i+1}')
                print(f'Optimal PCs index = {i}')
                print(f'Deviation = {tmp_dev}')
                
                opt_n_pcs_linear_dev = i
                opt_n_pcs_linear_dev_list.append(i)
                
                plt.scatter(x[i], y[i], s=15, color='red')
                plt.axvline(x[i], color='red', linestyle='--', linewidth=1)
                
                thr_hit=True
        
        if thr_hit==False:
            print(f'Optimal nPCs = {len(x)}')
            print(f'Optimal PCs index = {len(x)-1}')
            print(f'Deviation = {tmp_dev}')
            
            opt_n_pcs_linear_dev = len(x)-1
            opt_n_pcs_linear_dev_list.append(len(x)-1)
            
            plt.scatter(x[i], y[i], s=15, color='red')
            plt.axvline(x[i], color='red', linestyle='--', linewidth=1)
            
            thr_hit=True
            
        plt.xlabel('PCs Ranked')
        plt.ylabel('Cumulative Variance')
        plt.title(f'Dev Thr = {dev_thr}, i = {opt_n_pcs_linear_dev} | nPCs = {opt_n_pcs_linear_dev+1}')
        
        plt.ylim([y[0]*0.9, y[-1]*1.1])

        if save:
            fn =  + f'{output_dir}/PCA_linear_deviation_optimization_dev_thr_{dev_thr}.png'
            plt.tight_layout()
            print(fn)
            plt.savefig(fn, dpi=400)

        plt.show()
    
    fig, ax = plt.subplots(figsize=(6,4))
    
    plt.plot(dev_thr_list, opt_n_pcs_linear_dev_list)
    plt.scatter(dev_thr_list, opt_n_pcs_linear_dev_list, color='red', s=10)
    
    ax.invert_xaxis()
    plt.ylabel('Optimal PCs')
    plt.xlabel('Deviation Threshold')

    if save:
        fn = f'{output_dir}/PCA_linear_deviation_optimization_dev_thr_stability.png'
        plt.tight_layout()
        print(fn)
        plt.savefig(fn, dpi=400)
    
    fig, ax = plt.subplots(figsize=(6,4))
    
    plt.plot(dev_thr_list[9:], opt_n_pcs_linear_dev_list[9:])
    plt.scatter(dev_thr_list[9:], opt_n_pcs_linear_dev_list[9:], color='red', s=10)
    
    ax.invert_xaxis()
    plt.ylabel('Optimal PCs')
    plt.xlabel('Deviation Threshold')

    if save:
        fn = f'{output_dir}/PCA_linear_deviation_optimization_dev_thr_stability_2.png'
        plt.tight_layout()
        print(fn)
        plt.savefig(fn, dpi=400)


### MAST & GSEA HELPER FUNCTIONS ###
# Define functions for running MAST and GSEA
def my_run_mast(cond1, cond2, infile=None, data_field='logX', hvg_filter='TRUE', sample_covariate='FALSE', covar_field='sample', cell_filter='FALSE',  
                readonly=False, wait=True, force=False, plot=True, num_cores=30,
                use_filter=None, xlabel_cutoff=0.1, groupName='condition', reverse=False,
                sortBy='scaled_rank_score', quiet=True, mast_outdir=None, outname_str="", num_label=30):

    if use_filter is None:
        filterExt=''
        filterGenes=None
        filterName=''
    else:
        filterExt=use_filter['ext']
        filterGenes=use_filter['genes']
        filterName=' ( ' + use_filter['name'] + ')'
        
    if groupName == 'geneExpression':
        mast_outdir = f'{output_dir}/MAST-geneExp/{cond1}'.replace('cellType ', '_')
        outname = f'{cond1}_geneExp{cond2}_{outname_str}'.replace(' ', '_')
    else:
        if mast_outdir == None:
            mast_outdir = f'{output_dir}/{cond1}_vs_{cond2}/MAST'.replace(' ', '_')
            
    outname=f'{cond1}_vs_{cond2}_{outname_str}'.replace(' ', '_')
            
    stdoutfile=f"{mast_outdir}/stdout{outname_str}.txt"
    outfile = f"{mast_outdir}/{outname}.csv"
    plotfile = f"{mast_outdir}/{outname}{filterExt}.png".replace(' ', '_')
    
    if not readonly:
        if (infile is None):
            print("No infile (h5ad) given! Cannot run MAST.")
            return(None)
        if not path.exists(outfile):
            if not force:
                print(f"{outfile} not found, setting force to True")
            force=True
        os.system(f'{git_dir}/MAST_new.R /root/scripts/')
        cmd=f'mkdir -p "{mast_outdir}"'
        os.system(cmd)
        outname_arg=""
        
        if outname_str != "":
            outname_arg=f"--outname-str \"{outname_str}\""
            
        cmd=f'R --vanilla --args --infile "{infile}" --data_field "{data_field}" --outdir "{mast_outdir}" --hvg_filter "{hvg_filter}" \
                 --sample_covariate "{sample_covariate}" --covar_field {covar_field} --cell_filter "{cell_filter}" --groups {groupName} --comp-groups "{cond1}" "{cond2}"'
        cmd += f' --numcores {num_cores} {outname_arg} < /root/scripts/MAST_new.R'
        
        run_in_background(cmd, stdoutfile, force=force, wait=wait, quiet=quiet)
        
    if not path.exists(outfile):
        if (readonly):
            if not reverse:
                return(run_mast(cond2, cond1, infile=infile, readonly=True, plot=plot, use_filter=use_filter,
                               xlabel_cutoff=xlabel_cutoff, groupName=groupName, reverse=True, sortBy=sortBy))
            else:
                print("outfile " + outfile + " does not exist, returning")
        else:
            print(cmd + " has not finished, returning.")
        return(None)
    
    mastResults = read_mast_results(outfile, filterGenes, filterExt, reverse=reverse, sortBy=sortBy)
    
    if plot:
        plot_mast_results(mastResults, title='MAST results ' + ' ' + cond1 + ' vs ' + cond2 + filterName, 
                          num_label=num_label, arrows=True, fontsize=15, plot_outfile=plotfile, xlabel_cutoff=xlabel_cutoff)
    return(mastResults)

def read_mast_results(filename, filter_genes=None, filter_label='_filtered', reverse=False,
                      sortBy='scaled_rank_score'):
    print("Reading " + filename)
    mastResults = pd.read_csv(filename)
    mastResults.rename(index=str, columns={"primerid": "gene", "coef": "log2FC", 'Pr(>Chisq)':'p',
                                          'Pr..Chisq.':'p'}, inplace=True, errors='ignore')
    pmin= np.min(np.array([x for x in mastResults.p if x!=0]))
    fdrmin = np.min(np.array([x for x in mastResults.fdr if x!=0]))
    #mastResults.log2FC[np.isnan(mastResults.log2FC)] = np.nanmax(mastResults.log2FC)
    mastResults = mastResults[(np.isnan(mastResults['ci.lo']) == False) & (np.isnan(mastResults['ci.hi']) == False)]
    mastResults.drop(['Unnamed: 0'], axis=1, inplace=True)
    mastResults.set_index('gene', drop=True, inplace=True)
    mastResults['bonferroni'] = mastResults['p']*mastResults.shape[0]
    bonmin = np.min(np.array([x for x in mastResults.bonferroni if x!=0]))
    mastResults.loc[mastResults.p==0,'p'] = pmin
    mastResults.loc[mastResults.fdr==0,'fdr'] = fdrmin
    mastResults.loc[mastResults.bonferroni==0,'bonferroni'] = bonmin
    mastResults.loc[mastResults.bonferroni > 1,'bonferroni'] = 1
    mastResults['rank_score'] = -10*np.log10(mastResults['bonferroni'])*np.sign(mastResults['log2FC'])
    mastResults['FC'] = 2.0**mastResults['log2FC']
    mastResults['scaled_rank_score'] = mastResults['rank_score']*np.abs(mastResults['log2FC'])
    mastResults['abs_scaled_rank_score'] = np.abs(mastResults['scaled_rank_score'])
    if reverse:
        mastResults['log2FC'] = -mastResults['log2FC']
        mastResults['FC'] = 1.0/mastResults['FC']
        mastResults['ci.hi'] = -mastResults['ci.hi']
        mastResults['ci.lo'] = -mastResults['ci.lo']
        mastResults['rank_score'] = -mastResults['rank_score']
        mastResults['scaled_rank_score'] = -mastResults['scaled_rank_score']
    mastResults = mastResults.sort_values(by=sortBy, ascending=False)
    if filter_genes is not None:
        mastResults = mastResults.loc[[g for g in mastResults.index if g in filter_genes]]
        filename=os.path.splitext(filename)[0] + filter_label + '.csv'
        mastResults[['p', 'log2FC', 'FC', 'fdr', 'ci.hi','ci.lo','bonferroni','rank_score', 'scaled_rank_score']].to_csv(filename)
    return(mastResults)

def make_volcano_plot(df, title='MAST volcano plot', plot_outfile = None, 
                      ptype='bonferroni', mlog10_thresh=-np.log10(0.1), log2FC_thresh = np.log2(1.1),
                      num_label=15, arrows=True, 
                      fontsize=17,
                      figsize=(12,10),
                      axis_lines=False,
                      xlabel_cutoff=0.1,
                      minY=-1,
                      ycol='fdr',
                      xcol='log2FC',
                      xlabel='$log_2(FC)$',
                      ylabel='$-log{10}(p_{adj})$',
                      labelcol='index', 
                      label_list=None,
                      color_var=None,
                      cmap=None,
                      s=5, color='r'):
    
    print("ptype:", ptype)
    print("ycol:", ycol)
    # Identify significant genes to highlightlabelle
    target_num_label = num_label
    num_label=500
    while num_label > target_num_label:
        f = ((df[ycol] < 0.05) & ((df[xcol] > xlabel_cutoff) | (df[xcol] < -xlabel_cutoff)) & (df[ptype] < 0.5))
        num_label = np.sum(f)
        if (num_label > target_num_label):
            xlabel_cutoff = xlabel_cutoff+0.01
    print(f"num_label={num_label}, xlabel_cutoff={xlabel_cutoff}")

    x = df[xcol].to_numpy()
    
    y = df[ptype].to_numpy()
    if (np.sum(y==0) > 0):
        y[y==0] = np.min(y[y!=0])/2
    y = -np.log10(abs(df[ptype].to_numpy()))*np.sign(df[ptype]).to_numpy()
    
    if label_list==None:
        sig_ind = np.where((np.abs(x) > log2FC_thresh) * (np.abs(y) > mlog10_thresh))[0]    
        f = sig_ind
        
    else:
        if labelcol=='index':
            sig_ind = np.array(df.index.map(lambda x: x in label_list)).astype(bool)
        else:
            sig_ind = np.array(df[labelcol].map(lambda x: x in label_list)).astype(bool)
        f = np.where((np.abs(x) > log2FC_thresh) * (np.abs(y) > mlog10_thresh))[0]  
        
    fig, ax = plt.subplots(figsize=figsize)
    
    # Make Volcano plot
    plt.scatter(x, y, s=s, c='k')
    
    if color_var == None:
        plt.scatter(x[sig_ind], y[sig_ind], s=s, c=color)
    else:
        if cmap == None:
            print('No cmap provided.')
        else:
            plt.scatter(x[sig_ind], y[sig_ind], s=s, c=[cmap[x] for x in df.iloc[sig_ind][color_var].values])
        
    plt.title(title, fontsize=14)
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)
    plt.xlabel(xlabel, size=14, weight='normal')
    plt.ylabel(ylabel, size=14, weight='normal')
    #maxX = np.nanmax(np.append(np.abs(x), 1.5))
    #maxY = np.nanmax(np.append(y, 350))
    maxX = np.nanmax(np.abs(x))
    maxY = np.nanmax(y)
    plt.xlim(-maxX-0.75, maxX+0.75)
    plt.ylim(minY, maxY*1.33)
    #plt.grid(b=None)
    sns.despine()
    ax=plt.gca()
    ax.grid(False)
    
    if axis_lines:
        plt.vlines(x=0, ymin=minY, ymax=maxY+100, linestyles='dashed', color='grey')
        plt.hlines(y=0, xmin=-maxX-1, xmax=maxX+1, linestyles='dashed', color='grey')
    
    #f = (df['fdr'] < 0.01) & ((df['log2FC'] > 0.2) | (df['log2FC'] < -0.2))
    if labelcol == 'index':
        z = sorted(zip(y[f], x[f], df.index[f]), reverse=True)
    else:
        z = sorted(zip(y[f], x[f], df[labelcol][f]), reverse=True)
        
    
    if label_list==None:
        if num_label > 0:
    #        z=sorted(zip(y,x,df.index), reverse=True)[:num_label]
            texts = []
            for i in range(len(z)):
                texts.append(ax.text(z[i][1], z[i][0], z[i][2], fontsize=fontsize))
            if (arrows):
                niter=adjust_text(texts, x=x, y=y, 
                                  #force_text=(0.5,0.5),
                                  #force_points=(0.5,0.5),
                                  #expand_text=(2,2),
                                  #expand_points=(2,2),
                                  #autoalign='x',
                                  precision=0.001,
                                  arrowprops=dict(arrowstyle='-|>', color='gray', lw=0.5))
            else:
                niter = adjust_text(texts, x=x, y=y, force_text=0.05)
            #print("niter=" + str(niter))
        
    else:
        z = [x for x in z if x[2] in label_list]
            
        #print(z)
        
        texts = []
        for i in range(len(z)):
            texts.append(ax.text(z[i][1], z[i][0], z[i][2], fontsize=fontsize))
        if (arrows):
            niter=adjust_text(texts, x=x, y=y, 
                              #force_text=(0.5,0.5),
                              #force_points=(0.5,0.5),
                              #expand_text=(2,2),
                              #expand_points=(2,2),
                              #autoalign='x',
                              precision=0.001,
                              arrowprops=dict(arrowstyle='-|>', color='gray', lw=0.5))
        else:
            niter = adjust_text(texts, x=x, y=y, force_text=0.05)

    
    # SAVE FIGURE
    if plot_outfile is not None:
        d = os.path.dirname(plot_outfile)
        if not os.path.exists(d):
            os.makedirs(d)
        plt.savefig(plot_outfile, bbox_inches='tight', dpi=400)
        print("Wrote " + plot_outfile)
    plt.show()
    

# Import mouse to human gene mapping dictionary
# The mToH_mapping dictionary is hosted as a pkl file here: https://github.com/LaughneyLab/Lung_Histological_Transformation/tree/main/data/
with open('/workdir/varmus_single_cell/data/mToH_mapping.pkl', 'rb') as file:
    mToH_mapping = pickle.load(file)

# Define GSEA helper functions
def run_gsea(rank, output_root, fdr_cutoff=0.25, label='GSEA run', force=False, wait=False,
            gmtfile=f'', return_command_only=False,
            readonly=False, gene_map=None):
    label = label.replace(" ", "_")
    
    # look to see if output folder already exists
    if return_command_only and not force:
        f = glob.glob(output_root + "/**/*gsea_report*pos*tsv", recursive=True)
        if len(f):
            return None
    if not readonly:
        Path(output_root).mkdir(parents=True, exist_ok=True)
    rnkFile=f'{output_root}/input.rnk'
    if rank is not None:
        if gene_map is not None:
            rank = pd.DataFrame(rank)
            rank['hgene'] = [("" if i.upper() in gene_map.values() else i.upper()) if i not in gene_map else gene_map[i] for i in rank.index]            
            droprows = (rank['hgene'] == '')
            print(f'Dropping {sum(droprows)} results with no human gene mapping')
            rank = rank.loc[~droprows]
            rank.set_index('hgene', inplace=True)
        else:
            rank.index = [i.upper() for i in rank.index]
        rank.to_csv(rnkFile, sep='\t', header=False, index=True)
    cmd=f'{git_dir}/GSEA_4.1.0/gsea-cli.sh GSEAPreranked -rnk {rnkFile} -gmx {gmtfile} -collapse No_Collapse -mode Max_probe -norm meandiv'
    cmd += f' -nperm 10000 -scoring_scheme weighted -rpt_label {label} -create_svgs true -include_only_symbols true'
    cmd += f' -make_sets true -plot_top_x 20 -rnd_seed 888 -set_max 1500 -set_min 1 -zip_report false -out {output_root}'
    if return_command_only:
        return(cmd)
    stdoutfile=f'{output_root}/{label}_stdout.txt'
    if not readonly:
        run_in_background(cmd,f'{output_root}/{label}_stdout.txt',  wait=wait, force=force, quiet=True)
    # recover information from run
    f = glob.glob(output_root + "/**/*gsea_report*pos*tsv", recursive=True)
    #print(f, output_root + "/**/*gsea_report*pos*tsv")
    if f is None:
        raise RuntimeError(
            'seqc.JavaGSEA was not able to recover the output of the Java '
            'executable. This likely represents a bug.')
    f.sort(key=os.path.getmtime)
    f = f[0]
    f2 = f.replace('_pos_', '_neg_')
    if not path.exists(f2):
        raise RuntimeError("Error finding neg result file for GSEA " + f2)
    names = ['name', 'size', 'es', 'nes', 'p', 'fdr_q', 'fwer_p', 'rank_at_max', 'leading_edge']
    pos = pd.read_csv(f, sep='\t', infer_datetime_format=False, parse_dates=False).iloc[:, :-1]
    pos.drop(['GS<br> follow link to MSigDB', 'GS DETAILS'], axis=1, inplace=True)
    
    neg = pd.read_csv(f2, sep='\t', infer_datetime_format=False, parse_dates=False).iloc[:, :-1]
    neg.drop(['GS<br> follow link to MSigDB', 'GS DETAILS'], axis=1, inplace=True)
    pos.columns, neg.columns = names, names
    pos['direction'] = 'pos'
    neg['direction'] = 'neg'
    aa = pos.sort_values('fdr_q', ascending=True).fillna(0)
    bb = neg.sort_values('fdr_q', ascending=True).fillna(0)
    aa = aa[aa['fdr_q'] <= fdr_cutoff]
    bb = bb[bb['fdr_q'] <= fdr_cutoff]
    outcsv = f'{output_root}/{label}.csv'.replace(' ', '_')
    aa.to_csv(outcsv, sep='\t', index=False)
    with open(outcsv, 'a') as f:
        f.write('\n')
    bb.to_csv(outcsv, sep='\t', index=False, mode='a')
    return pd.concat([pos, neg])


def gsea_linear_scale(data: pd.Series) -> pd.Series:
        """scale input vector to interval [-1, 1] using a linear scaling
        :return correlations: pd.Series, data scaled to the interval [-1, 1]
        """
        data = data.copy()
        data -= np.min(data, axis=0)
        data /= np.max(data, axis=0) / 2
        data -= 1
        return data