#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import cellrank as cr
import scvelo as scv

# Plotting
import matplotlib.pyplot as plt


#read in subset
adata = sc.read('EGFRSub.h5ad')

#run cytotrace

adata.X = adata.raw.X
adata.layers["spliced"] = adata.X
adata.layers["unspliced"] = adata.X
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

from cellrank.kernels import CytoTRACEKernel

ctk = CytoTRACEKernel(adata).compute_cytotrace()

#save adata w/ cyto results
adata.write('EGFRSub_wCyto.h5ad')
