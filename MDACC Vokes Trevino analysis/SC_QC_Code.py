#!/usr/bin/env python
# coding: utf-8

#import packages

import numpy as np
import pandas as pd
import scanpy as sc
import os


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
pd.options.display.max_columns = None



# # Read Data




#get folder names
folders = os.listdir(directory)
print(folders)
if '.DS_Store' in folders:
    folders.pop(folders.index('.DS_Store'))
print(f"\nThere are {len(folders)} samples.\n")





#Remove GEX from sample names
sampleindex = []
for f in folders:
    if 'GEX' in f:
        sampleindex.append(f.split('_')[0])
    else:
        sampleindex.append(f)

new = []
for s in sampleindex:
    if 'GEX' in s:
        hold =s.split('GEX')[0]
        if len(hold) < 3:
            new.append('JH0'+hold)
        else:
            new.append('JH'+hold)
    else:
        new.append(s)
sampleindex = new


print(sampleindex)


#read data
adata1 = sc.read_10x_h5(directory + '/'+ folders[0] + '/' +'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata1.var_names_make_unique()
adata2 = sc.read_10x_h5(directory + '/'+ folders[1] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata2.var_names_make_unique()
adata3 = sc.read_10x_h5(directory + '/'+ folders[2] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata3.var_names_make_unique()
adata4 = sc.read_10x_h5(directory + '/'+ folders[3] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata4.var_names_make_unique()
adata5 = sc.read_10x_h5(directory + '/'+ folders[4] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata5.var_names_make_unique()
adata6 = sc.read_10x_h5(directory + '/'+ folders[5] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata6.var_names_make_unique()
adata7 = sc.read_10x_h5(directory + '/'+ folders[6] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata7.var_names_make_unique()
adata8 = sc.read_10x_h5(directory + '/'+ folders[7] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata8.var_names_make_unique()
adata9 = sc.read_10x_h5(directory + '/'+ folders[8] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata9.var_names_make_unique()
adata10 = sc.read_10x_h5(directory + '/'+ folders[9] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata10.var_names_make_unique()
adata11 = sc.read_10x_h5(directory + '/'+ folders[10] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata11.var_names_make_unique()
adata12 = sc.read_10x_h5(directory + '/'+ folders[11] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata12.var_names_make_unique()
adata13 = sc.read_10x_h5(directory + '/'+ folders[12] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata13.var_names_make_unique()
adata14 = sc.read_10x_h5(directory + '/'+ folders[13] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata14.var_names_make_unique()
adata15 = sc.read_10x_h5(directory + '/'+ folders[14] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata15.var_names_make_unique()
adata16 = sc.read_10x_h5(directory + '/'+ folders[15] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata16.var_names_make_unique()
adata17 = sc.read_10x_h5(directory + '/'+ folders[16] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata17.var_names_make_unique()
adata18 = sc.read_10x_h5(directory + '/'+ folders[17] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata18.var_names_make_unique()
adata19 = sc.read_10x_h5(directory + '/'+ folders[18] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata19.var_names_make_unique()
adata20 = sc.read_10x_h5(directory + '/'+ folders[19] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata20.var_names_make_unique()
adata21 = sc.read_10x_h5(directory + '/'+ folders[20] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata21.var_names_make_unique()
adata22 = sc.read_10x_h5(directory + '/'+ folders[21] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata22.var_names_make_unique()
adata23 = sc.read_10x_h5(directory + '/'+ folders[22] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata23.var_names_make_unique()
adata24 = sc.read_10x_h5(directory + '/'+ folders[23] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata24.var_names_make_unique()
adata25 = sc.read_10x_h5(directory + '/'+ folders[24] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata25.var_names_make_unique()
adata26 = sc.read_10x_h5(directory + '/'+ folders[25] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata26.var_names_make_unique()
adata27 = sc.read_10x_h5(directory + '/'+ folders[26] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata27.var_names_make_unique()
adata28 = sc.read_10x_h5(directory + '/'+ folders[27] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata28.var_names_make_unique()
adata29 = sc.read_10x_h5(directory + '/'+ folders[28] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata29.var_names_make_unique()
adata30 = sc.read_10x_h5(directory + '/'+ folders[29] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata30.var_names_make_unique()
adata31 = sc.read_10x_h5(directory + '/'+ folders[30] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata31.var_names_make_unique()
adata32 = sc.read_10x_h5(directory + '/'+ folders[31] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata32.var_names_make_unique()
adata33 = sc.read_10x_h5(directory + '/'+ folders[32] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata33.var_names_make_unique()
adata34 = sc.read_10x_h5(directory + '/'+ folders[33] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata34.var_names_make_unique()
adata35 = sc.read_10x_h5(directory + '/'+ folders[34] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata35.var_names_make_unique()
adata36 = sc.read_10x_h5(directory + '/'+ folders[35] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata36.var_names_make_unique()
adata37 = sc.read_10x_h5(directory + '/'+ folders[36] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata37.var_names_make_unique()
adata38 = sc.read_10x_h5(directory + '/'+ folders[37] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata38.var_names_make_unique()
adata39 = sc.read_10x_h5(directory + '/'+ folders[38] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata39.var_names_make_unique()
adata40 = sc.read_10x_h5(directory + '/'+ folders[39] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata40.var_names_make_unique()
adata41 = sc.read_10x_h5(directory + '/'+ folders[40] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata41.var_names_make_unique()
adata42 = sc.read_10x_h5(directory + '/'+ folders[41] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata42.var_names_make_unique()
adata43 = sc.read_10x_h5(directory + '/'+ folders[42] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata43.var_names_make_unique()
adata44 = sc.read_10x_h5(directory + '/'+ folders[43] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata44.var_names_make_unique()
adata45 = sc.read_10x_h5(directory + '/'+ folders[44] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata45.var_names_make_unique()
adata46 = sc.read_10x_h5(directory + '/'+ folders[45] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata46.var_names_make_unique()
adata47 = sc.read_10x_h5(directory + '/'+ folders[46] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata47.var_names_make_unique()
adata48 = sc.read_10x_h5(directory + '/'+ folders[47] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata48.var_names_make_unique()
adata49 = sc.read_10x_h5(directory + '/'+ folders[48] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata49.var_names_make_unique()
adata50 = sc.read_10x_h5(directory + '/'+ folders[49] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata50.var_names_make_unique()
adata51 = sc.read_10x_h5(directory + '/'+ folders[50] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata51.var_names_make_unique()
adata52 = sc.read_10x_h5(directory + '/'+ folders[51] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata52.var_names_make_unique()
adata53 = sc.read_10x_h5(directory + '/'+ folders[52] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata53.var_names_make_unique()
adata54 = sc.read_10x_h5(directory + '/'+ folders[53] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata54.var_names_make_unique()
adata55 = sc.read_10x_h5(directory + '/'+ folders[54] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata55.var_names_make_unique()
adata56 = sc.read_10x_h5(directory + '/'+ folders[55] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata56.var_names_make_unique()
adata57 = sc.read_10x_h5(directory + '/'+ folders[56] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata57.var_names_make_unique()
adata58 = sc.read_10x_h5(directory + '/'+ folders[57] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata58.var_names_make_unique()
adata59 = sc.read_10x_h5(directory + '/'+ folders[58] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata59.var_names_make_unique()
adata60 = sc.read_10x_h5(directory + '/'+ folders[59] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata60.var_names_make_unique()
adata61 = sc.read_10x_h5(directory + '/'+ folders[60] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata61.var_names_make_unique()
adata62 = sc.read_10x_h5(directory + '/'+ folders[61] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata62.var_names_make_unique()
adata63 = sc.read_10x_h5(directory + '/'+ folders[62] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata63.var_names_make_unique()
adata64 = sc.read_10x_h5(directory + '/'+ folders[63] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata64.var_names_make_unique()
adata65 = sc.read_10x_h5(directory + '/'+ folders[64] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata65.var_names_make_unique()
adata66 = sc.read_10x_h5(directory + '/'+ folders[65] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata66.var_names_make_unique()
adata67 = sc.read_10x_h5(directory + '/'+ folders[66] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata67.var_names_make_unique()
adata68 = sc.read_10x_h5(directory + '/'+ folders[67] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata68.var_names_make_unique()
adata69 = sc.read_10x_h5(directory + '/'+ folders[68] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata69.var_names_make_unique()
adata70 = sc.read_10x_h5(directory + '/'+ folders[69] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata70.var_names_make_unique()
adata71 = sc.read_10x_h5(directory + '/'+ folders[70] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata71.var_names_make_unique()
adata72 = sc.read_10x_h5(directory + '/'+ folders[71] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata72.var_names_make_unique()
adata73 = sc.read_10x_h5(directory + '/'+ folders[72] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata73.var_names_make_unique()
adata74 = sc.read_10x_h5(directory + '/'+ folders[73] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata74.var_names_make_unique()
adata75 = sc.read_10x_h5(directory + '/'+ folders[74] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata75.var_names_make_unique()
adata76 = sc.read_10x_h5(directory + '/'+ folders[75] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata76.var_names_make_unique()
adata77 = sc.read_10x_h5(directory + '/'+ folders[76] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata77.var_names_make_unique()
adata78 = sc.read_10x_h5(directory + '/'+ folders[77] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata78.var_names_make_unique()
adata79 = sc.read_10x_h5(directory + '/'+ folders[78] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata79.var_names_make_unique()
adata80 = sc.read_10x_h5(directory + '/'+ folders[79] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata80.var_names_make_unique()
adata81 = sc.read_10x_h5(directory + '/'+ folders[80] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata81.var_names_make_unique()
adata82 = sc.read_10x_h5(directory + '/'+ folders[81] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata82.var_names_make_unique()
adata83 = sc.read_10x_h5(directory + '/'+ folders[82] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata83.var_names_make_unique()
adata84 = sc.read_10x_h5(directory + '/'+ folders[83] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata84.var_names_make_unique()
adata85 = sc.read_10x_h5(directory + '/'+ folders[84] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata85.var_names_make_unique()
adata86 = sc.read_10x_h5(directory + '/'+ folders[85] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata86.var_names_make_unique()
adata87 = sc.read_10x_h5(directory + '/'+ folders[86] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata87.var_names_make_unique()
adata88 = sc.read_10x_h5(directory + '/'+ folders[87] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata88.var_names_make_unique()
adata89 = sc.read_10x_h5(directory + '/'+ folders[88] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata89.var_names_make_unique()
adata90 = sc.read_10x_h5(directory + '/'+ folders[89] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata90.var_names_make_unique()
adata91 = sc.read_10x_h5(directory + '/'+ folders[90] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata91.var_names_make_unique()
adata92 = sc.read_10x_h5(directory + '/'+ folders[91] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata92.var_names_make_unique()
adata93 = sc.read_10x_h5(directory + '/'+ folders[92] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata93.var_names_make_unique()
adata94 = sc.read_10x_h5(directory + '/'+ folders[93] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata94.var_names_make_unique()
adata95 = sc.read_10x_h5(directory + '/'+ folders[94] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata95.var_names_make_unique()
adata96 = sc.read_10x_h5(directory + '/'+ folders[95] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata96.var_names_make_unique()
adata97 = sc.read_10x_h5(directory + '/'+ folders[96] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata97.var_names_make_unique()
adata98 = sc.read_10x_h5(directory + '/'+ folders[97] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata98.var_names_make_unique()
adata99 = sc.read_10x_h5(directory + '/'+ folders[98] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata99.var_names_make_unique()
adata100 = sc.read_10x_h5(directory + '/'+ folders[99] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata100.var_names_make_unique()
adata101 = sc.read_10x_h5(directory + '/'+ folders[100] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata101.var_names_make_unique()
adata102 = sc.read_10x_h5(directory + '/'+ folders[101] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata102.var_names_make_unique()
adata103 = sc.read_10x_h5(directory + '/'+ folders[102] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata103.var_names_make_unique()
adata104 = sc.read_10x_h5(directory + '/'+ folders[103] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata104.var_names_make_unique()
adata105 = sc.read_10x_h5(directory + '/'+ folders[104] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata105.var_names_make_unique()
adata106 = sc.read_10x_h5(directory + '/'+ folders[105] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata106.var_names_make_unique()
adata107 = sc.read_10x_h5(directory + '/'+ folders[106] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata107.var_names_make_unique()
adata108 = sc.read_10x_h5(directory + '/'+ folders[107] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata108.var_names_make_unique()
adata109 = sc.read_10x_h5(directory + '/'+ folders[108] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata109.var_names_make_unique()
adata110 = sc.read_10x_h5(directory + '/'+ folders[109] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata110.var_names_make_unique()
adata111 = sc.read_10x_h5(directory + '/'+ folders[110] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata111.var_names_make_unique()
adata112 = sc.read_10x_h5(directory + '/'+ folders[111] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata112.var_names_make_unique()
adata113 = sc.read_10x_h5(directory + '/'+ folders[112] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata113.var_names_make_unique()
adata114 = sc.read_10x_h5(directory + '/'+ folders[113] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata114.var_names_make_unique()
adata115 = sc.read_10x_h5(directory + '/'+ folders[114] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata115.var_names_make_unique()
adata116 = sc.read_10x_h5(directory + '/'+ folders[115] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata116.var_names_make_unique()
adata117 = sc.read_10x_h5(directory + '/'+ folders[116] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata117.var_names_make_unique()
adata118 = sc.read_10x_h5(directory + '/'+ folders[117] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata118.var_names_make_unique()
adata119 = sc.read_10x_h5(directory + '/'+ folders[118] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata119.var_names_make_unique()
adata120 = sc.read_10x_h5(directory + '/'+ folders[119] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata120.var_names_make_unique()
adata121 = sc.read_10x_h5(directory + '/'+ folders[120] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata121.var_names_make_unique()
adata122 = sc.read_10x_h5(directory + '/'+ folders[121] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata122.var_names_make_unique()
adata123 = sc.read_10x_h5(directory + '/'+ folders[122] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata123.var_names_make_unique()
adata124 = sc.read_10x_h5(directory + '/'+ folders[123] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata124.var_names_make_unique()
adata125 = sc.read_10x_h5(directory + '/'+ folders[124] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata125.var_names_make_unique()
adata126 = sc.read_10x_h5(directory + '/'+ folders[125] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata126.var_names_make_unique()
adata127 = sc.read_10x_h5(directory + '/'+ folders[126] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata127.var_names_make_unique()
adata128 = sc.read_10x_h5(directory + '/'+ folders[127] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata128.var_names_make_unique()
adata129 = sc.read_10x_h5(directory + '/'+ folders[128] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata129.var_names_make_unique()
adata130 = sc.read_10x_h5(directory + '/'+ folders[129] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata130.var_names_make_unique()
adata131 = sc.read_10x_h5(directory + '/'+ folders[130] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata131.var_names_make_unique()
adata132 = sc.read_10x_h5(directory + '/'+ folders[131] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata132.var_names_make_unique()
adata133 = sc.read_10x_h5(directory + '/'+ folders[132] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata133.var_names_make_unique()
adata134 = sc.read_10x_h5(directory + '/'+ folders[133] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata134.var_names_make_unique()
adata135 = sc.read_10x_h5(directory + '/'+ folders[134] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata135.var_names_make_unique()
adata136 = sc.read_10x_h5(directory + '/'+ folders[135] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata136.var_names_make_unique()
adata137 = sc.read_10x_h5(directory + '/'+ folders[136] + '/' + 'cellbender_10x_pbmc_filtered.h5')#,genome='GRCh38')
adata137.var_names_make_unique()


#combine into single adata
adata = adata1.concatenate(adata2,adata3,adata4,adata5,adata6,adata7,adata8,adata9,adata10,adata11,adata12,
                           adata13,adata14,adata15,adata16,adata17,adata18,adata19,adata20,adata21,adata22,
                           adata23,adata24,adata25,adata26,adata27,adata28,adata29,adata30,adata31,adata32,
                           adata33,adata34,adata35,adata36,adata37,adata38,adata39,adata40,adata41,adata42,
                           adata43,adata44,adata45,adata46,adata47,adata48,adata49,adata50,adata51,adata52,
                           adata53,adata54,adata55,adata56,adata57,adata58,adata59,adata60,adata61,adata62,
                           adata63,adata64,adata65,adata66,adata67,adata68,adata69,adata70,adata71,adata72,
                           adata73,adata74,adata75,adata76,adata77,adata78,adata79,adata80,adata81,adata82,
                           adata83,adata84,adata85,adata86,adata87,adata88,adata89,adata90,adata91,adata92,
                           adata93,adata94,adata95,adata96,adata97,adata98,adata99,adata100,adata101,adata102,
                           adata103,adata104,adata105,adata106,adata107,adata108,adata109,adata110,adata111,adata112,
                           adata113,adata114,adata115,adata116,adata117,adata118,adata119,adata120,adata121,adata122,
                           adata123,adata124,adata125,adata126,adata127,adata128,adata129,adata130,adata131,adata132,
                           adata133,adata134,adata135,adata136,adata137,
                           batch_key = 'batch',batch_categories=sampleindex)


#create dictionary of sample names and assign observation name
BATCHsubannot = {}
for i in range(0,len(folders)):
    #BATCHsubannot[str(i)] = folders[i]
    BATCHsubannot[sampleindex[i]] = folders[i]
print(BATCHsubannot)

# add a new `.obs` column called `BatchLegend` by mapping clusters to annotation using pandas `map` function
adata.obs['BatchLegend'] = adata.obs['batch'].map(BATCHsubannot).astype('category')


# # Initial Quality Control 


#filter out cells with less than 200 genes in each data set
sc.pp.filter_cells(adata,min_genes=200)


# # Add Meta Data


Meta_df = pd.read_excel(metadirectory + 'Sequenced Samples Sample Database.xlsx')



sample_names = Meta_df['specimen_id'].astype(str).apply(str.lower).tolist()



for i in range(0,len(sample_names)):
    if sample_names[i][0] == 'j':
        sample_names[i] = sample_names[i] + '_gex'




#sample_names[sample_names.index('lung -normal 10')] = 'lung-normal-10'
sample_names[sample_names.index('biopsy-4-liver')] = 'biopsy4'
#sample_names[sample_names.index('lung - tumor 9')] = 'lung-tumor-10'
sample_names[sample_names.index('biopsy-1')] = 'biopsy1'
sample_names[sample_names.index('jh105_atac_gex')] = 'jh105_atac'
sample_names[sample_names.index('jh067_atac_gex')] = 'jh067_atac'



for i in range(0,len(sample_names)):
    if 'jh-' in sample_names[i]:
        sample_names[i] = sample_names[i].split('-')[0] + sample_names[i].split('-')[1]



sample_names[sample_names.index('jh208_gex')] = 'jh208_exon14_gex'
sample_names[sample_names.index('jh205_gex')] = 'jh205_exon14_gex'
sample_names[sample_names.index('jh235_gex')] = 'jh235_exon14_gex'





#check sample names match names in folders
nfound = 0
names = []
for folder in folders:
    if folder.lower() not in sample_names:
        print(f"{folder.lower()} NOT FOUND")
        names.append(folder)
    else:
        #print(folder)
        
        nfound = nfound + 1
print(nfound)
print(names.sort())



#add patient number
patient_no = Meta_df['PATIENT_NO'].tolist()
annot = {}
for folder in folders:
    if folder.lower() in sample_names:
         annot[folder]=patient_no[sample_names.index(folder.lower())]

adata.obs['PATIENT_NO'] = adata.obs['BatchLegend'].map(annot).astype('category')



#add driver
#Driver_classification_condensed	
driver = Meta_df['Driver_classification_condensed'].tolist()

annot = {}
for folder in folders:
    if folder.lower() in sample_names:
         annot[folder]=driver[sample_names.index(folder.lower())]

adata.obs['DRIVER'] = adata.obs['BatchLegend'].map(annot).astype('category')




#add sequencing batch

Meta_df['sequencing_batch'].value_counts(dropna=False)
seqbatch = Meta_df['sequencing_batch'].tolist()

annot = {}
for folder in folders:
    if folder.lower() in sample_names:
         annot[folder]=seqbatch[sample_names.index(folder.lower())]

#print(annot)            
    
annot['JH033_GEX'] = 'Le_Vokes_Batch_3'
adata.obs['SEQBATCH'] = adata.obs['BatchLegend'].map(annot).astype('category')





#combine 2020 batch 
# Add SEQBATCH condensed
for key in annot.keys():
    if '2020' in str(annot[key]):
        annot[key] = 'EGFR_2020'
adata.obs['SEQBATCHcond'] = adata.obs['BatchLegend'].map(annot).astype('category')





#add 2020 batch vs New Batches
for key in annot.keys():
    if '2020' in str(annot[key]):
        annot[key] = 'EGFR_2020'
    else:
        annot[key] = 'Le_Vokes'
adata.obs['NVOBATCH'] = adata.obs['BatchLegend'].map(annot).astype('category')





#treatment condensed
treatment = Meta_df['Treatment_condensed'].tolist()

annot = {}
for folder in folders:
    if folder.lower() in sample_names:
         annot[folder]=treatment[sample_names.index(folder.lower())]

adata.obs['TREATMENT'] = adata.obs['BatchLegend'].map(annot).astype('category')


# # Continue Quality Control 




#determine % mitochondrial 
adata.var['mt'] = adata.var_names.str.startswith('MT-') #annotate the group 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p = False, inplace = True)





#determine % ribosomal
adata.var['ribo'] = adata.var_names.str.startswith(("RPS", "RPL")) #annotate the group 
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p = False, inplace = True)




#determine % mitochondrial 
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))
sc.pp.calculate_qc_metrics(adata, qc_vars=['hb'], percent_top=None, log1p = False, inplace = True)


#filter cells nUMI < 1000 
sc.pp.filter_cells(adata,min_counts=1000)




#filter cells with less than 500 genes
adata = adata[adata.obs.n_genes_by_counts > 500,:]




#filter cells with mitochondrial % > 20%
adata = adata[adata.obs.pct_counts_mt < 20,:]




#filter genes found in less than 10 cells
sc.pp.filter_genes(adata, min_cells = 10)




#filter cell with a complexity < 0.8
adata = adata[np.log10(adata.obs.n_genes_by_counts)/np.log10(adata.obs.total_counts) > 0.8,:]


# # Predict and remove Doublets w/Scrublet



import scrublet as scr
scrub = scr.Scrublet(adata.X)
adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
#adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.call_doublets(threshold=0.8)
#scrub.plot_histogram()

sum(adata.obs['predicted_doublets'])



scrub.plot_histogram()




# add in column with singlet/doublet instead of True/False

adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)

#filter doublets

adata = adata[adata.obs['doublet_info'] == 'False',:]


# # Normalization



adata.layers["counts"] = adata.X

#Normalize gene counts
sc.pp.normalize_total(adata, target_sum=1e6, exclude_highly_expressed=True, max_fraction = 0.05)

# Logarithmize the data.
sc.pp.log1p(adata,base=2)

#Store raw data, for finding markers in each cluster and other analysis
adata.raw = adata




#save raw counts for INFERCNV output
adata.to_df(layer="counts")



#store QC data
adata.write('/2024_02_23_adata_afterQC.h5ad')

