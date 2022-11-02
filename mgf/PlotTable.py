#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py

Functions:
    # plot_table():      Plot a table with results of the gap-filled sums

"""

#%% Initialization
### Imports from python
import matplotlib.pyplot as plt

### Own imports
import mgf

#%%
def plot_table(sums, ini, file_path='', run_number='', suffix=''):
    # plot_table():     Plot a table with sums and errors
    # sums:             Results dataframe with sums
    # ini:              Settings from ini file
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run
    # suffix:           Suffix to filename
#%%   
    numCols = sums.shape[1]
    fig_width = numCols*1
    fig_height = 5.0
    fsize = 12
    fig, plx = plt.subplots(figsize=(fig_width, fig_height))
    
    # Hide plot and axes
    fig.patch.set_visible(False)
    plx.axis('off')

    #vDraw table
    table = plx.table(cellText=sums.values, colLabels=sums.columns, rowLabels=sums.index, colColours=['lightgray']*numCols, loc='center')
    w, h = table[0,1].get_width(), table[0,1].get_height()
    table.add_cell(0, -1, w,h, text=sums.index.name, facecolor='lightgray')
    table.scale(1.3, 1.8)
 
    # Code to adjust height and passing
    # cellDict = table.get_celld()
    # for i in range(-1,sums.shape[1]):
    #     cellDict[(0,i)].set_height(0.3)
    #     for j in range(1,sums.shape[0]+1):
    #         cellDict[(j,i)].set_height(0.2)
    #         cellDict[(j,i)].PAD = 0.2

    mgf.utl.save_plot('sums', ini, file_path, run_number, suffix=suffix, folder='_res')
    
            