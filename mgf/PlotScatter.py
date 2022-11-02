#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py

Functions:
    # plot_scatter():   Scatter plot of data or residuals
    

    
Notes:
    For statistics and plots like R2, the model (MR) has to be on x-axis.
    Code inspired from Joe Kington:
        https://stackoverflow.com/questions/7941207/is-there-a-function-to-make-scatterplot-matrices-in-matplotlib

"""

#%% Initialization
### Imports from python
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats #linear regression

### Own imports
import mgf

"""
For TESTING:
    scenario='hhs';
"""

#%%
def plot_scatter(data, ini, scenario, file_path='', run_number='', suffix=''):
    # plot_scatter():   Scatter plot of data or residuals (residual scatter plot if 'res' in variable suffix)
    # data:             Dataframe with timestamp index, flux measurements, (maybe meteo,) and filled fluxes for various techniques
    # ini:              Settings from ini file
    # scenario:         Name of artificial gap filling scenario
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run
    # suffix:           Suffix to filename
    
#%%
    # Set data for plot
    flux_org = mgf.utl.name_flux(ini)
    col_ToPlot = data.filter(regex=scenario).columns
    data_sub = data[np.append(flux_org,col_ToPlot.values)][data[flux_org].notnull()]#[1:10]
    if 'res' in suffix:
        # Substract first column to get residuals
        data_res = data_sub.iloc[:,1:].subtract(data_sub[flux_org], axis=0)
        data_res[flux_org] = data_sub[flux_org]
        data_sub = data_res
    plt_min = data_sub.min().min()
    plt_max = data_sub.max().max()
    if 'res' in suffix: #even axis limits
        plt_min = -abs(max(plt_min, plt_max, key=abs))
        plt_max = abs(max(plt_min, plt_max, key=abs))
        
    # Figure settings
    numCols = 5 #max 3 for portrait, max 5 for landscape
    if (col_ToPlot.shape[0] <= 12): numCols = 3
    numRows = int(np.ceil(col_ToPlot.shape[0] / numCols))
    fig_square = 2.5 #Each box, 3.8, 5.0 ???
    fig_m_size = 3 #Marker size 4, 5
    fig_width = numCols*fig_square
    fig_height = numRows*fig_square
#
#    fig, plx = plt.subplots(figsize=(fig_width, fig_height))
#    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
#    fig.subplots_adjust(hspace=0.01)

    # Draw scatter plot
    fig, sub = plt.subplots(nrows=numRows, ncols=numCols, figsize=(fig_width, fig_height))
    fig.subplots_adjust(hspace=0.03, wspace=0.03)
    fig.patch.set_facecolor('white')
    #plt.rc('ytick',labelsize=12) #! Screws up labels of all plots ... !
    #plt.rc('xtick',labelsize=12)

    # Labels (y_0, ... need to be up here)
    fsize=12
    if 'res' in suffix: 
        y_0,y_1 = 0, 0 #for red line
        x_label = str('Observed fluxes for '+ini.FluxGas+' ('+ini.FluxUnit+')')
        y_label = str('Residuals of predicted fluxes for '+ini.FluxGas+' ('+ini.FluxUnit+')')
    else:
        y_0,y_1 = plt_min,plt_max #for red line
        x_label = str('Predicted fluxes for '+ini.FluxGas+' ('+ini.FluxUnit+')')
        y_label = str('Observed fluxes for '+ini.FluxGas+' ('+ini.FluxUnit+')')
    #cols_mid = int(numCols/2)#Middle of columns --> Not needed anymore with fig.supxlabel
    #sub[numRows-1,cols_mid].set_xlabel(x_label, fontsize=fsize, weight='bold')
    x_label_y = 0.03 #default y=0.01
    y_label_x = 0.02 #default x=0.02
    if (numCols==5): x_label_y = 0.05; y_label_x = 0.07
    fig.supxlabel(x_label, fontsize=fsize, weight='bold',y=x_label_y)
    fig.supylabel(y_label, fontsize=fsize, weight='bold',x=y_label_x)
    fig.suptitle('Scatterplots of artficial \''+scenario+'\' scenarios', x=0.5, y=0.95, fontsize=fsize+4)
    
    #Plot data_sub, line, label
    for x_i in range(0,numRows):
        for y_i in range(0,numCols):
            if (x_i*numCols+y_i) < col_ToPlot.shape[0]:
                sub[x_i,y_i].plot([plt_min,plt_max],[y_0,y_1], color='red', linestyle='solid', lw=0.5)
                sub[x_i,y_i].set_ylim(plt_min, plt_max) #needed since for residuals not limit of red line
                if 'res' in suffix:
                    y = data_sub[col_ToPlot[x_i*numCols+y_i]]
                    x = data_sub[flux_org]  #Observed on x-axis

                else:
                    x = data_sub[col_ToPlot[x_i*numCols+y_i]] #Model on x-axis (MR)
                    y = data_sub[flux_org]
                t_b = x.notnull() & y.notnull()
                sub[x_i,y_i].plot(x, y, linestyle='none', marker='o', markersize=fig_m_size, markerfacecolor='0.3', markeredgecolor='None', alpha=.5)
                slope, intercept, r_value, p_value, std_err = stats.linregress(x[t_b],y[t_b])
                sub[x_i,y_i].annotate(col_ToPlot[x_i*numCols+y_i], (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
                #{:04.2f}'.format(
                sub[x_i,y_i].annotate('R2='+'{:04.2f}'.format(r_value**2), (0.05, 0.9), xycoords='axes fraction', ha='left', fontsize=fsize-2)#, fontweight='bold')
                sub[x_i,y_i].tick_params(labelsize=fsize-2)
            else:
                sub[x_i,y_i].plot([plt_min,plt_max],[y_0,y_1], color='red', linestyle='solid', lw=0.5)
                sub[x_i,y_i].set_ylim(plt_min, plt_max) #needed since for residuals not limit of red line

    #Set ticks and axis properties
    for ax in sub.flat:
        # Hide all ticks and labels
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
                              
        # Set up ticks only on one side for the "edge" subplots...
        if ax.get_subplotspec().is_first_col():
            ax.yaxis.set_ticks_position('left')
        if ax.get_subplotspec().is_last_col():
            ax.yaxis.set_ticks_position('right')
        if ax.get_subplotspec().is_first_row():
            ax.xaxis.set_ticks_position('top')
        if ax.get_subplotspec().is_last_row():
            ax.xaxis.set_ticks_position('bottom')
            
    # Turn on the proper x or y axes ticks.
    for y_i in range(0,numCols):
        sub[0,y_i].xaxis.set_visible(True)
        sub[numRows-1,y_i].xaxis.set_visible(True)
    for x_i in range(0,numRows):
        sub[x_i,0].yaxis.set_visible(True)
        sub[x_i,numCols-1].yaxis.set_visible(True)

    
    if suffix != '': suffix = '_' + suffix 
    mgf.utl.save_plot('psc', ini, file_path, run_number, suffix=scenario+suffix)
    