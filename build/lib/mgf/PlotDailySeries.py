#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py
    
Functions:
    # plot_series():    Plot time series of fluxes, residuals, three types of sums (data, gaps, all)

"""
#%% Initialization

### Imports from python
import matplotlib.pyplot as plt
import matplotlib.dates as pmd
import numpy as np

### Own imports
import mgf

"""  
For TESTING:
    scenario='real'; suffix = 'daily'
"""
#%%
def plot_daily(data, ini, scenario='real', file_path='', run_number='', suffix=''):
    # plot_daily():     Plot daily time series of fluxes, residuals, three types of sums (data, gaps, all)
    # data:             Dataframe with timestamp index, flux measurements, (maybe meteo,) and filled fluxes for various techniques
    # ini:              Settings from ini file
    # scenario:         Artificial gap filling scenario
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run
    # suffix:           Suffix to filename
#%%
    # Init stuff for plot
    flux_org = mgf.utl.name_flux(ini)
    factorHH = float(ini.ConvFactor)

    # Set data for plot
    col_ToPlot = data.filter(regex=scenario).columns
    data_sub = data[np.append(flux_org,col_ToPlot.values)]#[2000:3000]
    data_1d = data_sub.resample('D',closed='right',label='left').sum() * factorHH #Needs to be left for daily... 

    # # These subsets of data makes sense for 'hhs' and give more information than only 'real' --> not used any more but maybe helpful later
    #     # Points during measurements #original flux data for filled real gaps, artificial gaps for 'hhs' scenario
    #     data_meas_sub = data_sub[data_sub[flux_org].notnull()]
    #     data_1d_meas = data_meas_sub.resample('D',closed='right',label='left').sum() * factorHH 
    #     # Points during gaps #exist for 'hhs' and 'real'
    #     data_gaps_sub = data_sub[data_sub[flux_org].isnull()]
    #     data_1d_gaps = data_gaps_sub.resample('D',closed='right',label='left').sum() * factorHH

    # Plot properties
    numCols = 1
    numRows = 2
    fig_width = 15 #12.0
    fig_height = fig_width / 2 * numRows #20.0
    if ('ms' in suffix): fig_width = 12; fig_height = 9
    fig, sub = plt.subplots(nrows=numRows, ncols=numCols, figsize=(fig_width, fig_height))
    fig.subplots_adjust(hspace=0.03)
    fig.patch.set_facecolor('white')
    fsize=12

    # Set colors
    labels = data_1d.columns.str.replace('_'+scenario,'') #Flux column (labels[0]) gets explicit assignment of basic blue in plot
    lineColors, techColumns, usedKeys = mgf.utl.set_colors('gft_rainbow',labels[1:])

    # Time series daily sums
    sub[0].set_prop_cycle('color', lineColors)
    sub[0].plot(data_1d.iloc[:,1:], '+', markersize=4)
    #sub[0].plot(data_1d[flux_org], 'b_', markersize=10)
    sub[0].bar(data_1d.index, data_1d[flux_org], width=.5, align='center', fill=False, edgecolor='blue', linewidth=1)
    if ('ms' not in suffix): #If not for ms, annotate on plots
        sub[0].annotate('Daily sums of observed fluxes (blue bars) plus filled real gaps', 
                        (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
    sub[0].tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True, bottom=False, top=True, left=True, right=True)
 
    # Cumulative sum of sums
    sub[1].set_prop_cycle('color', lineColors)
    sub[1].plot(data_1d.iloc[:,1:].cumsum(skipna=True), '-', linewidth=0.5) # linestyle='solid'
    flux=data_1d[flux_org] #Just for better naming in legend
    flux.name = 'F_' + ini.FluxGas
    sub[1].plot(data_1d[flux_org].cumsum(skipna=True).ffill(), 'b-', linewidth=2) #ffill() seems to be only needed in older versions of py
 
    #sub[1].legend(np.append(col_ToPlot.values,flux_org), loc='lower left')
    if ('ms' not in suffix): #If not for ms, annotate on plots
        sub[1].annotate('Cumulative sum of observed fluxes (blue line) plus filled real gaps', (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
    sub[1].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=True, bottom=True, top=False, left=True, right=True)

    # Own legend: Adding scatter plot to print a legend
    #sub[1].scatter([], [], color='blue', label='Flux')
    for i_k, key in enumerate(usedKeys):
        sub[1].scatter([], [], color=techColumns[key], label=key)
    sub[1].legend(loc='best', frameon=True, handlelength=0.75)

    # Format the date ticks, see https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/date.html
    months = pmd.MonthLocator()  # every month
    months_fmt = pmd.DateFormatter('%Y-%m')
    days = pmd.DayLocator()
    sub[numRows-1].xaxis.set_major_locator(months)
    sub[numRows-1].xaxis.set_major_formatter(months_fmt)
    sub[numRows-1].xaxis.set_minor_locator(days)
    

    # Set title and label
    fig.suptitle('Time series plots of daily sums', x=0.5, y=0.95, fontsize=fsize+4)
    ytitle = ini.FluxGas + ' ('+ini.ConvSums.replace('period','day')+')'
    sub[0].set_ylabel(ytitle, fontsize=fsize, weight='bold')
    ytitle = ini.FluxGas + ' ('+ini.ConvSums+')'
    sub[1].set_ylabel(ytitle, fontsize=fsize, weight='bold')
    sub[numRows-1].set_xlabel('Time (days)', fontsize=fsize, weight='bold')
    
    if suffix != '': suffix = '_' + suffix 
    mgf.utl.save_plot('pds', ini, file_path, run_number, suffix='daily_'+scenario+suffix)    
