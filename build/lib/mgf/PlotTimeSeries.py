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
import matplotlib.lines as mlines
import numpy as np

### Own imports
import mgf

"""
TESTING:
    scenario='hhs';
"""

#%%
def plot_series(data, ini, scenario, file_path='', run_number='', suffix=''):
    # plot_series():    Plot time series of fluxes, residuals, three types of sums (data, gaps, all)
    # data:             Dataframe with timestamp index, flux measurements, (maybe meteo,) and filled fluxes for various techniques
    # ini:              Settings from ini file
    # scenario:         Artificial gap filling scenarioo
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run
    # suffix:           Suffix to filename
    
#%%
  # Set data for plot
    flux_org = mgf.utl.name_flux(ini)
    col_ToPlot = data.filter(regex=scenario).columns
    data_sub = data[np.append(flux_org,col_ToPlot.values)]#[2000:3000]
    # All data points with measured data, i.e. artificial gaps
    data_meas_all = data_sub.copy()
    data_meas_all.loc[~data[flux_org].notnull(),] = np.nan #Subset filled with NAs but with same length as original data
    # All data points with measured data, i.e. artificial gaps, but only with complete set of GFT
    data_meas_cplt = data_sub.copy()
    v_cplt_meas = data[flux_org].notnull() & mgf.utl.check_completeness(data, scenario)
    data_meas_cplt.loc[~v_cplt_meas,] = np.nan #Subset filled with NAs but with same length as original data
    # All data points at gaps
    data_gaps_all = data_sub.copy()
    v_gaps_all = ~data[flux_org].notnull()
    data_gaps_all.loc[~v_gaps_all,] = np.nan #Subset filled with NAs but with same length as original data
    

    # Plot properties
    numCols = 1
    numRows = 2
    fig_width = 15 #12.0
    fig_height = fig_width / 2 * numRows #20.0
    fig, sub = plt.subplots(nrows=numRows, ncols=numCols, figsize=(fig_width, fig_height))
    fig.subplots_adjust(hspace=0.03)
    fig.patch.set_facecolor('white')
    fsize=12

    # Set colors
    labels = data_meas_all.columns.str.replace('_'+scenario,'')
    lineColors, techColumns, usedKeys = mgf.utl.set_colors('gft_rainbow',labels[1:])
    
    if suffix == 'fluxes_all':
        # Time series fluxes
        sub[0].set_prop_cycle('color', lineColors)
        sub[0].plot(data_meas_all.iloc[:,1:], '+', markersize=1)
        sub[0].plot(data_meas_all[flux_org], 'bo', markersize=2)
        sub[0].annotate('Observed and predicted fluxes at artifical gaps', (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
        sub[0].tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True, bottom=False, top=True, left=True, right=True)
           
        # Time series residuals
        sub[1].set_prop_cycle('color', lineColors)
        sub[1].plot(data_meas_all.iloc[:,1:].subtract(data_meas_all[flux_org], axis=0), '+', markersize=1)
        sub[1].plot(data_meas_all[flux_org].subtract(data_meas_all[flux_org], axis=0), 'bo', markersize=2)
        sub[1].annotate('Residuals of fluxes at artifical gaps', (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
        sub[1].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=True, bottom=True, top=False, left=True, right=True)

        # Set title and label
        fig.suptitle('Time series plots of artificial \''+scenario+'\' scenarios\n (all artificial gaps, even if not filled by some of the GFTs)', x=0.5, y=0.95, fontsize=fsize+4)
        ytitle = ini.FluxGas+' ('+ini.FluxUnit+')'
        sub[0].set_ylabel(ytitle, fontsize=fsize, weight='bold')
        ytitle = ini.FluxGas+' ('+ini.FluxUnit+')'
        sub[1].set_ylabel(ytitle, fontsize=fsize, weight='bold')

    elif suffix == 'fluxes_cplt':
        # Time series fluxes
        sub[0].set_prop_cycle('color', lineColors)
        sub[0].plot(data_meas_cplt.iloc[:,1:], '+', markersize=1)
        sub[0].plot(data_meas_cplt[flux_org], 'bo', markersize=2)
        sub[0].annotate('Observed and predicted fluxes at artifical gaps', (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
        sub[0].tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True, bottom=False, top=True, left=True, right=True)
           
        # Time series residuals
        sub[1].set_prop_cycle('color', lineColors)
        sub[1].plot(data_meas_cplt.iloc[:,1:].subtract(data_meas_cplt[flux_org], axis=0), '+', markersize=1)
        sub[1].plot(data_meas_cplt[flux_org].subtract(data_meas_cplt[flux_org], axis=0), 'bo', markersize=2)
        sub[1].annotate('Residuals of fluxes at artifical gaps', (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
        sub[1].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=True, bottom=True, top=False, left=True, right=True)
        
        # Set title and label
        fig.suptitle('Time series plots of artificial \''+scenario+'\' scenarios\n (only artificial gaps filled with complete suite of GFTs)', x=0.5, y=0.95, fontsize=fsize+4)
        ytitle = ini.FluxGas+' ('+ini.FluxUnit+')'
        sub[0].set_ylabel(ytitle, fontsize=fsize, weight='bold')
        ytitle = ini.FluxGas+' ('+ini.FluxUnit+')'
        sub[1].set_ylabel(ytitle, fontsize=fsize, weight='bold')

    elif suffix == 'fluxes_real': #Real gaps filled with 'hhs' scenarios
    
        # Time series fluxes
        sub[0].set_prop_cycle('color', lineColors)
        sub[0].plot(data_sub.iloc[:,1:], '+', markersize=1)
        sub[0].plot(data_sub[flux_org], 'bo', markersize=2)
        sub[0].annotate('Observed and predicted fluxes at all data points', (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
        sub[0].tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True, bottom=False, top=True, left=True, right=True)
  
        # Time series fluxes
        sub[1].set_prop_cycle('color', lineColors)
        sub[1].plot(data_gaps_all.iloc[:,1:], '+', markersize=1)
        sub[1].plot(data_gaps_all[flux_org], 'bo', markersize=2)
        sub[1].annotate('Predicted fluxes at real gaps', (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
        sub[1].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=True, bottom=True, top=False, left=True, right=True)

        # Set title and label
        fig.suptitle('Time series plots of artificial \''+scenario+'\' scenarios and real gaps', x=0.5, y=0.95, fontsize=fsize+4)
        ytitle = ini.FluxGas+' ('+ini.FluxUnit+')'
        sub[0].set_ylabel(ytitle, fontsize=fsize, weight='bold')
        ytitle = ini.FluxGas+' ('+ini.FluxUnit+')'
        sub[1].set_ylabel(ytitle, fontsize=fsize, weight='bold')

    elif suffix == 'cums_cplt':     
        # Cumulative sums of data points with measurements
        #filled_meas = data_meas_all.notnull().sum(axis=1) == data_meas_all.shape[1]
        sub[0].set_prop_cycle('color', lineColors)
        sub[0].plot(float(ini.ConvFactor)*data_meas_cplt.iloc[:,1:].cumsum(skipna=True), '-', linewidth=0.5) # linestyle='solid'
        sub[0].plot(float(ini.ConvFactor)*data_meas_cplt[flux_org].cumsum(skipna=True), 'b-', linewidth=2)
        sub[0].annotate('Cumulative observed and predicted fluxes at artifical gaps', (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
        sub[0].tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True, bottom=False, top=True, left=True, right=True)
       
        # Cumulative sums of residuals
        sub[1].set_prop_cycle('color', lineColors)
        sub[1].plot(float(ini.ConvFactor)*data_meas_cplt.iloc[:,1:].subtract(data_meas_cplt[flux_org], axis=0).cumsum(skipna=True), '-', linewidth=0.5) # linestyle='solid'
        sub[1].plot(float(ini.ConvFactor)*(data_meas_cplt[flux_org]-data_meas_cplt[flux_org]).cumsum(skipna=True), 'b-', linewidth=2)
        sub[1].annotate('Cumulative sums of residuals at artifical gaps', (0.5, 0.1), xycoords='axes fraction', ha='center', va='center', fontsize=fsize-2, fontweight='bold')
        sub[1].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=True, bottom=True, top=False, left=True, right=True)
  
        # Set title and label
        fig.suptitle('Time series plots of artificial \''+scenario+'\' scenarios\n (only artificial gaps filled with complete suite of GFTs)', x=0.5, y=0.95, fontsize=fsize+4)
        ytitle = ini.FluxGas+' ('+ini.ConvSums+')'
        sub[0].set_ylabel(ytitle, fontsize=fsize, weight='bold')
        ytitle = ini.FluxGas+' ('+ini.ConvSums+')' 
        sub[1].set_ylabel(ytitle, fontsize=fsize, weight='bold')
        
    else:
        print("!!! Interrupt plot_series() plot. The following suffix not found:", suffix)
        return
                
    # Own legend: Adding scatter plot to print a legend
    #sub[2].scatter([], [], color='blue', label='Flux')
    for i_k, key in enumerate(usedKeys):
        sub[1].scatter([], [], color=techColumns[key], label=key)
    sub[1].legend(loc='best', frameon=True, handlelength=0.75)

    # format the date ticks, see https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/date.html
    months = pmd.MonthLocator()  # every month
    months_fmt = pmd.DateFormatter('%Y-%m')
    days = pmd.DayLocator()
    sub[numRows-1].xaxis.set_major_locator(months)
    sub[numRows-1].xaxis.set_major_formatter(months_fmt)
    sub[numRows-1].xaxis.set_minor_locator(days)
    sub[numRows-1].set_xlabel('Time (half-hours)', fontsize=fsize, weight='bold')
   
    if suffix != '': suffix = suffix + '_'
    mgf.utl.save_plot('pts', ini, file_path, run_number, suffix=suffix+scenario)    
