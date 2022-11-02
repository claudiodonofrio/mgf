#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py

Functions:
    # boot_min_max():       Calculate min/max across result frames to have the same axis limits
    # draw_boxes():         Draw boxes with bootstrapping results for each error measure
    # plot_bootstat_single(): Plot bootstrapping results for single column
    # plot_bootstat_double(): Plot bootstrapping results for double columns
    # plot_bootstat_triple(): Plot bootstrapping results for triple columns

Notes:
    Box plot code based on from matplotlib gallery (https://matplotlib.org/gallery/statistics/boxplot_demo.html)
    For statistics and plots like R2, the model (MR) has to be on x-axis.

"""
#%% Initialization
### Imports from python
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
import pandas as pd
### Own imports
import mgf

"""
For TESTING:
    boot_1=boot_hhs_ft; error_type='Bias'; scens=['hhs']; colors='gft_rainbow'
    boot_1=boot_hhs_ft; boot_2=boot_days_ft; scens=['hhs','days']; times=['full-time','full-time']; colors='gft_rainbow'
    boot_1=boot_hhs_ft; boot_2=boot_hhs_dt; boot_3=boot_hhs_nt; scens=['hhs','hhs','hhs']; times=['day-time','night-time','full-time']
    data_plt = boot_1.loc['Bias']; error_type='Bias'; scenario='hhs'; colors=['limegreen', 'teal']; fsize=12; zoom=np.nan
    colors='gft_rainbow'
"""

#%%
def boot_min_max(boot_1, boot_2='', boot_3=''):
    # boot_min_max():    Calculate min/max across result frames to have the same axis limits
    # boot_1,_2,_3:      Result data frame(s) with bootstrapping results (should be same size)
    
    # Quick and dirty solution: variables used globally... 
    def_errors = ['Bias','SDev','R2']
    global all_max, all_min
    all_min={}
    all_max={}
    for i_e in def_errors:
        if not isinstance(boot_2, pd.DataFrame): #then single column
            all_min[i_e] = min(boot_1.loc[i_e].min())
            all_max[i_e] = max(boot_1.loc[i_e].max())
        elif not isinstance(boot_3, pd.DataFrame): #then double columns
            all_min[i_e] = min(boot_1.loc[i_e].min().min(),boot_2.loc[i_e].min().min())
            all_max[i_e] = max(boot_1.loc[i_e].max().max(),boot_2.loc[i_e].max().max())
        else: #Triple column
            all_min[i_e] = min(boot_1.loc[i_e].min().min(),boot_2.loc[i_e].min().min(),boot_3.loc[i_e].min().min())
            all_max[i_e] = max(boot_1.loc[i_e].max().max(),boot_2.loc[i_e].max().max(),boot_3.loc[i_e].max().max())
    
#%%
def draw_boxes(plx, data_plt, error_type, scenario, colors, fsize, unit):
    # draw_boxes():     Draw boxes with bootstrapping results for each error measure
    # plx:              Plot frame
    # data_plt:         Dataframe with data to box plot
    # error_type:       'Bias', 'SDev', or 'R2', (currently not used 'RMSE', 'Residuals')
    # scenario:         Name of gap filling scenario
    # colors:           Color of bars, can be fixed (colors='gft_rainbow') or multiple colors (colors=['limegreen', 'teal'])
    # fsize:            Basic font size (12)
    # unit:             Name of unit (for labels)

    # Set limits of plots
    y_min = np.nan
    y_max = np.nan
    if error_type == 'Bias':
        lim = max(abs(all_min['Bias']),abs(all_max['Bias']))
        y_range = +lim - (-lim)
        y_min = -lim - 0.01*y_range
        y_max = +lim - 0.01*y_range
        error_text = 'Bias'+' ('+unit+')'
    if error_type == 'SDev':
        y_range = all_max['SDev'] - all_min['SDev']
        y_min = all_min['SDev'] - 0.01*y_range
        y_max = all_max['SDev'] + 0.01*y_range
        error_text = 'SDev'+' ('+unit+')'      
    if error_type == 'R2':
        y_min = 0.0 - 0.02
        y_max = 1.0 + 0.02
        error_text = 'R-squared'
#    if error_type == 'RMSE':
#        y_range = all_max['RMSE'] - all_min['RMSE']
#        y_min = all_min['RMSE'] - 0.01*y_range
#        y_max = all_max['RMSE'] + 0.01*y_range
#        error_text = 'RMSE'
#    if error_type == 'Residuals':
#        lim = max(abs(data_plt.min().min()),abs(data_plt.max().max()))
#        y_min = -lim - 0.2
#        y_max = +lim + 0.2
#        error_text = 'Residuals'

    #Transpose data
    data_plt_t = data_plt.transpose()
    numBoxes = data_plt_t.shape[0]
    data_array = data_plt_t.values.tolist()
    
    #Plot data to boxplots
    flier=dict(markeredgecolor='gray', markerfacecolor='none', marker='o', markersize=4, alpha=0.5)
    meds=dict(color='black')
    # Standard box plot properties: quartiles, whiskers at [10, 90] percentiles, whisker at 1.5 quartile, all outliers
    bp = plx.boxplot(data_array, notch=0, vert=1, flierprops=flier, medianprops=meds, whis=[10, 90], zorder=2) #sym='+', 
    bp = plx.boxplot(data_array, notch=0, vert=1, flierprops=flier, medianprops=meds, whis=[10, 90], zorder=2, showfliers=False) #sym='+', 

    # Set the axes ranges and axes labels
    plx.set_ylabel(error_text, fontsize=fsize, weight='bold')
  #  plx.tick_params(labelsize=fsize)
    labels = data_plt_t.index.str.replace('_'+scenario,'')
    plx.set_xticks(range(1,len(labels)+1)) #Needed to avoid error: ValueError: The number of FixedLocator locations (32), usually from a call to set_ticks, does not match the number of ticklabels (16).
    plx.set_xticklabels(labels, rotation=45, ha='right', fontsize=fsize)
   # plx.xticks(labels, rotation=45, ha='right', fontsize=fsize)
    if np.isfinite(y_min): plx.set_ylim(y_min, y_max)
    
    # Horizontal grid behind plot objects
    plx.yaxis.grid(True, linestyle='-', which='major', color='gray', linewidth=0.2) #color='lightgrey', alpha=0.5)
    plx.axhline(0, linestyle='-', color='black', linewidth=0.5, zorder=1)
    plx.set_axisbelow(True)

    # Set colors for techniques or user input
    boxColors = mgf.utl.set_colors(colors, labels)[0]
    
    # Fill boxes with colors and additional markers for mean and [5,95] percentiles
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = np.column_stack([boxX, boxY])
        # Alternate between provided colors
        k = i % len(boxColors)
        boxPolygon = Polygon(boxCoords, facecolor=boxColors[k], zorder=2)
        plx.add_patch(boxPolygon)
        
        # Now use position of the medians to draw other measures (mean and percentiles)
        med = bp['medians'][i]
        # With horizontal alignment and in the center of each box
        plx.plot([np.average(med.get_xdata())], [np.average(data_array[i])], color='k', marker='*', zorder=3)
        # Percentiles
        plx.plot([np.average(med.get_xdata())], [np.percentile(data_array[i],95)], color='k', marker='o', markersize=4, zorder=3)
        plx.plot([np.average(med.get_xdata())], [np.percentile(data_array[i],5)], color='k', marker='o', markersize=4, zorder=3)

#%%   
def plot_bootstats_single(boot_1, ini, scens, times, colors, file_path='', run_number='', suffix=''): 
    # plot_bootstats_single(): Plot bootstrapping results for single column
    # boot_1            Result data frame(s) with bootstrapping results
    # ini:              Settings from ini file   
    # scens:            Scenario for 1 (array)
    # times:            Time of day setting for 1 (array)
    # colors:           Color of bars, can be fixed (colors='gft_rainbow') or multiple colors (colors=['limegreen', 'teal'])
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run
#%%     
    # Figure settings
    numBoxes = boot_1.shape[1]
    num_columns = 1
    numErrors = 3
    fig_height = numErrors * 5.0
    box_spacing = 0.4 #Here figure width varies and box spacing constant
    fig_width = num_columns * numBoxes * box_spacing
    if (fig_width <= 5.0): fig_width = 5.0
    fsize = 12
    fig, plx = plt.subplots(figsize=(fig_width, fig_height))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    fig.subplots_adjust(hspace=0.01)
    fig.patch.set_facecolor('white')
    boot_min_max(boot_1) #Determine min max overall single scenario
    
    # Plot data
    scenario = scens[0]
    ax1 = plt.subplot(311)
    ax1.xaxis.set_visible(False)
    draw_boxes(ax1, data_plt = boot_1.loc['R2'], error_type='R2', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax2 = plt.subplot(312)
    ax2.xaxis.set_visible(False)
    draw_boxes(ax2, data_plt = boot_1.loc['SDev'], error_type='SDev', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax3 = plt.subplot(313)
    draw_boxes(ax3, data_plt = boot_1.loc['Bias'], error_type='Bias', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)

    # Plot title and lable
    fig.suptitle('Bootstrapping performances\n for '+ini.FluxGas, x=0.5, y=0.97, fontsize=fsize+4)
    ax1.set_title(scens[0]+'\n('+times[0]+')', fontsize=fsize+2, weight='bold')
    ax3.set_xlabel('Gap-filling techniques', fontsize=fsize, weight='bold')

    # Finally, add a basic legend
    #fig.text(0.84, 0.07, 'Quality      qc=0',  backgroundcolor=boxColors[0], color='black', weight='normal', size=fsize)
    #fig.text(0.84, 0.02, 'Quality qc=0&1', backgroundcolor=boxColors[1], color='white', weight='normal', size=fsize)
    ##For top position, use 0.75, 0.85
#%%
    if suffix != '': suffix = '_' + suffix 
    mgf.utl.save_plot('pba', ini, file_path, run_number, suffix='single_'+scens[0]+'_'+times[0]+suffix) # pba - plot of bootstrap analysis


#%%   
def plot_bootstats_double(boot_1, boot_2, ini, scens, times, colors, file_path='', run_number='', suffix=''):    
    # plot_bootstats_triple():     Plot bootstrapping results for double columns
    # boot_1,_2:        Result data frame(s) with bootstrapping results (should be same size)
    # ini:              Settings from ini file   
    # scens:            Scenarios for 1, 2
    # times:            Time of day setting for 1, 2
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run
    # suffix:           Suffix to filename
    
#%%   
    # Figure settings
    numBoxes = boot_1.shape[1]
    num_columns = 2
    numErrors = 3
    fig_height = numErrors * 5.0
    #Figure dimensions:
    #   figure height=16, actual size: 13", probably -3" because of bbox_inches=tight, in ms: 9.25"
    #   figure width=8.5, actual size: 8.5", in ms 6.25"
    #   Bigger figures lead to smaller lines --> wanted, i.e. height = 20 or width = 12 are good sizes
    fig_width = 15.0 #Here figure width constant and box spacing varies
    box_spacing = fig_width / num_columns / numBoxes    
    fsize = 12
    fig, plx = plt.subplots(figsize=(fig_width, fig_height))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    fig.subplots_adjust(hspace=0.00, wspace=0.03)
    fig.patch.set_facecolor('white')
    boot_min_max(boot_1, boot_2) #Determine min max overall all scenarios
    if mgf.__MS__ == True and ini.FluxGas == 'CO2':
        all_min['SDev'] = min(all_min['SDev'],1.012306654-0.1)
    
    # Plot data boot_1
    data_res = boot_1
    scenario = scens[0]
    ax11 = plt.subplot(321)
    ax11.xaxis.set_visible(False)
    draw_boxes(ax11, data_plt = data_res.loc['R2'], error_type='R2', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax12 = plt.subplot(323)
    ax12.xaxis.set_visible(False)
    draw_boxes(ax12, data_plt = data_res.loc['SDev'], error_type='SDev', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax13 = plt.subplot(3,2,5)
    draw_boxes(ax13, data_plt = data_res.loc['Bias'], error_type='Bias', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)

    # Plot data boot_2
    data_res = boot_2
    scenario = scens[1]
    ax21 = plt.subplot(322)
    ax21.xaxis.set_visible(False)
    ax21.yaxis.tick_right()
    ax21.yaxis.label.set_visible(False)
    #ax21.yaxis.set_label_position('right') #or show on righthand side
    draw_boxes(ax21, data_plt = data_res.loc['R2'], error_type='R2', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax22 = plt.subplot(324)
    ax22.xaxis.set_visible(False)
    ax22.yaxis.tick_right()
    ax22.yaxis.label.set_visible(False)
    draw_boxes(ax22, data_plt = data_res.loc['SDev'], error_type='SDev', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax23 = plt.subplot(326)
    ax23.yaxis.tick_right()
    ax23.yaxis.label.set_visible(False)
    draw_boxes(ax23, data_plt = data_res.loc['Bias'], error_type='Bias', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    
    #Add results of GFC to CO2 for ms (mean,P_10,P_90 for DE2000)
    if mgf.__MS__ == True and ini.FluxGas == 'CO2':
        e_R2 = 0.710813355
        e_SDev = 2.380878059
        e_Bias = 0.030531425
        ax11.errorbar(7.5,e_R2, yerr=[[e_R2 - 0.4158596],[0.9001154 - e_R2]], marker='*', color='red', ms=12, ecolor='red', elinewidth=2, capsize=4)
        ax12.errorbar(7.5,e_SDev, yerr=[[e_SDev - 1.012306654],[3.621322929 - e_SDev]], marker='*', color='red', ms=12, ecolor='red', elinewidth=2, capsize=4)
        ax13.errorbar(7.5,e_Bias, yerr=[[e_Bias - -0.1526872],[0.1822482 - e_Bias]], marker='*', color='red', ms=12, ecolor='red', elinewidth=2, capsize=4)
        
    # Plot title and lable
    fig.suptitle('Bootstrapping performances for '+ini.FluxGas, x=0.5, y=0.95, fontsize=fsize+4)
    ax11.set_title(scens[0]+'\n('+times[0]+')', fontsize=fsize+2, weight='bold')
    ax21.set_title(scens[1]+'\n('+times[1]+')', fontsize=fsize+2, weight='bold')
    ax23.set_xlabel('Gap-filling techniques', fontsize=fsize, weight='bold', x=0)
#%%
    if suffix != '': suffix = '_' + suffix
    scens_name = scens[0]
    if scens[0] != scens[1]: scens_name = 'both'
    mgf.utl.save_plot('pba', ini, file_path, run_number, suffix='double_'+scens_name+'_'+times[0]+suffix)   


#%%   
def plot_bootstats_triple(boot_1, boot_2, boot_3, ini, scens, times, colors, file_path='', run_number='', suffix=''):    
    # plot_bootstats_triple():     Plot bootstrapping results for triple columns
    # boot_1,_2,_3:     Result data frame(s) with bootstrapping results (should be same size)
    # ini:              Settings from ini file   
    # scens:            Scenarios of 1, 2, 3
    # times:            Time of day setting for 1, 2, 3
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run
#%%   
    # Figure settings
    numBoxes = boot_1.shape[1]
    num_columns = 3
    numErrors = 3
    fig_height = numErrors * 5.0
    fig_width = 15.0 #Here figure width constant and box spacing varies
    box_spacing = fig_width / num_columns / numBoxes    
    fsize = 12
    fig, plx = plt.subplots(figsize=(fig_width, fig_height))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    fig.subplots_adjust(hspace=0.00, wspace=0.03)
    fig.patch.set_facecolor('white')
    boot_min_max(boot_1, boot_2, boot_3) #Determine min max overall all scenarios
    
    # Plot data boot_1
    data_res = boot_1
    scenario = scens[0]
    ax11 = plt.subplot(331)
    ax11.xaxis.set_visible(False)
    draw_boxes(ax11, data_plt = data_res.loc['R2'], error_type='R2', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax12 = plt.subplot(334)
    ax12.xaxis.set_visible(False)
    draw_boxes(ax12, data_plt = data_res.loc['SDev'], error_type='SDev', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax13 = plt.subplot(337)
    draw_boxes(ax13, data_plt = data_res.loc['Bias'], error_type='Bias', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)

    # Plot data boot_2
    data_res = boot_2
    scenario = scens[1]
    ax21 = plt.subplot(332)
    ax21.xaxis.set_visible(False)
    #ax21.yaxis.set_visible(False) #Also turns off grid :(
    ax21.yaxis.label.set_visible(False) # Turn off yaxis name
    ax21.yaxis.set_ticklabels([]) # Turn of yaxis numbers
    ax21.yaxis.set_ticks_position('none') # Turn off ticks
    draw_boxes(ax21, data_plt = data_res.loc['R2'], error_type='R2', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax22 = plt.subplot(335)
    ax22.xaxis.set_visible(False)
    #ax22.yaxis.set_visible(False) #Also turns off grid :(
    ax22.yaxis.label.set_visible(False) # Turn off yaxis name
    ax22.yaxis.set_ticklabels([]) # Turn of yaxis numbers
    ax22.yaxis.set_ticks_position('none') # Turn off ticks
    draw_boxes(ax22, data_plt = data_res.loc['SDev'], error_type='SDev', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax23 = plt.subplot(338)
    #ax23.yaxis.set_visible(False) #Also turns off grid :(
    ax23.yaxis.label.set_visible(False) # Turn off yaxis name
    ax23.yaxis.set_ticklabels([]) # Turn of yaxis numbers
    ax23.yaxis.set_ticks_position('none') # Turn off ticks
    draw_boxes(ax23, data_plt = data_res.loc['Bias'], error_type='Bias', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)

    # Plot data boot_3
    data_res = boot_3
    scenario = scens[2]
    ax31 = plt.subplot(333)
    ax31.xaxis.set_visible(False)
    ax31.yaxis.tick_right()
    ax31.yaxis.label.set_visible(False)
    draw_boxes(ax31, data_plt = data_res.loc['R2'], error_type='R2', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax32 = plt.subplot(336)
    ax32.xaxis.set_visible(False)
    ax32.yaxis.tick_right()
    ax32.yaxis.label.set_visible(False)
    draw_boxes(ax32, data_plt = data_res.loc['SDev'], error_type='SDev', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)
    ax33 = plt.subplot(339)
    ax33.yaxis.tick_right()
    ax33.yaxis.label.set_visible(False)
    draw_boxes(ax33, data_plt = data_res.loc['Bias'], error_type='Bias', scenario=scenario, colors=colors, fsize=fsize, unit=ini.FluxUnit)

    # Plot title and lable
    fig.suptitle('Bootstrapping performances for '+ini.FluxGas, x=0.5, y=0.95, fontsize=fsize+4)
    ax11.set_title(scens[0]+'\n('+times[0]+')', fontsize=fsize+2, weight='bold')
    ax21.set_title(scens[1]+'\n('+times[1]+')', fontsize=fsize+2, weight='bold')
    ax31.set_title(scens[2]+'\n('+times[2]+')', fontsize=fsize+2, weight='bold')
    ax23.set_xlabel('Gap-filling techniques', fontsize=fsize, weight='bold')

#%%
    if suffix != '': suffix = '_' + suffix
    times_name = times[0]
    if times[0] == 'full-time' and times[1] == 'day-time' and times[2] == 'night-time': times_name = 'fdn'
    mgf.utl.save_plot('pba', ini, file_path, run_number, suffix='triple_'+scens[0]+'_'+times_name+suffix)
