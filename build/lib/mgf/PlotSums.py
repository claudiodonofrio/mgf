#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py

Functions:
    # plot_sums():      Plot fluxes and errors summed up over the whole period (for all techniques or ensemble only)
    # plot_sums_ens():  Plot fluxes and errors summed up over the whole period of all gft (left plot) and ensemble gft (right plot)
        
Notes:
    Some code inspired here
        http://codewithmax.com/2018/03/17/plotting-error-bars-in-python-using-matplotlib-and-numpy-random/
"""

#%% Initialization
### Imports from python
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

### Own imports
import mgf

"""
For TESTING:
    colors='gft_rainbow'
"""

#%%
def plot_sums(sums, ini, colors, file_path='', run_number='', suffix=''):
    # plot_sums():      Plot fluxes and errors summed up over the whole period (for all techniques or ensemble only)
    # sums:             Results dataframe with sums
    # ini:              Settings from ini file
    # colors:           Color of bars, can be fixed colors (colors='gft_rainbow' or 'gft_gray') or multiple colors (colors=['limegreen', 'teal'])
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run
    # suffix:           Suffix to filename

 #%%   
    num_dots = sums.shape[0]
    fig_width = num_dots*0.5
    fig_height = 5.0
    fsize = 12
    
    #Calculate y limits     
    y_min = min(min(sums['SumTotal']-sums['ErrorTotal']), min(sums['SumObs']))
    y_max = max(max(sums['SumTotal']+sums['ErrorTotal']), max(sums['SumObs']))
    y_range = +y_max - y_min
    y_min = y_min - 0.05*y_range
    y_max = y_max + 0.05*y_range
    
    fig, plx = plt.subplots(figsize=(fig_width, fig_height))
    labels = sums.index
    fig.patch.set_facecolor('white')
    fig.suptitle('Aggregated fluxes (sums)', x=0.5, y=0.95, fontsize=fsize+2)
    dotColors, techColumns, usedKeys = mgf.utl.set_colors(colors,labels)
    c_bline = 'blue'
    if colors=='gft_gray': c_bline='black'
    # To have different color points
    for i_l, label in enumerate(labels):
        plx.errorbar(i_l,sums['SumTotal'][label], yerr=sums['ErrorTotal'][label], marker='o', color='k', ms=8, mfc=dotColors[i_l],
                     ecolor='0.4', elinewidth=1, capsize=4)
        plx.errorbar(i_l,sums['SumTotal'][label], yerr=sums['BiasGaps'][label], marker=' ', 
                     ecolor='black', elinewidth=1, capsize=4)
    #Add baseline with sum of observed data
    plx.plot([-0.5,num_dots-1+0.5], [min(sums['SumObs']), max(sums['SumObs'])], linestyle='dashed', color=c_bline, linewidth=2)
    
    # Own legend: Adding scatter plot to print a legend
    if suffix != 'ens': #Do not add legend for ens plot
        for i_k, key in enumerate(usedKeys):
            plx.scatter([], [], edgecolors='k', facecolors=techColumns[key], label=key)
        plx.legend(loc='best', frameon=True)

    #Add box for ensemble to plot
    if suffix == 'ens':
        # Add outer rectangle (x,y),w,h
        r_x = -0.5
        r_y = sums['LowerCI'].min(axis=0)
        r_w = num_dots
        r_h = sums['UpperCI'].max(axis=0) - r_y
        plx.add_patch(Rectangle((r_x,r_y),r_w,r_h, facecolor='0.95'))
        # Add inner rectangle (x,y),w,h
        r_x = -0.5
        r_y = sums['SumTotal'].min(axis=0)
        r_w = num_dots
        r_h = sums['SumTotal'].max(axis=0) - r_y
        plx.add_patch(Rectangle((r_x,r_y),r_w,r_h, facecolor='0.85'))
        
    # Set axis labels and ticks
    xlabels = sums.index
    ytitle = ini.FluxGas + ' ('+ini.ConvSums+')'
    plx.set_xticks(range(num_dots))
    plx.set_xticklabels(xlabels, rotation=45, ha='right')
    plx.set_xlabel('Gap filling techniques', fontsize=fsize, weight='bold')
    plx.set_ylim(y_min, y_max)
    plx.set_ylabel(ytitle, fontsize=fsize, weight='bold')
    plx.yaxis.grid(True, linestyle='-', which='major', color='gray', linewidth=0.2)
    
#%%
    if suffix != '': suffix = '_' + suffix 
    mgf.utl.save_plot('pse', ini, file_path, run_number, suffix='real'+suffix)

#%%    
def plot_sums_ens(sums, sums_ens, res_ens, ini, colors, file_path='', run_number='', suffix=''):
    # plot_sums_ens():  Plot fluxes and errors summed up over the whole period of all gft (left plot) and ensemble gft (right plot)
    # sums:             Results dataframe with all sums
    # sums_ens:         Dataframe with sums of ensemble gft only
    # res_ens:          Dataframe with ensemble results
    # ini:              Settings from ini file
    # colors:           Color of bars, can be fixed colors (colors='gft_rainbow' or 'gft_gray') or multiple colors (colors=['limegreen', 'teal'])
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run
    # suffix:           Suffix to filename

 #%%   
 
    # Plot properties
    num_dots = sums.shape[0]
    num_dots_ens = sums_ens.shape[0]
    fig_width = (num_dots + num_dots_ens)*0.35
    fig_height = 5.0
    fsize = 12
    numCols = 2
    numRows = 1

    # Make plot
    fig, sub = plt.subplots(nrows=numRows, ncols=numCols, figsize=(fig_width, fig_height),
                            gridspec_kw={'width_ratios': [num_dots, num_dots_ens]})
    fig.subplots_adjust(hspace=0.03, wspace=0.03) #Needs to be after gridspec command...
    fig.patch.set_facecolor('white')
    #fig.suptitle('Aggregated fluxes (sums) of ensemble', x=0.5, y=0.95, fontsize=fsize+2)
    
    #Calculate y limits     
    y_min = min(min(sums['SumTotal']-sums['ErrorTotal']), min(sums['SumObs']))
    y_max = max(max(sums['SumTotal']+sums['ErrorTotal']), max(sums['SumObs']))
    y_range = +y_max - y_min
    y_min = y_min - 0.05*y_range
    y_max = y_max + 0.05*y_range
    # Make labels and colors
    labels = sums.index
    dotColors, techColumns, usedKeys = mgf.utl.set_colors(colors,labels)
    c_bline = 'blue'
    if colors=='gft_gray': c_bline='black'
    # Own legend: Adding scatter plot to print a legend
    for i_k, key in enumerate(usedKeys):
        sub[0].scatter([], [], edgecolors='k', facecolors=techColumns[key], label=key)
    sub[0].legend(loc='best', frameon=True)
    # To have different color points
    for i_l, label in enumerate(labels):
        sub[0].errorbar(i_l,sums['SumTotal'][label], yerr=sums['ErrorTotal'][label], marker='o', color='k', ms=8, mfc=dotColors[i_l],
                     ecolor='0.4', elinewidth=1, capsize=4)
        sub[0].errorbar(i_l,sums['SumTotal'][label], yerr=sums['BiasGaps'][label], marker=' ', 
                     ecolor='black', elinewidth=1, capsize=4)
    #Add baseline with sum of observed data
    sub[0].plot([-0.5,num_dots-1+0.5], [min(sums['SumObs']), max(sums['SumObs'])], linestyle='dashed', color=c_bline, linewidth=2)
    
    # Repeat for ensemble sum
    #Calculate y limits     
    y_min_ens = min(min(sums_ens['SumTotal']-sums_ens['ErrorTotal']), min(sums_ens['SumObs']))
    y_max_ens = max(max(sums_ens['SumTotal']+sums_ens['ErrorTotal']), max(sums_ens['SumObs']))
    y_min_ens = min(sums_ens['SumTotal']-sums_ens['ErrorTotal'])
    y_max_ens = max(sums_ens['SumTotal']+sums_ens['ErrorTotal'])
    y_range_ens = +y_max_ens - y_min_ens
    y_min_ens = y_min_ens - 0.05*y_range_ens
    y_max_ens = y_max_ens + 0.05*y_range_ens
    # Make labels and colors
    labels = sums_ens.index
    dotColors, techColumns, usedKeys = mgf.utl.set_colors(colors,labels)
    c_bline = 'blue'
    if colors=='gft_gray': c_bline='black'
    for i_l, label in enumerate(labels):
        sub[1].errorbar(i_l,sums_ens['SumTotal'][label], yerr=sums['ErrorTotal'][label], marker='o', color='k', ms=8, mfc=dotColors[i_l],
                     ecolor='0.4', elinewidth=1, capsize=4)
        sub[1].errorbar(i_l,sums_ens['SumTotal'][label], yerr=sums['BiasGaps'][label], marker=' ', 
                     ecolor='black', elinewidth=1, capsize=4)
    #Add baseline with sum of observed data
   #sub[1].plot([-0.5,num_dots_ens-1+0.5], [min(sums_ens['SumObs']), max(sums_ens['SumObs'])], linestyle='dashed', color=c_bline, linewidth=2)
    
    #Add box for ensemble to plot
    # Add outer rectangle (x,y),w,h
    r_x = -0.5
    r_y = res_ens['LowerCI']
    r_w = num_dots_ens
    r_h = res_ens['UpperCI'] - r_y
    sub[1].add_patch(Rectangle((r_x,r_y),r_w,r_h, facecolor='0.95'))
    # Add inner rectangle (x,y),w,h
    r_x = -0.5
    r_y = res_ens['LowerTot']
    r_w = num_dots_ens
    r_h = res_ens['UpperTot'] - r_y
    sub[1].add_patch(Rectangle((r_x,r_y),r_w,r_h, facecolor='0.85'))   
    
    # #Add text with statistics
    # # place a text box in upper left in axes coords
    # textstr = '\n'.join((
    # r'$\mathrm{UpperCI}:%.2f$' % (res_ens['UpperCI'], ),
    # r'$\mathrm{UpperUnc}:%.2f$' % (res_ens['UpperUnc'], ),
    # r'$\mathrm{UpperTot}:%.2f$' % (res_ens['UpperTot'], )))
    # props = dict(alpha=0.5)
    # sub[1].text(0.05, 0.95, textstr, transform=sub[1].transAxes, fontsize=12, verticalalignment='top', bbox=props)


    # Set axis labels and ticks
    xlabels = sums.index
    ytitle = ini.FluxGas + ' ('+ini.ConvSums+')'
    sub[0].set_xticks(range(num_dots))
    sub[0].set_xticklabels(xlabels, rotation=45, ha='right')
    sub[0].set_xlabel('Gap filling techniques', fontsize=fsize, weight='bold', x=1)
    sub[0].set_ylim(y_min, y_max)
    sub[0].set_ylabel(ytitle, fontsize=fsize, weight='bold')
    sub[0].yaxis.grid(True, linestyle='-', which='major', color='gray', linewidth=0.2)
    
    # Set axis labels and ticks for ensemble sum
    xlabels = sums_ens.index
    sub[1].set_xticks(range(num_dots_ens))
    sub[1].set_xticklabels(xlabels, rotation=45, ha='right')
    sub[1].set_ylim(y_min_ens, y_max_ens)
 #   sub[1].yaxis.grid(True, linestyle='-', which='major', color='gray', linewidth=0.2)
    sub[1].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True, bottom=True, top=False, left=False, right=True)
    
#%%
    if suffix != '': suffix = '_' + suffix 
    mgf.utl.save_plot('pse', ini, file_path, run_number, suffix='double_real'+suffix)
