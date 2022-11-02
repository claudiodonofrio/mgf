#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py

Functions:
    # plot_gaps():      Plot the distribution of gaps in the dataset
    
"""

#%% Initialization
### Imports from python
import matplotlib.pyplot as plt
import numpy as np

### Own imports
import mgf

"""
For TESTING:
    col_name = 'Flux_' + ini.FluxGas
    
"""

#%%
def plot_gaps(data, ini, col_name, file_path='', run_number=''):
    # plot_gaps():      Plot the distribution of gaps in the dataset
    # data:             Dataframe with timestamp index and flux measurements with real gaps
    # ini:              Settings from ini file
    # col_name:         Column name with gappy data (meteo or original or filled fluxes)
    # file_path:        File path to data directory
    # run_number:       Number of gap filling run

#%%    
    df = mgf.utl.distrGaps(data, col_name).drop([0])
    numBars = df.shape[0]
    #vline_pos = sum(df['length'] <= 12)
    fig_width = numBars * 0.5 
    fig_height = 5.0
    fsize = 12
    bar_width = 0.35
    m_size=16
    rot_grad = 0
    if mgf.__MS__ == True and ini.FluxGas == 'CO2':
        fig_width = 12.5
        m_size = 8
        rot_grad = 90
            
    #y_min = min(min(df['freq']),0)
    #y_max = max(max(df['freq']),100)
    
    # Some extended lists for drawing the cumulative gap percentage
    L = list(np.arange(numBars)-bar_width)
    L.extend(list(np.arange(numBars)))
    L.extend(list(np.arange(numBars)+bar_width))
    L.sort()
    P = list(df['perc_data_sums'])
    P.extend(list(df['perc_data_sums']))
    P.extend(list(df['perc_data_sums']))
    P.sort()
    
    fig, plx1 = plt.subplots(figsize=(fig_width, fig_height))
    fig.patch.set_facecolor('white')
    plx2 = plx1.twinx()
    #plt.axvline(x=vline_pos-1+bar_width, color='grey', linewidth=plx1.spines['top'].get_linewidth()) 

    plx1.bar(np.arange(numBars)-bar_width/2, df['freq'], width=bar_width, color='gray')
    plx2.plot(L, P,  color='black', linewidth=1, linestyle='dotted')
    plx2.plot(np.arange(numBars), df['perc_data_sums'],  color='black', marker='_', markersize=m_size, linestyle='')
    plx2.bar(np.arange(numBars)+bar_width/2, df['perc_data'], width=bar_width, color='red')
    
    plx1.set_xticks(range(numBars))
    plx1.set_xticklabels(df['length'].astype(int), rotation = rot_grad)
    plx1.set_xlabel('Length of gaps (hh)', fontsize=fsize, weight='bold')
    plx1.set_ylabel('Gap frequency', fontsize=fsize, weight='bold')
    plx2.set_ylim(0, 105)
    plx2.set_ylabel('Gap amount (%)', fontsize=fsize, weight='bold')
    
    # Make description
    tot_data=1000
    tot_gaps=600
    l1 = 'Gap distribution for ' + ini.FluxGas
    l2 = '\n# of fluxes: \t' + '\t' + str(sum(data[col_name].notnull()))
    l3 = '\n# of gaps: \t' + '\t' + str(sum(data[col_name].isnull()))
    l4 = '\n# in total: \t' + '\t' + ' ' + str(data[col_name].size)
    textstr = (l1 + l2 + l3 + l4).expandtabs()

    #Latex code: '$\sigma=%.2f$'

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)

    # place a text box in upper left in axes coords
    plx1.text(0.35, 0.95, textstr, transform=plx1.transAxes, fontsize=fsize, verticalalignment='top', bbox=props)
#%%
    mgf.utl.save_plot('pgd', ini, file_path, run_number, suffix='real')