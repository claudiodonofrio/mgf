#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py
    
Functions:
    Many small utility functions for globally used settings, data handling, file handling, plot colors, ...,
    for a detailed description see below.

"""
#%% Initialization
### Imports from python
import pandas as pd
import numpy as np
from decimal import Decimal #For correct rounding in OS X
import os
import sys
import shutil
import matplotlib.pyplot as plt

### Options
pd.options.display.max_rows = 20

#%% Settings/Naming of globally used variables

def name_flux(ini):
    # name_flux():      Set name of original flux to be filled
    # ini:              Settings from ini file
    flux_org = 'FluX_' + ini.FluxGas #Capital X to hopefully have unique name
    return flux_org
    
def name_run(ini, run_number)    :
    # name_run():       Set name of gap fillin run
    # ini:              Settings from ini file
    # run_number:       Number of gap filling run    
    return ini.FluxGas + '_' + run_number

def name_scen():
    # name_scen():      Return fixed names of scenarios
    scenarios = ['hhs','days'] #Artificial gap scenarios: single half-hours and single full days
    return scenarios

def set_boot():
    # set_boot():       Set bootstrap repetitions
    rep = 999; #number of repititions, 99 for testing 999 for final run
    # Percent of the data to be used for bootstrapping
    percent = 50; #compromise between 0 and 100 and typical percentage of gaps
    seed = 99; #seed for random number generator
    return rep, percent, seed

def set_gapl_S(): 
    # set_gapl_S():     Set length of short gap to distinguish short gaps
    # Used for gap filling with linear interpolation and for calculating uncertainties
    gapl_S = 12; # 12 hhs = 6 hours
    return gapl_S

def set_time_of_day(data, ini, time_of_day):
    # set_time_of_day(): Keep fluxes of time of day and set the rest to NaN using daylight information
    # data:             Dataframe with timestamp index, flux measurements, (maybe meteo,) and filled fluxes for various techniques
    # ini:              Settings from ini file
    # time_of_day:      'dt' or 'nt'  
    data_copy = data.copy()
    flux_org = name_flux(ini)
    if (time_of_day == 'dt'):
        data_copy.loc[data_copy[ini.LightVar] <= float(ini.LightThres), flux_org] = np.nan
    elif (time_of_day == 'nt'):
        data_copy.loc[data_copy[ini.LightVar] > float(ini.LightThres), flux_org] = np.nan
    else:
        print('!* String for time of day not specified correctly:', time_of_day)
    return data_copy
    
#%% Avoiding rounding problems

def round_dec(number, digits):
    # round_dec():      Round numbers explicitly
    # number:           Number to round
    # digits:           Number of digits
    r_number = float(round(Decimal(number),digits))
    # OS X: Some times rounding errors like 26.001999999999999 instead of 26.002 (for round() as well as np.round() )
    # Hard to reproduce since occurs within variables and not if numbers are provided explicitly
    # For other notes on rounding, see here: https://stackoverflow.com/questions/13479163/round-float-to-x-decimals
    return r_number
    
#%% Count gap distribution

def countGaps(data, col_name):
    # countGaps():      Return vector with gap length counter
    # data:             Pandas data frame
    # col_name:         Name of column with gaps (NAs) to count
    #
    # Cumulative sum of gaps
    cum_sum =  data[col_name].isnull().cumsum()
    # Series of cumulative maxima per non-gap
    cum_max = (data[col_name].notnull()*cum_sum).cummax()
    # Gap length counter
    count_gaps = cum_sum - cum_max
    return count_gaps

def countGapsTot(data, col_name):
    # countGLength():   Return vector with total gap length
    # data:             Pandas data frame
    # col_name:         Name of column with gaps (NAs) to count
    #    
    # Group by gaps
    cum_group = (data[col_name].notnull() != data[col_name].notnull().shift(1)).cumsum() * data[col_name].isnull()
    # Gap total length
    count_length = data.groupby(cum_group)[col_name].transform('size') * data[col_name].isnull()
    return count_length

    #Notes: Code inspired from:
    #https://stackoverflow.com/questions/32850185/change-value-if-consecutive-number-of-certain-condition-is-achieved-in-pandas
    # .shift(1) compares current row's condition to next row's which has effect of marking 'True' where a run of values begins.

def distrGaps(data, col_name):
    # distrGaps():      Return data frame how gaps are distributed (length, frequency, total number of hhs)
    # data:             Pandas data frame
    # col_name:         Name of column with gaps (NAs) to count
    #  
    gaplength = countGapsTot(data, col_name)
    gaps = pd.DataFrame(gaplength.value_counts()).sort_index()
    gaps['length'] = gaps.index
    gaps['freq'] = gaps[col_name]/gaps['length']
    gaps['amount'] = gaps['length'] * gaps['freq'] 
    gaps['sums'] = gaps[col_name].cumsum() - gaps.iloc[0,0]
    gaps['perc_gaps'] = gaps[col_name]/sum(gaps.iloc[1:,0])*100 #percentage of gaps from all gaps
    gaps['perc_gaps_sums'] = gaps['sums']/sum(gaps.iloc[1:,0])*100
    gaps['perc_data'] = gaps[col_name]/data[col_name].size*100 #percentage of gaps in data
    gaps['perc_data_sums'] = gaps['sums']/data[col_name].size*100
    return gaps 

#%% Data handling

def gen_fcol(data, ini, col_name):
    # gen_fcol():       Generate flux column with original data (with or w/o flagging) to be filled
    # data:             Dataframe with timestamp index, flux measurements, (maybe meteo,) and filled fluxes for various techniques
    # ini:              Settings from ini file

    if ini.FluxFlag == 'none': # No flag provided
        # Copy and if -9999, set to nan
        data.insert(loc=0, column=col_name, value=np.where(data[ini.FluxColumn] == -9999, np.nan, data[ini.FluxColumn]))
        print("New column created with fluxes to fill for", col_name, 'with no flagging')
    else:
        # Only use these fluxes up to FlagMax quality, set rest to nan
        data.insert(loc=0, column=col_name, value=np.where(data[ini.FluxFlag] <= float(ini.FlagMax), data[ini.FluxColumn], np.nan))
        print("New column created with fluxes to fill for", col_name, 'with flagging (', ini.FluxFlag, '<=', ini.FlagMax, ')')


def check_completeness(data, scenario):
    # check_completeness(): Function that checks if the artificial gaps are filled with the complete suite of techniques
    #                       (at all positions, i.e. measurements and/or real gaps)
    # data:                 Dataframe with timestamp index, flux measurements, (maybe meteo,) and filled fluxes for various techniques
    # scenario:             Artificial gap scenario
    #
    #Returns a boolean vector
    return data.filter(regex=scenario).notnull().sum(axis=1) == data.filter(regex=scenario).shape[1]

#%% Logging to screen and file

class Logger(object):    
    def __init__(self,log_file):
   #def __init__(self, log_file):
        self.terminal = sys.stdout
        #self.log = open('/Users/amoffat/Documents/MGF/_MGF_Runs/DE-BaF_NrX//201910171038/log_201910171038.log','a')
        self.log = log_file

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        #pass    
        
def start_to_log(file_path, run_number):
    # start_to_log():   Start logging to screen and file
    # file_path:        File path to data directory
    # run_number:       Number of gap-filling run
    
    stdout_org = sys.stdout
    # Log to screen and file
    run_path = file_path + '/' + run_number
    log_file = open(run_path+'/log_'+run_number+'.log', 'a')
    sys.stdout = Logger(log_file)
    # For only logging to log file: sys.stdout = log_file

    return stdout_org, log_file

def end_to_log(stdout_org, log_file):
    # end_to_log():   End logging to screen and file
    # stdout_org and log_file from start_to_log()
    
    sys.stdout = stdout_org
    log_file.close()
    
#Inspired from here:
#        https://stackoverflow.com/questions/14906764/how-to-redirect-stdout-to-both-file-and-console-with-scripting
#        https://stackoverflow.com/questions/7152762/how-to-redirect-print-output-to-a-file-using-python
        
#%% File handling (loading and saving)
    
def create_dir(path):
    # create_dir():     Create directory
    # path:             Path and name of directory
    try: #make run directory
        os.mkdir(path)
    except OSError:
        print("Could not create dir:", path)
    else:
        print ("Directory created:", path)

def init_mgf_run(ini, file_path, run_number, pkg_vers):
    # init_mgf_run():   Initialize MGF run
    # ini:              Settings from ini file
    # file_path:        File path to data directory
    # pkg_vers:         Provide version of MGF package #TESTING: pkg_vers=mgf.__version__

    # Load flux data to be filled
    data = pd.read_csv(str(file_path+'/'+ini.FileData), sep=ini.FileSeparator, encoding='latin1', parse_dates=['DateTime'], index_col=['DateTime'])
    print('*** Fill gaps for:', ini.FileData, 'with', data.shape[0], 'data rows which equals', data.shape[0] / 48, 'days. ***' )

    # Check if time stamps are unique
    print('Checks: Time stamps are unique:', data.index.is_unique)
    
    # Check if full days
    print('Checks: Length of time stamps are a multiple of complete days:', (data.shape[0] % 48 == 0))
    
    # Check if first / last half-hour
    last = data.shape[0]-1
    check_firsthalf = (data.index[0].hour == 0) & (data.index[0].minute == 30)
    print('Checks: First day starts at first half-hour (Fluxnet convention 0:30):', check_firsthalf)
    check_lasthalf = (data.index[last].hour == 0) & (data.index[last].minute == 0)
    print('Checks: Last day ends at mid night (Fluxnet convention 0:00 next day):', check_lasthalf)   
    
    # Check for complete days
    index = pd.date_range(start=data.index[0], end=data.index[last], freq='30min')
    check_complete = data.index.equals(index)
    print('Checks: Time stamps consist of continuous half-hours:', check_complete)
    
    # Generate best fluxes to be filled column
    flux_org = name_flux(ini)
    check_unique = (flux_org not in data.columns)
    print('Checks: Name of new to be filled column', flux_org, 'is unique:', check_unique)
    gen_fcol(data, ini, flux_org)
    
    # Count gaps in data, 
    gaps = countGapsTot(data, flux_org)
    check_maxsize = max(gaps) <= 10*48
    print('Checks: Longest gap in fluxes is smaller than ten days:', check_maxsize, '(Max gap size:', max(gaps), 'hhs)')
    
    if (sum([check_firsthalf, check_lasthalf, check_complete, check_unique, check_maxsize]) < 5):
        sys.exit("!!! At least one of the file checks failed, see log for more information.")
    
    # Create sub directories
    run_path = file_path + '/' + run_number
    create_dir(run_path+'/_in')
    create_dir(run_path+'/_mgf')
    create_dir(run_path+'/_res')
    create_dir(run_path+'/_plots')
    
    #Copy original files to /_in
    shutil.copy2(file_path+'/'+ini.FileData, file_path+'/'+run_number+'/_in/')
    shutil.copy2(file_path+'/'+ini.FileIni, file_path+'/'+run_number+'/_in/')
    if (ini.FileModels != 'none'):
        shutil.copy2(file_path+'/'+ini.FileModels, file_path+'/'+run_number+'/_in/')
    print('*** Original data file copied to', './'+run_number+'/_in/'+ini.FileData, '***')
    
    #Write ini with additional run information
    ini['RunNumber'] = run_number
    ini['CodeVersion'] = pkg_vers
    ini.to_csv(file_path+'/'+run_number+'/_mgf/ini_'+run_number+'.ini', float_format='%g', na_rep='NA', header=True)
    print('*** Ini settings loaded and saved to', './'+run_number+'/_mgf/ini_'+run_number+'.ini', '***')
    
    # Last checks
    print('Checks: Number of gaps in all columns:')
    print(" ".join(data.isnull().sum().to_string().replace('\n',', ').split()))
    file_name = save_results(data, ini, file_path, run_number, prefix='../_mgf/gf_FluxData')
    print('*** Original flux data loaded and saved to: ', file_name, '***')
        
    # Return data and run number (no modeling results since these only needed to be copied)
    return data

def load_ini(ini_file, file_path):
    # load_ini():       Load original ini with settings for new run
    # ini_file:         Name of ini file
    # file_path:        File path to data directory
    
    # Read in ini-file 
    ini_data = pd.read_csv(file_path+'/'+ini_file, encoding='latin1', index_col=['Variable'])
    ini = ini_data.Settings
    print('*** Loaded ini:', file_path+'/'+ini_file)
    
    # Add information of file name
    ini['FileIni'] = ini_file #Original ini file

    #Returns
    return ini

def reload_ini(file_path, run_number):
    # reload_ini():     Load ini with settings for existing run
    # ini_file:         Name of ini file
    # file_path:        File path to data directory

    # Read in ini-file 
    ini_file = './'+run_number+'/_mgf/ini_'+run_number+'.ini'
    ini_data = pd.read_csv(file_path+'/'+ini_file, encoding='latin1', index_col=['Variable'])
    
    ini = ini_data.Settings
    print('*** Loaded ini:', file_path+'/'+ini_file)
    
    #Returns
    return ini

def save_results(results, ini, file_path, run_number='', prefix='', suffix='', folder='_res'):
    # save_results():   Save results to .csv file
    # results:          Some dataframe with results to be saved
    # prefix:           Prefix to filename
    # suffix:           Suffix to filename
    # file_path:        File path to data directory
    # run_name:         Name of gap filling run plus <optional> some additional information

    # Set prefix / suffix
    if prefix != '':
        prefix = prefix+'_'
    if suffix != '':
        suffix = '_'+suffix

    # Set run name
    run_name = prefix + ini.FluxGas + '_' + run_number + suffix
    
    # Shorten run_name e.g. from plot_bootstats()
    run_name = run_name.replace('full-time','ft')
    run_name = run_name.replace('night-time','nt')
    run_name = run_name.replace('day-time','dt')
    
    # Save csv file
    run_path = file_path + '/' + run_number
    file_name = folder + '/' + run_name+  '.csv'
    
    if prefix == 'boot_': #double-index-structure pandas
        results.to_csv(run_path + '/' + file_name, float_format='%g', na_rep='NA',index = [0,1])
    else: #Simple data array
        results.to_csv(run_path + '/' + file_name, float_format='%g', na_rep='NA', index=True)
    
    print('*** Results saved to', file_name)
    return file_name

def load_data(ini, file_path, run_number, prefix):
    # load_data():  Read in data (e.g. fluxes, meteo, models, scenarios, filled real gaps)
    # See save_results() for a variable description
    run_path = file_path + '/' + run_number
    file_name = run_path + '/_mgf/mgf_' + prefix + '_' + ini.FluxGas + '_' + run_number + '.csv'        
    data = pd.read_csv(file_name, sep=',', encoding='latin1', parse_dates=['DateTime'], index_col=['DateTime'])
    print('*** Loaded all data (fluxes, meteo, models, scenarios, filled real gaps) for run', ini.FluxGas, run_number)
    
    #Returns
    return data

def load_err(ini, file_path, run_number, suffix):
    # load_err():   Read in error results
    # See save_results() for a variable description
    run_path = file_path + '/' + run_number
    file_name = run_path + '/_mgf/berr_' + ini.FluxGas + '_' + run_number + '_' + suffix + '.csv'        
    data = pd.read_csv(file_name, sep=',', encoding='latin1', index_col='Technique')
    print('*** Loaded bootstrap error estimates for run', ini.FluxGas, run_number)
    
    #Returns
    return data

def load_bootres(ini, file_path, run_number, suffix):
    # load_bootres():   Read in bootstrapping results
    # See save_results() for a variable description
    run_path = file_path + '/' + run_number
    file_name = run_path+'/_mgf/boot_'+ ini.FluxGas + '_' + run_number + '_' + suffix + '.csv'      
    data = pd.read_csv(file_name, sep=',', encoding='latin1', index_col=[0,1])
    print('*** Loaded bootstrap results of', ini.FluxGas, run_number, suffix)
    
    #Returns
    return data
    
def show_table(results):
    # show_table():     Plot results in a table
    # results:          Dataframe with results

    fig, ax = plt.subplots()
    # Hide axes
    ax.axis('off')
    ax.axis('tight')

    table = ax.table(cellText=results.astype('float64').round(2).values, colLabels=results.columns, loc='center')
    table.scale(1.3,1.3)
    fig.tight_layout()

    plt.show()

def save_plot(plt_type, ini, file_path, run_number, suffix='', folder='_plots'):
    # save_plot():      Save plot in .pdf file format
    # file_path:        File path to data directory
    # run_name:         Name of gap filling run plus <optional> some additional information
    # suffix:           Suffix to filename
    # folder:           Name of plot directory
    
    if file_path != '':
        if suffix != '':
           suffix = '_'+suffix
           
        # Set run name
        run_name = plt_type + '_' + ini.FluxGas + '_' + run_number + suffix
        
        # Shorten run_name
        run_name = run_name.replace('full-time','ft')
        run_name = run_name.replace('night-time','nt')
        run_name = run_name.replace('day-time','dt')
        
        # Save plot to pdf
        run_path = file_path + '/' + run_number
        file_name = folder +'/' + run_name + '.pdf'
        plt.savefig(run_path+'/'+file_name, format='pdf', bbox_inches='tight')
        plt.show()
        print('Plot saved to:', file_name)
    else:
        plt.show() #Only show plot
    

#%%
def set_colors(colors, labels):
    # set_colors():      Set colors for techniques or user input
    # colors:            Setting or colors of user
    # labels:            Index with column names
    #
    #TESTING:
    #    colors=['limegreen', 'teal']; # Poster #['xkcd:light blue', 'xkcd:cerulean']
    #    colors='gft_rainbow'
    #    colors='gft_gray'
    #    colors=['whitesmoke'] #['silver'] #['gray']
    
    techLabels = labels.to_series().str.split('_', n = 1, expand=True)[0]
    boxColors = []
    techColumns = {}
    if (colors == 'gft_rainbow'): # Color by techniques: rainbow
        NUM_COLORS = 20
        cm = plt.get_cmap('gist_rainbow')
        color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
        
        techColumns['IP'] = color[1]
        techColumns['WDM'] = color[3]
        techColumns['FDA'] = color[5]
        techColumns['MDA'] = color[8]
        techColumns['MDC'] = color[11]
        techColumns['LUT'] = color[12]
        techColumns['ANN'] = color[17]
        techColumns['Model'] = color[19]
        techColumns['Default'] = color[19]
        for key in techLabels:
            boxColors.append(techColumns.get(key, techColumns['Default']))
    elif (colors == 'gft_gray'): # Color by techniques: gray
        techColumns['IP'] = '1.0'
        techColumns['WDM'] = '0.9'
        techColumns['MDA'] = techColumns['WDM']
        techColumns['FDA'] = techColumns['WDM']
        techColumns['MDC'] = techColumns['WDM']
        techColumns['LUT'] = '0.7'
        techColumns['ANN'] = '0.4'
        techColumns['Model'] = techColumns['ANN']
        techColumns['Default'] = techColumns['ANN']
        for key in techLabels:
            boxColors.append(techColumns.get(key, techColumns['Default']))
    else: #Colors provided by user
        boxColors = colors
        
    # Check which keys where used
    usedKeys=pd.unique(techLabels)
            
    return boxColors, techColumns, usedKeys

"""
# Draw rainbow colors
#https://stackoverflow.com/questions/8389636/creating-over-20-unique-legend-colors-using-matplotlib/8391452#8391452   
NUM_COLORS = 20

cm = plt.get_cmap('gist_rainbow')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
for i in range(NUM_COLORS):
    ax.plot(np.arange(10)*(i+1))

#fig.savefig('moreColors.png')
plt.show()
"""