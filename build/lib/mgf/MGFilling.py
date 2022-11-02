#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py

Functions:
    # run_GFT():        Run gap-filling techniques (GFT), i.e. routines to fill gaps with multiple techniques
    # inspect_FF():     Inspect filled fluxes (FF), i.e. load gap-filled fluxes and (if available) modelled fluxes, summarize results and plot
    # bootstrap_FF():   Bootstrap scenarios from filled fluxes (FF) at data points of artificials gaps
    # analyse_BS():     Analyse bootstrap (BS) results including aggregating fluxes (sums) with uncertainties
    # pick_GE():        Pick good GFTs for GFT ensemble
    
"""
#%% Initialization
### Imports from python
import pandas as pd
import numpy as np
import datetime as dt #for timeit
import sys

### Own imports
sys.path.append('/Users/antje/W_Projects/MGF/_MultiGapFilling')
import mgf

### Options
pd.options.display.max_rows = 20

"""
# For reloading library (after changes):
    import importlib #for reloading libraries
    importlib.reload(mgf)
    In contrast to %autoreload, importlib.reload() also reset global variables set in module. In most cases, it is what you want.

# Automatic reloading of libraries
    #%load_ext autoreload
    %autoreload 2
    
# For TESTING
    ini_file='tNr_BaF_2016_v03.ini'; file_path='/Users/antje/W_Projects/MGF/_MultiGapFilling/example/DE-BaF_tNr';
    run_number='202208121237';
    good_GFTs='FDA_hh6|MDA_hh5|LUT_V1V2_d7|LUT_V1V2V3_d7|ANN_CRU|ANN_all'
    
"""

#%%    
def run_GFT(ini_file, file_path):
    # run_GFT():        Run gap-filling techniques, i.e. routines to fill gaps with multiple techniques
    #                   Note that filling the data points of real gaps within the 'hhs' scenarios means filling the real gaps
    #                   If not all gaps could be filled by a specific GFT, the remaining gaps are filled with default technique as specified in the ini-file
    # ini_file:         Name of ini file
    # file_path:        File path to data directory

    # Set run number
    now = dt.datetime.now()   
    run_number = str(now.year)+str(now.month).zfill(2)+str(now.day).zfill(2)+str(now.hour).zfill(2)+str(now.minute).zfill(2)
    
    # Start to log to file
    run_path = file_path + '/' + run_number
    mgf.utl.create_dir(run_path)
    stdout_org, log_file = mgf.utl.start_to_log(file_path, run_number)
    
    # Read in ini-file and data
    ini = mgf.utl.load_ini(ini_file, file_path)
        
    # Load and init data for gap filling
    data = mgf.utl.init_mgf_run(ini, file_path, run_number, mgf.__version__)
        
    # Set up scenario names
    scenarios=mgf.utl.name_scen()
    
    # Define LUT techniques
    LUT_MGFs = mgf.gft.def_LUTs(var_1=ini.LUTVar_1, range_1=float(ini.LUTRange_1), var_2=ini.LUTVar_2, range_2=float(ini.LUTRange_2), 
                     var_3=ini.LUTVar_3, range_3=float(ini.LUTRange_3), var_light=ini.LightVar, var_thres=float(ini.LightThres))
    LUT_MDS = mgf.gft.def_LUT_MDS(var_1=ini.LUT2Var_1, range_1=float(ini.LUT2Range_1), var_2=ini.LUT2Var_2, range_2=float(ini.LUT2Range_2), 
                     var_3=ini.LUT2Var_3, range_3=float(ini.LUT2Range_3))
    LUTs = pd.concat([LUT_MGFs,LUT_MDS], axis=1)

    techniques = LUTs.loc['Technique'].values.tolist() 
    print('*** The following techniques will be applied:', techniques + ['IP_lin','IP_mov'], '***')

    # Fill gaps with various techniques    
    flux_org = mgf.utl.name_flux(ini)
    t_start = dt.datetime.now()
    # Repeat for hhs and days scenarios
    #TESTING: mgf.gft.fill_LUT(data, flux_org, LUTs['LUT_V1_d3'], 'hhs') 
    for i_s in (scenarios):
    #TESTING: i_s='hhs'; i_s='days'
        mgf.gft.calc_IPs(data, flux_org, 'IP_lin', i_s)
        mgf.gft.calc_IPs(data, flux_org, 'IP_mov', i_s) 
        for i_t in (techniques):
        #TESTING: i_t='LUT_V1_d3';
            mgf.gft.fill_LUT(data, flux_org, LUTs[i_t], i_s) 
        print('Checks: GFT run time from ', t_start, ' to ', dt.datetime.now())
    
    mgf.utl.save_results(data, ini, file_path, run_number, prefix='../_mgf/mgf_Scenarios')
    print('*** Gap filling of scenarios with several techniques finished for', ini.FluxGas,', run number:',  run_number, '***')
    
    # End to log to file
    mgf.utl.end_to_log(stdout_org, log_file)

    return(run_number) #Name of gap filling run results

#%%    
def inspect_FF(file_path, run_number):
    # inspect_FF():     Inspect filled fluxes, i.e. load gap-filled fluxes and (if available) modelled fluxes, summarize results and plot
    # file_path:        File path to data directory
    # run_number:       Number of gap-filling run

    # Start to log to file
    stdout_org, log_file = mgf.utl.start_to_log(file_path, run_number)
    
    # Read in ini-file
    ini = mgf.utl.reload_ini(file_path=file_path, run_number=run_number)
    
    # Read in run_GFT results (fluxdata plus scenarios)
    run_path = file_path + '/' + run_number
    file_name = run_path+'/_mgf/mgf_Scenarios_' + ini.FluxGas + '_' + run_number + '.csv'        
    data = pd.read_csv(file_name, sep=',', encoding='latin1', parse_dates=['DateTime'], index_col=['DateTime']) 
    print('*** Loaded results of GFT scenarios for :', data.filter(regex="(_hhs)$").columns.values)
    
    # Read in gap filling results from models if available 
    if (ini.FileModels != 'none'):
        # Load model data
        data_mod = pd.read_csv(file_path + '/' + ini.FileModels, sep=',', encoding='latin1', parse_dates=['DateTime'], index_col=['DateTime'])
        check_index = data_mod.index.equals(data.index)
        print('Checks: Timestamps of model file equal time stamps of data file:', check_index)
        if (check_index == False):
            print('!!! Timestamp indices of flux data and model data are not equal:', data_mod.index.difference(data.index) )
            print('!!! Model data will not be loaded.')
        else:
            print('*** Loaded model results with columns:', data_mod.columns.values)
            mgf.utl.save_results(data_mod, ini, file_path, run_number, prefix='mgf_Models', folder='_mgf')
            data = pd.concat([data, data_mod], axis=1)
    else:
        print('*** No model results to load. ***')
    
    # Fill real data gaps with the various techniques using the 'hhs' scenarios on the real gaps.
    col_ToUse = data.filter(regex="(?<!_days)$").columns
    data_filled = data[col_ToUse.values].copy()
    techniques = data.filter(regex="(_hhs)$").columns.values
    flux_org = mgf.utl.name_flux(ini)
    for i_t in techniques:
        data_filled[i_t] = data_filled.apply(lambda row: row[i_t] if np.isnan(row[flux_org]) else row[flux_org], axis=1)
        data_filled[i_t] = data_filled.apply(lambda row: row[ini.DefGFT] if np.isnan(row[i_t]) else row[i_t], axis=1)
        
    #Save data frame with gap-filled real gaps                                             
    data_filled.rename(columns=lambda x: x.replace('_hhs','_real'), inplace=True)
    mgf.utl.save_results(data_filled, ini, file_path, run_number, prefix='mgf_FilledReal')
    
    #Save all data
    data_all = pd.concat([data, data_filled.filter(regex="(_real)$")], axis=1)
    mgf.utl.save_results(data_all, ini, file_path, run_number, prefix='mgf_AllData', folder='_mgf')
    print('*** All data (fluxes, models, scenarios, filled real gaps) saved ***')
    
    #*** The next steps could be in a seperate routine loading all data *** 
    # Make description of data_allset
    descr = mgf.ffa.make_descr(data_all, ini)
    mgf.utl.save_results(descr, ini, file_path=file_path, run_number=run_number, prefix='descr', suffix='ini')
    
    # Gap distribution plots
    mgf.pgd.plot_gaps(data_all, ini, flux_org, file_path, run_number)
    
    # Scatterplots
    mgf.psc.plot_scatter(data_all, ini, 'hhs', file_path, run_number)
    mgf.psc.plot_scatter(data_all, ini, 'days', file_path, run_number)
    
    # Residuen scatterplots
    mgf.psc.plot_scatter(data_all, ini, 'hhs', file_path, run_number, suffix='resid')
    mgf.psc.plot_scatter(data_all, ini, 'days', file_path, run_number, suffix='resid')
        
    # Scatterplot subset of data_all e.g. for manuscript
    if mgf.__TEST__ == True: # (run only manually when testing)
        data_all_sub = data_all.filter(regex=(flux_org+'|IP_lin|IP_mov|FDA_hh6|MDA_hh5|MDC_d3|LUT_V1_d3|LUT_V1V2_d3|ANN'))
        mgf.psc.plot_scatter(data_all_sub, ini=ini, scenario='hhs', file_path=file_path, run_number=run_number, suffix='sub')

    # Time series plots (fluxes, residuals, cum sums...)
    mgf.pts.plot_series(data_all, ini, 'hhs',  file_path, run_number, suffix='fluxes_cplt')
    mgf.pts.plot_series(data_all, ini, 'days',  file_path, run_number, suffix='fluxes_cplt')
    mgf.pts.plot_series(data_all, ini, 'hhs',  file_path, run_number, suffix='fluxes_real') #Not applicable for 'days' since dealing with real gaps...
    mgf.pts.plot_series(data_all, ini, 'hhs',  file_path, run_number, suffix='cums_cplt')
    mgf.pts.plot_series(data_all, ini, 'days',  file_path, run_number, suffix='cums_cplt')
    if mgf.__TEST__ == True: #Might be confusing and the other plots are sufficient
        mgf.pts.plot_series(data_all, ini, 'hhs',  file_path, run_number, suffix='fluxes_all') 
        mgf.pts.plot_series(data_all, ini, 'days',  file_path, run_number, suffix='fluxes_all')
    # Time series plot of daily sums
    mgf.pds.plot_daily(data_all, ini, 'real',  file_path, run_number, suffix='')

    # End to log to file
    mgf.utl.end_to_log(stdout_org, log_file)
    
#%%
def bootstrap_FF(file_path, run_number):
    # bootstrap_FF():   Bootstrap scenarios from filled fluxes (at data points of artificials gaps)
    # file_path:        File path to data directory
    # run_number:       Number of gap-filling run
    
    # Start to log to file
    stdout_org, log_file = mgf.utl.start_to_log(file_path, run_number)
    
    # Read in ini-file and data
    ini = mgf.utl.reload_ini(file_path=file_path, run_number=run_number)
    data = mgf.utl.load_data(ini, file_path, run_number, prefix='AllData')
    flux_org = mgf.utl.name_flux(ini)

    # Bootstrap filled artificial gaps
    # No filter, i.e. full-time
    print('>>> Bootstrapping full-time:')
    boot_hhs_ft = mgf.ffa.bootstrap_arti(data=data, col_obs=flux_org, scenario='hhs')
    boot_days_ft = mgf.ffa.bootstrap_arti(data=data, col_obs=flux_org, scenario='days')
    mgf.utl.save_results(boot_hhs_ft, ini, file_path, run_number, prefix='boot', suffix='ft_hhs', folder='_mgf')
    mgf.utl.save_results(boot_days_ft, ini, file_path, run_number, prefix='boot', suffix='ft_days', folder='_mgf')
       
    # Filter day-time
    print('>>> Bootstrapping day-time:')
    data_dt = mgf.utl.set_time_of_day(data, ini, 'dt')
    boot_hhs_dt = mgf.ffa.bootstrap_arti(data=data_dt, col_obs=flux_org, scenario='hhs')
    boot_days_dt = mgf.ffa.bootstrap_arti(data=data_dt, col_obs=flux_org, scenario='days')
    mgf.utl.save_results(boot_hhs_dt, ini, file_path, run_number, prefix='boot', suffix='dt_hhs', folder='_mgf')
    mgf.utl.save_results(boot_days_dt, ini, file_path, run_number, prefix='boot', suffix='dt_days', folder='_mgf')
        
    # Filter night-time by setting observed value to gap (hence not avialable for bootstrapping)
    print('>>> Bootstrapping night-time:')
    data_nt = mgf.utl.set_time_of_day(data, ini, 'nt')
    boot_hhs_nt = mgf.ffa.bootstrap_arti(data=data_nt, col_obs=flux_org, scenario='hhs')
    boot_days_nt = mgf.ffa.bootstrap_arti(data=data_nt, col_obs=flux_org, scenario='days')
    mgf.utl.save_results(boot_hhs_nt, ini, file_path, run_number, prefix='boot', suffix='nt_hhs', folder='_mgf')
    mgf.utl.save_results(boot_days_nt, ini, file_path, run_number, prefix='boot', suffix='nt_days', folder='_mgf')
    
    # End to log to file
    mgf.utl.end_to_log(stdout_org, log_file)
    
#%%
def analyse_BS(file_path, run_number):
    # analyse_BS():     Analyse bootstrap results including aggregating fluxes (sums) with uncertainties
    # file_path:        File path to data directory
    # run_number:       Number of gap-filling run
    
    # Start to log to file
    stdout_org, log_file = mgf.utl.start_to_log(file_path, run_number)

    # Read in ini-file
    ini = mgf.utl.reload_ini(file_path=file_path, run_number=run_number)
    
    # Load boot data
    # Full-time
    boot_hhs_ft = mgf.utl.load_bootres(ini, file_path, run_number, suffix='ft_hhs')
    boot_days_ft = mgf.utl.load_bootres(ini, file_path, run_number, suffix='ft_days')
    # Day-time
    boot_hhs_dt = mgf.utl.load_bootres(ini, file_path, run_number, suffix='dt_hhs')
    boot_days_dt = mgf.utl.load_bootres(ini, file_path, run_number, suffix='dt_days')
    # Night-time
    boot_hhs_nt = mgf.utl.load_bootres(ini, file_path, run_number, suffix='nt_hhs')
    boot_days_nt = mgf.utl.load_bootres(ini, file_path, run_number, suffix='nt_days')


   # Plot bootstrapping statistics for fulltime as double plot with hhs & days in one plot
    mgf.pbs.plot_bootstats_double(boot_hhs_ft, boot_days_ft, ini=ini, scens=['hhs','days'], times=['full-time','full-time'], colors='gft_rainbow', file_path=file_path, run_number=run_number)
    # Plot bootstrapping statistics for ft, dt, nt in one plot of all techniques
    mgf.pbs.plot_bootstats_triple(boot_hhs_ft, boot_hhs_dt, boot_hhs_nt, ini=ini, scens=['hhs','hhs','hhs'], times=['full-time','day-time','night-time'], colors='gft_rainbow', file_path=file_path, run_number=run_number)
    mgf.pbs.plot_bootstats_triple(boot_days_ft, boot_days_dt, boot_days_nt, ini=ini, scens=['days','days','days'], times=['full-time','day-time','night-time'], colors='gft_rainbow', file_path=file_path, run_number=run_number)

    # Plot subset of data and special plots e.g. for manuscript
    if mgf.__TEST__ == True: # (run only manually when testing)
        boot_days_ft_sub = boot_days_ft.filter(regex=('IP|MDA_hh5|MDC_d3|MDC_d7|ANN'))
        boot_hhs_ft_sub = boot_hhs_ft.filter(regex=('IP|MDA_hh5|MDC_d3|MDC_d7|ANN'))
        boot_days_dt_sub = boot_days_dt.filter(regex=('IP|MDA_hh5|MDC_d3|MDC_d7|ANN'))
        boot_days_nt_sub = boot_days_nt.filter(regex=('IP|MDA_hh5|MDC_d3|MDC_d7|ANN'))
        # Triple box stats plot
        mgf.pbs.plot_bootstats_triple(boot_days_ft_sub, boot_days_dt_sub, boot_days_nt_sub, ini=ini, scens=['days','days','days'], times=['full-time','day-time','night-time'], colors='gft_gray', file_path=file_path, run_number=run_number, suffix='sub_g')
        # Double box stats plot
        mgf.pbs.plot_bootstats_double(boot_hhs_ft_sub, boot_days_ft_sub, ini=ini, scens=['hhs','days'], times=['full-time','full-time'], colors='gft_gray', file_path=file_path, run_number=run_number, suffix='sub_g')
        # Single plots with different color settings
        mgf.pbs.plot_bootstats_single(boot_days_ft_sub, ini=ini, scens=['days'], times=['full-time'], colors='gft_rainbow', file_path=file_path, run_number=run_number, suffix='sub')       
        mgf.pbs.plot_bootstats_single(boot_days_nt_sub, ini=ini, scens=['days'], times=['night-time'], colors='gft_gray', file_path=file_path, run_number=run_number, suffix='sub_g')
        mgf.pbs.plot_bootstats_single(boot_days_dt_sub, ini=ini, scens=['days'], times=['day-time'], colors=['limegreen', 'teal'], file_path=file_path, run_number=run_number, suffix='sub_2col')

    # Calculate and save bootstrap errors
    err_ft = mgf.ffa.calc_errors(boot_hhs_ft,boot_days_ft)
    mgf.utl.save_results(err_ft, ini, file_path, run_number, prefix='berr', suffix='ft', folder='_mgf') #be = bootstrap errors
    err_dt = mgf.ffa.calc_errors(boot_hhs_dt,boot_days_dt)
    mgf.utl.save_results(err_dt, ini, file_path, run_number, prefix='berr', suffix='dt', folder='_mgf')
    err_nt = mgf.ffa.calc_errors(boot_hhs_nt,boot_days_nt)
    mgf.utl.save_results(err_nt, ini, file_path, run_number, prefix='berr', suffix='nt', folder='_mgf')
 
    # Table with summarized results
    err_table = mgf.ffa.make_table(ini, err_ft, err_dt, err_nt)
    mgf.utl.save_results(err_table, ini, file_path, run_number, prefix='errors')

    # Calculate flux sums, save and plot
    data = mgf.utl.load_data(ini, file_path, run_number, prefix='AllData')
    sums = mgf.ffa.calc_sums(data, ini, err_ft=err_ft)
    mgf.utl.save_results(sums, ini, file_path, run_number, prefix='sums', suffix='all')
    mgf.utl.show_table(sums)
    # Plot sums in figure and table
    mgf.pse.plot_sums(sums, ini, colors='gft_rainbow', file_path=file_path, run_number=run_number)
    if mgf.__TEST__ == True: # As gray color plot (run only manually when testing)
        mgf.pse.plot_sums(sums, ini, colors='gft_gray', file_path=file_path, run_number=run_number, suffix='g')
    mgf.pta.plot_table(sums, ini, file_path=file_path, run_number=run_number, suffix='all')
           
    # End to log to file
    mgf.utl.end_to_log(stdout_org, log_file)

#%%    
def pick_GE(file_path, run_number, good_GFTs=''):
    # pick_GE():        Pick good GFTs for GFT ensemble
    # file_path:        File path to data directory
    # run_number:       Number of gap-filling run
    # good_GFTs:        List of GFT names separated with | in a string for regex, i.e. 'MDA_hh5|MDC_d3|MDC_d7'
    
    # Start to log to file
    stdout_org, log_file = mgf.utl.start_to_log(file_path, run_number)

    # Read in ini-file
    ini = mgf.utl.reload_ini(file_path=file_path, run_number=run_number)
    
    # Set picks for manuscript
    if good_GFTs == '': #not provided, then pick GFTs of manuscript
        if (ini.FluxGas=='NH3'): good_GFTs = 'IP|WDM|FDA|MDA|LUT_V1_d7|LUT_V1V2_d7|LUT_V1V2V3_d7|ANN' 
        
    # Calculate flux sums, save and plot for subset
    # Load and filter data
    data = mgf.utl.load_data(ini, file_path, run_number, prefix='AllData')
    data_ens = data.filter(regex=(mgf.utl.name_flux(ini)+'|'+good_GFTs))
    err_ft = mgf.utl.load_err(ini, file_path, run_number, suffix='ft')
    # Calculate sums
    sums = mgf.ffa.calc_sums(data, ini, err_ft=err_ft)
    sums_ens = mgf.ffa.calc_sums(data_ens, ini, err_ft=err_ft)
    mgf.utl.save_results(sums_ens, ini, file_path, run_number, prefix='sums', suffix='ens')
    mgf.utl.show_table(sums_ens)
    mgf.pta.plot_table(sums_ens, ini, file_path=file_path, run_number=run_number, suffix='ens')
    res_ens = mgf.ffa.calc_ensemble(sums_ens, ini)
    mgf.utl.save_results(res_ens, ini, file_path, run_number, prefix='res', suffix='ens')
    
    # Plot
    mgf.pse.plot_sums(sums_ens, ini, colors='gft_rainbow', file_path=file_path, run_number=run_number, suffix='ens')
    mgf.pse.plot_sums_ens(sums, sums_ens, res_ens, ini, colors='gft_rainbow', file_path=file_path, run_number=run_number, suffix='ens')
    mgf.pds.plot_daily(data_ens, ini, 'real',  file_path, run_number, suffix='ens')
    # Plots with special formats for manuscript
    if mgf.__MS__ == True:
        mgf.pse.plot_sums_ens(sums, sums_ens, res_ens, ini, colors='gft_rainbow', file_path=file_path, run_number=run_number, suffix='ens_ms')
        mgf.pds.plot_daily(data_ens, ini, 'real',  file_path, run_number, suffix='ens_ms')
    print('*** Ensemble results of the aggregated fluxes:', res_ens)
    
    # End to log to file
    mgf.utl.end_to_log(stdout_org, log_file)
 
