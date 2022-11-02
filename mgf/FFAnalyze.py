#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py

Functions:
    # make_descr():     Make description of data and results
    # bootstrap_arti(): Bootstrap artificial scenarios from filled (predicted) half-hours
    # calc_errors():    Calculate error measures from bootstrapping results
    # format_errors():  Format strings for make_table()
    # make_table():     Combine error measures in table (bootstrapping results table, e.g. for manuscript)
    # calc_sums():      Calculate full period sums with errors from bootstrapping
    # calc_ensemble():  Calculate ensemble results from sums
    
"""
#%% Initialization
### Imports from python
import numpy as np
import pandas as pd

### Own imports
import mgf

### Options
pd.options.display.max_rows = 20

#%% Make dataframe with various properties of results (description)

def make_descr(data, ini):   
    # make_descr():     Make description of data and result properties (flux column, length, number of gaps, ...)
    # data:             Dataframe with timestamp index, flux measurements, (maybe meteo,) and filled fluxes for various techniques
    # ini:              Settings from ini file

    descr = {}
    # Settings from ini file
    descr = ini.to_dict()
    # Properties of filled fluxes
    flux_org = mgf.utl.name_flux(ini)
    descr['FluxOrg'] = flux_org
    descr['NumDays'] = mgf.utl.round_dec(data.shape[0]/48,2)
    descr['NumHHs'] = data.shape[0]
    descr['NumMeas'] = data[flux_org].notnull().sum()

    # Properties gaps
    descr['NumGaps'] = data[flux_org].isnull().sum()
    # Count small and long  gaps
    #Check totals
    if (descr['NumHHs'] != (descr['NumMeas'] + descr['NumGaps'])):
        print('!* Error in sum data points')
    # Count number of short and long real gaps
    gaps = mgf.utl.countGapsTot(data, flux_org)
    descr['NumGaps_short']  = sum(np.where((gaps != 0) & (gaps <= mgf.utl.set_gapl_S()), True, False))
    descr['NumGaps_long'] = sum(np.where((gaps != 0) & (gaps > mgf.utl.set_gapl_S()), True, False))
    if (descr['NumGaps'] != (descr['NumGaps_short'] + descr['NumGaps_long'])):
        print('!* Error in sum  gaps')
    descr['PercGaps'] = mgf.utl.round_dec((descr['NumGaps']/descr['NumHHs']*100),1)
     
    # Check for rows filled by all techniques (used for bootstrapping)
    v_cplt_all_hhs = data.filter(regex='hhs').notnull().sum(axis=1) == data.filter(regex='hhs').shape[1]
    v_cplt_all_days = data.filter(regex='days').notnull().sum(axis=1) == data.filter(regex='days').shape[1]
    data_dt = mgf.utl.set_time_of_day(data, ini, 'dt')
    data_nt = mgf.utl.set_time_of_day(data, ini, 'nt')
    descr['NumBoot_hhs_ft'] = sum(data[flux_org].notnull() & v_cplt_all_hhs)
    descr['NumBoot_hhs_dt'] = sum(data_dt[flux_org].notnull() & v_cplt_all_hhs)
    descr['NumBoot_hhs_nt'] = sum(data_nt[flux_org].notnull() & v_cplt_all_hhs)
    descr['NumBoot_days_ft'] = sum(data[flux_org].notnull() & v_cplt_all_days)
    descr['NumBoot_days_dt'] = sum(data_dt [flux_org].notnull() & v_cplt_all_days)
    descr['NumBoot_days_nt'] = sum(data_nt[flux_org].notnull() & v_cplt_all_days)
    
    descr = pd.Series(list(descr.values()), index=descr.keys(), name='Value')
    descr.index.name = 'Property'
    
    return descr

#%% Bootstrap statistical measures

def bootstrap_arti(data, col_obs, scenario):
    # bootstrap_arti(): Bootstrap artificial scenarios from filled (predicted) half-hours
    # data:             Dataframe with timestamp index, observed, and predicted fluxes for various techniques and models (plus other columns)
    # col_obs:          Data column with observed data (e.g. the measured fluxes)
    # scenario:         Artificial gap scenario, options: 'hhs','days' (#'mix' removed in _v14, Nov 2021 
    #                   ('mix' removed in _v14 (Nov 2021), used filled fluxes 50% of 'hhs' and 50% of 'days')
    # 
    #FOR TESTING: col_obs=mgf.utl.name_flux(ini); scenario='hhs'; scenario='days';

    #Data preparation and subsets
    col_tech = data.filter(regex=scenario).columns #Filter results   
     #Vector for artificial gaps filled with complete set of techniques at positions of measurements only, i.e. artificial gaps
    v_cplt_arti = data[col_obs].notnull() &  mgf.utl.check_completeness(data, scenario)

    print('>>> Number of artificial gaps used for bootstrapping \'', scenario,'\':', sum(v_cplt_arti))
    print('*! Missing artificial gaps due to (partially) incomplete techniques:', sum(data[col_obs].notnull() & ~mgf.utl.check_completeness(data, scenario)))

    #Settings for bootstrapping
    rep, percent, seed = mgf.utl.set_boot()
    np.random.seed(seed) # Set seed for random number generator with fixed number from mgf.utl.set_boot()
    print('>>> Settings for bootstrapping: repetitions =', rep, ', percent of data =', percent, ', seed =', seed)
    
    #Make result dataframe    
    index = pd.MultiIndex.from_product([['Bias', 'SDev', 'R2'], list(range(0,rep))])
    res = pd.DataFrame(index=index, columns=col_tech) #columns=pd.Index(col_tech, name='Technique'))
    
    #First repetitions, then techniques so that same subset for each technique
    for i_k in range(0,rep):
        # TESTING: rep=9; i_k=0
        # Pick indices from whole data range (i.e. subsets of gaps/non-gaps might vary)
        # Bootstrapping is typically performed with(!) replacement
        if scenario == 'hhs':
            num_set = int(data.shape[0] * percent / 100)
            boot_hh_all = np.random.choice(data.shape[0],num_set,replace=True)  #whole data range, multiple picks possible
            idx_hhs_all = data.index[boot_hh_all]
            idx_hhs_arti = idx_hhs_all[idx_hhs_all.isin(data[v_cplt_arti].index)] #subsetted to artificial gaps only
            data_boot = data[v_cplt_arti].loc[idx_hhs_arti]
            
        elif scenario == 'days': #Reprogrammed for picking full days
            num_days = int(data.shape[0]/48) # Available full days (no modulo)
            # Alternative idea: Could be limited to days with at least xx number of data points
            num_set = int(data.shape[0]/48 * percent / 100)
            # Position of artificial gap days
            boot_days_all = np.random.choice(num_days,num_set,replace=True)
            # Pick half-hours of day (moved to middle of time stamp for correct days)
            i_d = boot_days_all[0] # Initialize time index with first day
            day = data.index[0] - pd.DateOffset(hours=0.25) + pd.DateOffset(days=int(i_d))
            idx_cor = pd.date_range(start=day.strftime('%Y-%m-%d 00:15:00'), end=day.strftime('%Y-%m-%d 23:45:00'), freq='30min')
            for i_d in boot_days_all[1:]:
                day = data.index[0] - pd.DateOffset(hours=0.25) + pd.DateOffset(days=int(i_d))
                idx_cor = idx_cor.append(pd.date_range(start=day.strftime('%Y-%m-%d 00:15:00'), end=day.strftime('%Y-%m-%d 23:45:00'), freq='30min'))
                idx_days_all = idx_cor + pd.DateOffset(hours=0.25)
                # Half the days are not necessarily half the half-hours!!!
            idx_days_arti = idx_days_all[idx_days_all.isin(data[v_cplt_arti].index)] #subsetted to artificial gaps only
            data_boot = data[v_cplt_arti].loc[idx_days_arti] 
        
        # Calculate performance
        for i_t in col_tech:
            predicted = data_boot[i_t]
            observed = data_boot[col_obs]
            res[i_t].loc['Bias',i_k] = 1 / predicted.size * np.sum(predicted - observed) 
            #SDev with Laplace:
            res[i_t].loc['SDev',i_k] = (2**0.5) / predicted.size * np.sum(abs(predicted - observed)) 
            #SDev with Gauss:
            # res[i_t].loc['SDev',i_k] = ( sum(((predicted - observed) - np.mean(predicted - observed))**2) / predicted.size )**(0.5)
            res[i_t].loc['R2',i_k] = (sum((predicted-np.mean(predicted))*(observed-np.mean(observed))))**2 / sum((predicted-np.mean(predicted))**2) / sum((observed-np.mean(observed))**2)
            
    fn = lambda x: '{:04.2f}'.format(x) #format number to string
            
    for i_t in col_tech:
    #TESTING: i_t = col_tech[5];
        print(i_t, rep, percent, ', Bias:', fn(np.mean(res[i_t].loc['Bias'])), fn(np.std(res[i_t].loc['Bias'])),
              ', SDev:', fn(np.mean(res[i_t].loc['SDev'])), fn(np.std(res[i_t].loc['SDev'])), 
              ', R2:', fn(np.mean(res[i_t].loc['R2'])), fn(np.std(res[i_t].loc['R2'])))
    
    return res #Dataframe with result

#%% Calculate flux errors from bootstrapping results and sum up fluxes with errors
    
def calc_errors(boot_hhs,boot_days):
    # calc_errors():    Calculate error measures from bootstrapping results
    # boot_hhs:         Name of dataframe with 'hhs' bootstrap results
    # boot_days:        Name of dataframe with 'days' bootstrap results
    # 
    #FOR TESTING: boot_hhs=boot_hhs_ft; boot_days=boot_days_ft;

    col_tech_hhs = boot_hhs.columns.tolist()
    col_tech_days = boot_days.columns.tolist()
    errors = pd.DataFrame(columns=['Bias_10','Bias_90','SDev','SDev_SD','R2','R2_SD'], index=(col_tech_hhs+col_tech_days))
    errors.index.name = 'Technique'
    
    fr = lambda x: mgf.utl.round_dec(x,6) #round number to six digits
        
    for i_t in col_tech_hhs:
        #FOR TESTING: i_t = col_tech[0]
        errors['Bias_10'][i_t] = fr(np.percentile(boot_hhs[i_t].loc['Bias'],10))
        errors['Bias_90'][i_t] = fr(np.percentile(boot_hhs[i_t].loc['Bias'],90))
        errors['SDev'][i_t] = fr(np.mean(boot_hhs[i_t].loc['SDev']))
        errors['SDev_SD'][i_t] = fr(np.std(boot_hhs[i_t].loc['SDev']))
        errors['R2'][i_t] = fr(np.mean(boot_hhs[i_t].loc['R2']))
        errors['R2_SD'][i_t] = fr(np.std(boot_hhs[i_t].loc['R2']))
        
  
    for i_t in col_tech_days:
        #FOR TESTING: i_t = col_tech[0]
        errors['Bias_10'][i_t] = fr(np.percentile(boot_days[i_t].loc['Bias'],10))
        errors['Bias_90'][i_t] = fr(np.percentile(boot_days[i_t].loc['Bias'],90))
        errors['SDev'][i_t] = fr(np.mean(boot_days[i_t].loc['SDev']))
        errors['SDev_SD'][i_t] = fr(np.std(boot_days[i_t].loc['SDev']))
        errors['R2'][i_t] = fr(np.mean(boot_days[i_t].loc['R2']))
        errors['R2_SD'][i_t] = fr(np.std(boot_days[i_t].loc['R2']))
    
    return errors
    
def format_errors(errors):
    # format_errors():      Format errors as strings for make_table()
    # errors:               Dataframe with errors
    #
    #FOR TESTING: errors = err_ft;
    err_text = pd.DataFrame(columns=['Bias','SDev','R2'], index=errors.index)
    err_text.index.name = 'Technique'
    fn = lambda x: '{:04.4f}'.format(x) #format number to string

    for i_t in errors.index:
        # i_t = errors.index.values[0]
        err_text['Bias'][i_t] =  '(' + fn(errors['Bias_10'][i_t]) + ', ' + fn(errors['Bias_90'][i_t]) + ')'
        err_text['SDev'][i_t] =  fn(errors['SDev'][i_t]) + ' ±' + fn(errors['SDev_SD'][i_t])
        err_text['R2'][i_t] =  fn(errors['R2'][i_t]) + ' ±' + fn(errors['R2_SD'][i_t])

    return err_text

def make_table(ini, err_ft, err_dt, err_nt):
    # make_table():      Combine error measures and save to file (bootstrapping results table, e.g. for manuscript)
    
    tab_ft = format_errors(err_ft)
    tab_dt = format_errors(err_dt)
    tab_nt = format_errors(err_nt)
    
    err_table = pd.concat([tab_ft,tab_dt,tab_nt], keys=['ft','dt','nt'], names=['TimeBase',ini.FluxGas+' ('+ini.FluxUnit+')'])
    err_table =  err_table.unstack('TimeBase').stack('TimeBase')

    return err_table

#%% Calculate flux sums and their errors
    
def calc_sums(data, ini, err_ft):
    # calc_sums():      Calculate full period sums with errors from bootstrapping,
    #                   use 'hhs' for short gaps and 'days' for long gaps (mgf.utl.set_gapl_S() --> 12 hhs)
    # data:             Dataframe with timestamp index, fluxes, meteo and co, use 'real' columns
    # ini:              Settings from ini file
    # err_ft:           Errors (uncertainty estimates) from bootstrapping
    
    # Error estimates only work if techniques are available for hhs and days
    tech_hhs = data.filter(regex= 'hhs').columns.str.replace('_hhs','')
    tech_days = data.filter(regex= 'days').columns.str.replace('_days','')
    tech_both = tech_days[tech_days.isin(tech_hhs)]

    # Make sum data frame 
    sums = pd.DataFrame(index=tech_both, columns=['SumObs','SumFillReal','SumTotal','MissFillReal','RandomAll','BiasGaps','ErrorTotal','LowerCI','UpperCI'])
    sums.index.name = ini.FluxGas+' ('+ini.ConvSums+')'
    
    # Properties of original measured fluxes
    flux_org = mgf.utl.name_flux(ini)
    factorHH = float(ini.ConvFactor)
    
    fr = lambda x: mgf.utl.round_dec(x,4) #round number to four digits
        
    # Calculate sums
    for i_t in tech_both:
        #i_t = tech_both[2]
        # Sum over all measured (original) data points
        sums['SumObs'].at[i_t] = fr(data[flux_org].sum() * factorHH)
        # Sum over all filled real gaps        
        sums['SumFillReal'].at[i_t] = fr(data[i_t+'_real'].sum() * factorHH - sums['SumObs'].at[i_t])
        # Add up total sum
        sums['SumTotal'].at[i_t] = fr(sums['SumObs'][i_t] + sums['SumFillReal'][i_t])
        # Number of real gaps that could not be filled
        sums['MissFillReal'].at[i_t] = (data[flux_org].isnull() & data[i_t+'_hhs'].isnull()).sum() #Count only missings during the real gaps
        
        # Sum up uncertainties
        # Calculate random error from SDev for all data points from _hhs since this is the "optimal" model
        sums['RandomAll'].at[i_t] = fr(np.sqrt(data.shape[0] * (err_ft['SDev'].loc[i_t+'_hhs'])**2) * factorHH)
        # Calculate bias error for filled data points depending on gap size
        bias_hhs = max(abs(err_ft['Bias_10'][i_t+'_hhs']),abs(err_ft['Bias_90'][i_t+'_hhs']))
        bias_days = max(abs(err_ft['Bias_10'][i_t+'_days']),abs(err_ft['Bias_90'][i_t+'_days']))
        # Count number of short and long real gaps
        gaps = mgf.utl.countGapsTot(data, flux_org)
        num_S = sum(np.where((gaps != 0) & (gaps <= mgf.utl.set_gapl_S()), True, False))
        num_L = sum(np.where((gaps != 0) & (gaps > mgf.utl.set_gapl_S()), True, False))
        sums['BiasGaps'].at[i_t] = fr((num_S * bias_hhs + num_L * bias_days) * factorHH)
        # Add up the two error types
        sums['ErrorTotal'].at[i_t] = fr(sums['RandomAll'][i_t] + sums['BiasGaps'][i_t])
        
        # Sum confidence intervals
        sums['LowerCI'].at[i_t] = fr(sums['SumTotal'].at[i_t] - sums['ErrorTotal'].at[i_t])
        sums['UpperCI'].at[i_t] = fr(sums['SumTotal'].at[i_t] + sums['ErrorTotal'].at[i_t])
        
    return sums

#%% Calculate flux sums and their errors
    
def calc_ensemble(sums, ini):
    # calc_ensemble():  Calculate ensemble results from sums
    # sums:             Dataframe sums (aggregated fluxes) for each technique
    # ini:              Settings from ini file 
    
    fr = lambda x: mgf.utl.round_dec(x,1) #round number to one digits
        
    res_ens = pd.Series([], dtype=float)
    res_ens.index.name = 'EnsStats'
    res_ens.name = ini.FluxGas+' ('+ini.ConvSums+')'
            
    res_ens['UpperCI'] = fr(sums['UpperCI'].max(axis=0))
    res_ens['UpperUnc'] = fr(sums['UpperCI'].max(axis=0) - sums['SumTotal'].max(axis=0))
    res_ens['UpperTot'] = fr(sums['SumTotal'].max(axis=0))
    res_ens['Delta'] = fr(sums['SumTotal'].max(axis=0) - sums['SumTotal'].min(axis=0))
    res_ens['LowerTot'] = fr(sums['SumTotal'].min(axis=0))
    res_ens['LowerUnc'] = fr(sums['LowerCI'].min(axis=0) - sums['SumTotal'].min(axis=0))
    res_ens['LowerCI'] = fr(sums['LowerCI'].min(axis=0))
    res_ens['TotalCI'] = fr(sums['UpperCI'].max(axis=0) - sums['LowerCI'].min(axis=0))
       
    return res_ens