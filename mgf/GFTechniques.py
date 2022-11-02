#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Antje M. Lucas-Moffat
Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
Please cite the paper if using this code.

For a package description:
    See mgf/__init__.py
    
Functions:
    # gen_index():   Generate index to select data subset from timestamp index
    # def_LUTs():    Define dataframe with look-up tables (LUTs)
    # def_LUT_MDS(): Same as def_LUT() but for MDS only
    # fill_LUT():    Function to fill gaps with LUT
    # calc_IPs():    Calculate interpolation of gaps, two techniques implemented

Note:
    Diurnal interpolation techniques (WDM, FDA, MDA, MDC) are programmed the same way as look-up tables (LUTs).

"""
#%% Initialization
### Imports from python
import numpy as np
import pandas as pd

### Own imports
import mgf

### Options
pd.options.display.max_rows = 20
 
#%% Function for index generator

def gen_index(t_index, win_days, win_hhs, scenario='none'):
    # gen_index():      Generate index to select data subset from timestamp index
    # t_index:          Time index with time stamp at end of half-hour
    # win_hhs:          Window of ± half-hour blocks at same time of day 
    #                   Option win_hhs=-1 for full days
    #                   Option win_hhs=-6 for 3h-periods
    # win_days:         Window of ± adjacent days
    # scenario:         'none' - do nothing (i.e. all half-hours in wanted range)
    #                   'hhs' - wanted index but drop current half-hour (i.e. set to gap)
    #                   'days' - wanted indext but drop half-hours from current day (i.e. set full day to gap) 
    # Generator takes into account that flux data is right closed at the end of half-hour
    #
    # FOR TESTING: t_index=data.index[i_g]; win_hhs=5; win_days=0; scenario='hhs'
    
    index = []
    # Correct time stamp to middle of half-hour to avoid switch of periods/days
    t_index_cor = t_index - pd.DateOffset(hours=0.25)
    
    if win_hhs == -1: #Full days
        #Add preceding days
        start = t_index_cor - pd.DateOffset(days=win_days)
        #Add succeeding days
        end = t_index_cor + pd.DateOffset(days=win_days)
        index = pd.date_range(start=start.strftime('%Y-%m-%d 00:15:00'), end=end.strftime('%Y-%m-%d 23:45:00'), freq='30min')
        
    elif win_hhs == -6: #3h-periods
        jjj = pd.DataFrame(index=[t_index_cor],columns=[0]).resample('180min').asfreq()
        index_3h_start = jjj.index[0] + pd.DateOffset(hours=0.25) #Needs also correction
        index = pd.period_range(index_3h_start, periods=6, freq='30min').to_timestamp()
        index_3h = index
        for i_d in range(1,win_days+1):
            #Add preceding days
            index = index.append(index_3h - pd.DateOffset(days=i_d))
            index = index.append(index_3h + pd.DateOffset(days=i_d))

    else: # Half-hour blocks
        #Same day
        index.append(t_index_cor)
        for i_hh in range(1,win_hhs+1):
            index.append(t_index_cor-pd.DateOffset(hours=i_hh*0.5))
            index.append(t_index_cor+pd.DateOffset(hours=i_hh*0.5))
        #Multiple days
        for i_d in range(1,win_days+1):
            #Add preceding days
            ts_index=t_index_cor-pd.DateOffset(days=i_d)
            index.append(ts_index)
            for i_hh in range(1,win_hhs+1):
                index.append(ts_index-pd.DateOffset(hours=i_hh*0.5))
                index.append(ts_index+pd.DateOffset(hours=i_hh*0.5))
            #Add succeeding days
            ta_index=t_index_cor+pd.DateOffset(days=i_d)
            index.append(ta_index)
            for i_hh in range(1,win_hhs+1):
                index.append(ta_index-pd.DateOffset(hours=i_hh*0.5))
                index.append(ta_index+pd.DateOffset(hours=i_hh*0.5))
                
    # Attention: Mix of DatetimeIndex and Timestamp
    # Sort index and make sure it is DatetimeIndex (and not Timestamp)
    # For win_hhs > 23, unique() removes duplicates from index          
    index = pd.to_datetime(index).sort_values().unique() 
        
    if scenario == 'hhs':
        #Drop current half-hour
        index = index.drop(t_index_cor) #append() inplace but drop() not inplace?!
    elif scenario == 'days':
        #Drop current day
        t_start=t_index_cor.strftime('%Y-%m-%d 00:15:00')
        t_end=t_index_cor.strftime('%Y-%m-%d 23:45:00')
        index= index.drop(pd.date_range(t_start, t_end, freq='30min'), errors='ignore') #ignore since not all indices (of full day) might have been included

    # Correct time stamp back to end of half hour
    index = index + pd.DateOffset(hours=0.25) 
    return index #DatetimeIndex object


"""
#TESTING
    i_g=64; data.index[i_g]
    gen_index(data.index[i_g], win_days=0, win_hhs=2, scenario='days')
    gen_index(data.index[i_g], 1, 2, 'days')
    gen_index(data.index[i_g], 0, 2, 'hhs')
    gen_index(data.index[i_g], 0, 2, 'none')
"""

#%% Define LUT techniques including MDA, MDC, FDA, WDM

def def_LUTs(var_1='none', range_1=np.nan, var_2='none', range_2=np.nan, var_3='none', range_3=np.nan, var_light='none', var_thres=np.nan):
    # def_LUTs():           Define dataframe with look-up tables (LUTs)
    # var_1, range_1, ...:  Variable names and float ranges for LUT
    # var_light, var_thres: Variable name and threshold for daylight / nighttime
    #
    # LUTs data frame:
    # Technique:            Name of technique
    # WinDayStart:          Window size of ± days to start with
    # WinDayStep:           Number of days to increase window size
    # WinHHs:               Number of half-hours to interpolate
    # CondIFs:              if-conditions to look-up-table or to interpolate for others
    #
    # TESTING: var_1='none'; range_1=np.nan; var_2='none'; range_2=np.nan; var_3='none'; range_3=np.nan; var_light='none'; var_thres=np.nan
    # TESTING: var_light='Rg'; var_thres=5

    LUTs = pd.DataFrame(columns=[], index=['Technique','WinDayStart','WinDayStep','WinHHs','CondIFs'])
    if (var_light != 'none') & np.isfinite(var_thres):
        LUTs['WDM'] = ['WDM', 0, 1, -1, '((data.'+var_light+'[i_g] > '+str(var_thres)+') & (df_window.'+var_light+' > '+str(var_thres)+')) | '+
            '((data.'+var_light+'[i_g] <= '+str(var_thres)+') & (df_window.'+var_light+' <= '+str(var_thres)+'))']
    # For FDA, the data is weighted by data is weighted day/night for the whole win_days window and only per day for win_days=0
    LUTs['FDA_hh6'] = ['FDA_hh6', 0, 1, -6, 'df_window.'+var_light+' != True'] #Dummy, but could be real condition
    LUTs['MDA_hh5'] = ['MDA_hh5', 0, 1, 2, 'df_window.'+var_light+' != True'] #Dummy, but could be real condition
    LUTs['MDC_d3'] = ['MDC_d3', 3, 3, 0, 'df_window.'+var_light+' != True'] #Dummy, but could be real condition
    LUTs['MDC_d7'] = ['MDC_d7', 7, 7, 0, 'df_window.'+var_light+' != True'] #Dummy, but could be real condition
    if (var_1 != 'none') & np.isfinite(range_1):
        LUTs['LUT_V1_d3'] = ['LUT_V1_d3', 3, 3, -1, '(abs(df_window.'+var_1+' - data.'+var_1+'[i_g]) <= '+str(range_1)+')']
        LUTs['LUT_V1_d7'] = ['LUT_V1_d7', 7, 7, -1, '(abs(df_window.'+var_1+' - data.'+var_1+'[i_g]) <= '+str(range_1)+') ']
        if (var_2 != 'none') & np.isfinite(range_2) :
            LUTs['LUT_V1V2_d3'] = ['LUT_V1V2_d3', 3, 3, -1, '(abs(df_window.'+var_1+' - data.'+var_1+'[i_g]) <= '+str(range_1)+') & '+
               '(abs(df_window.'+var_2+' - data.'+var_2+'[i_g]) <= '+str(range_2)+')']
            LUTs['LUT_V1V2_d7'] = ['LUT_V1V2_d7', 7, 7, -1, '(abs(df_window.'+var_1+' - data.'+var_1+'[i_g]) <= '+str(range_1)+') & '+
                   '(abs(df_window.'+var_2+' - data.'+var_2+'[i_g]) <= '+str(range_2)+')']
            if (var_3 != 'none') & np.isfinite(range_3):
                LUTs['LUT_V1V2V3_d3'] = ['LUT_V1V2V3_d3', 3, 3, -1, '(abs(df_window.'+var_1+'- data.'+var_1+'[i_g]) <= '+str(range_1)+') & '+
                       '(abs(df_window.'+var_2+' - data.'+var_2+'[i_g]) <= '+str(range_2)+') & '+
                       '(abs(df_window.'+var_3+' - data.'+var_3+'[i_g]) <= '+str(range_3)+')']
                LUTs['LUT_V1V2V3_d7'] = ['LUT_V1V2V3_d7', 7, 7, -1, '(abs(df_window.'+var_1+' - data.'+var_1+'[i_g]) <= '+str(range_1)+') & '+
                       '(abs(df_window.'+var_2+' - data.'+var_2+'[i_g]) <= '+str(range_2)+') & '+
                       '(abs(df_window.'+var_3+' - data.'+var_3+'[i_g]) <= '+str(range_3)+')']
    return LUTs

def def_LUT_MDS(var_1='none', range_1=np.nan, var_2='none', range_2=np.nan, var_3='none', range_3=np.nan):
    # def_LUT_MDS():    Same as above but for defining MDS LUT setup only
    
    LUTs = pd.DataFrame(columns=[], index=['Technique','WinDayStart','WinDayStep','WinHHs','CondIFs'])
    if (var_1 != 'none') & np.isfinite(range_1) & (var_2 != 'none') & np.isfinite(range_2) & (var_3 != 'none') & np.isfinite(range_3):
        LUTs['LUT_MDS_d7'] = ['LUT_MDS_d7', 7, 7, -1, '(abs(df_window.'+var_1+' - data.'+var_1+'[i_g]) <= '+str(range_1)+') & '+
                           '(abs(df_window.'+var_2+' - data.'+var_2+'[i_g]) <= '+str(range_2)+') & '+
                           '(abs(df_window.'+var_3+' - data.'+var_3+'[i_g]) <= '+str(range_3)+')']
    elif (var_1 != 'none') & np.isfinite(range_1) & (var_2 != 'none') & np.isfinite(range_2) & (var_3 == 'none'): #If 3rd variable (VPD) is missing as in NH3 dataset
        LUTs['LUT_MDS_d7'] = ['LUT_MDS_d7', 7, 7, -1, '(abs(df_window.'+var_1+' - data.'+var_1+'[i_g]) <= '+str(range_1)+') & '+
                           '(abs(df_window.'+var_2+' - data.'+var_2+'[i_g]) <= '+str(range_2)+')']
    return LUTs

def fill_LUT(data, flux_org, methLUT, scenario):
    # fill_LUT():       Function to fill gaps with LUT techniques
    # data:             Dataframe with timestamp index, flux measurements, and maybe meteo
    # flux_org:         Data column with original measured fluxes to be filled/used and rest set to nan
    # methLUT:          Gap filling technique as defined in the dataframe 'LUTs' above
    # scenario:         Scenario for gap filling ('hhs', 'days', or 'none', see also gen_index())
    #
    # TESTING: data=data_co2; ini=co2.Settings; flux_org = gen_fcol(data,ini); methLUT=LUTs['WDM']; scenario='hhs'
    
    flux_fm = methLUT.Technique +'_'+ scenario #Flux filled with fill-technique
    print('>>> Fill technique:', flux_fm, '   (', methLUT.CondIFs,')')
    
    # Write gap filled values into new column
    data[flux_fm] = np.nan 
    win_days = methLUT.WinDayStart
    num_gaps = data[flux_fm].isnull().sum()
    num_gaps_org = data[flux_fm].isnull().sum()
    while num_gaps > 0:
        check_gaps = data[flux_fm].isnull()
        array_gaps = np.array(range(0,data.shape[0]))[np.array(check_gaps)]
        # Calculate LUT value
        #print('>', win_days, check_gaps.sum())
        for i_g in array_gaps:
        #TESTING: for i_g in array_gaps[0:5]:
            df_window = data[data.index.isin(gen_index(data.index[i_g], win_days=win_days, win_hhs=methLUT.WinHHs, scenario=scenario))]
            #--> Much faster to first reduce dataset to window BEFORE compound-if-query
            lut_entries = eval(methLUT.CondIFs)
            if (lut_entries.sum() >= 2): #if more than two flux values available (after MR only, difference to EF algorithm)
                data[flux_fm].iat[i_g] = np.mean(df_window[flux_org][lut_entries])

        # Increase window of days and break if too large        
        #Print remaining gaps for comparison with full MDS algorithm
        num_gaps =  data[flux_fm].isnull().sum()
        if (methLUT.Technique == 'LUT_MDS_d7'):
           print('>LUT_MDS: window size:', win_days, ', remaining gaps:', num_gaps, ' of', num_gaps_org,', in percent: ', '{:04.2f}'.format(num_gaps/num_gaps_org*100))

        win_days = win_days + methLUT.WinDayStep
        tot_days = int(data.shape[0] / 48.0)
        if (win_days > 0.5 * tot_days) & (num_gaps > 0):
           print('>! ', num_gaps, 'remaining gaps. Break since window size: ', win_days, ' is larger than half of total days: ', 0.5*tot_days, '.)')
           break

 #%% Define interpolation techniques
    
def calc_IPs(data, flux_org, technique, scenario):
    # calc_IPs():       Calculate interpolation of gaps, two techniques implemented
    # data:             Dataframe with timestamp index, flux measurements, and maybe meteo
    # flux_org:         Data column with original measured fluxed to be filled/used and rest set to nan
    # technique:        Interpolation technique 'IP_lin' or 'IP_mov'
    # scenario:         Scenario for gap filling ('hhs', 'days', or 'none', see also gen_index())
    #
    # TESTING: flux_org=mgf.utl.name_flux(ini); technique='IP_mov'; scenario='hhs'; scenario='days';
    
    flux_fm = technique +'_'+ scenario #Flux filled with fill-technique
    print('>>> Fill technique:', flux_fm)
    data[flux_fm]=np.nan
    # For linear interpolation, all real gaps could be filled at once.
    # Loop only need to simulate artificial gaps, each at a time.
    # Limit interpolation to certain day window to speed up processing.
    win_days = 1 #must be at least ±1 for linear interpolation and also for days scenario
    gapl_S = mgf.utl.set_gapl_S() #pre-definded length of short gaps
    num_gaps = data[flux_fm].isnull().sum()
    while num_gaps > 0:
        check_gaps = data[flux_fm].isnull()
        array_gaps = np.array(range(0,data.shape[0]))[np.array(check_gaps)]
        # Calculate IP value
        #print('>', win_days, check_gaps.sum())
        for i_g in array_gaps:
        #TESTING: for i_g in array_gaps[0:20]: # for i_g in array_gaps[16696:(16696+20)]: # i_g = 0
            # limit = maximum number of consecutive values set to (win_days*48) to ensures that gap does not range to end of ±win_days period
            df_window = data[data.index.isin(gen_index(data.index[i_g], win_days=win_days, win_hhs=-1, scenario='none'))] #Continuous index without scenario
            entries = df_window.index.isin(gen_index(data.index[i_g], win_days=win_days, win_hhs=-1, scenario=scenario)) #Index with scenario
            df_window.insert(0,'data',df_window[flux_org].where(entries))
            df_window.insert(1,'gaps',mgf.utl.countGapsTot(df_window,'data'))
           # Linear interpolation in two steps
            # Attention: 'limit' does not limit gap size but only how many values are interpolated and 'inside','outside' did not work as expected...
            if technique == 'IP_lin':
                if (df_window['gaps'].loc[data.index[i_g]] <= gapl_S) : #Gaps are filled if short, at least start and end are available
                # Interpolate with pandas.interpolate()
                    data[flux_fm].iat[i_g] = df_window['data'].interpolate(technique='time', axis=0, limit=gapl_S, limit_direction='both').loc[data.index[i_g]] #So cool, single line routine!  
#                elif (df_window['gaps'].loc[data.index[i_g]] < gapl_S) : #If within gapl_S interpolate with quarterly means
#                    df_window.is_copy = False
#                    # Calculate fixed 3h-mean and interpolate missing 3h-means
#                    df_window['fix_3h'] = df_window['data'].groupby(pd.Grouper(freq='180min',closed='right',label='left')).apply(lambda g: g.mean()).interpolate(technique='time', axis=0, limit=win_days, limit_direction='both') 
#                    # Propagate NEXT valid observation backwards to fill gap (i.e. pure propagation, no further interpolation)
#                    data[flux_fm].iat[i_g] = df_window['fix_3h'].bfill().loc[data.index[i_g]]
                else: # Otherwise use daily means
                    df_win_days = df_window['data'].groupby(pd.Grouper(freq='1d',closed='right',label='left')).apply(lambda g: g.mean()).interpolate(technique='time', axis=0, limit=win_days, limit_direction='both')
                    data[flux_fm].iat[i_g] = df_win_days.loc[(data.index[i_g]- pd.DateOffset(hours=0.25)).strftime('%Y-%m-%d')] ## Use shifted time stamp to middle of half-hour to avoid midnight switch
            elif technique == 'IP_mov': #Moving average interpolation
                if (df_window['gaps'].loc[data.index[i_g]] <= gapl_S):   # Only up to gapl_S hh as for IP_lin
                # Rolling = moving window with increasing window size.
                # Moving average (= mean) to fill wholes in time series (= interpolation).
                    roll_start = 5
                    roll_step = 2
                    min_hh = 2 #Minimum of available half-hours
                    window = roll_start
                    data[flux_fm].iat[i_g] = df_window['data'].rolling(window, center=True, min_periods=min_hh).mean().loc[data.index[i_g]]
                    #while np.isnan(data[flux_fm].iat[i_g]) & ((window*0.5) < (win_days*48)): # Since center=True option, half the window critical
                    while np.isnan(data[flux_fm].iat[i_g]):
                        window = window + roll_step
                        data[flux_fm].iat[i_g] = df_window['data'].rolling(window, center=True, min_periods=min_hh).mean().loc[data.index[i_g]]
                else:
                    min_days=2 #At least 2 for center of three days
                    df_win_days = df_window['data'].groupby(pd.Grouper(freq='1d',closed='right',label='left')).apply(lambda g: g.mean()).rolling(window=min_days+win_days, center=True, min_periods=min_days).mean()
                    data[flux_fm].iat[i_g] = df_win_days.loc[(data.index[i_g]- pd.DateOffset(hours=0.25)).strftime('%Y-%m-%d')] ## Use shifted time stamp to middle of half-hour to avoid midnight switch
        # Increase window of days and break if too large
        win_days = win_days + 1
        num_gaps =  data[flux_fm].isnull().sum()
        if (win_days > 10) & (num_gaps > 0):
           print('>!', num_gaps, 'remaining gaps. Break since window size: ', win_days, ' is larger than ±10 for interpolations.)')
           break

