#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Multiple Gap-Filling Tool
************************
@author: Antje M. Lucas-Moffat
Started: Aug 16 2018 11:35:08 
Updated: Aug 12 2022

Original code of:
    Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
!!! Please cite the paper if using this code !!!
Further details can be found therein.


Short description of scenarios and techniques
************************
Gap filling scenarios
    'hhs':  For each half-hour at a time, set to artificial gap.
    'days': For each day at a time, half-hours of whole day set to artificial gaps.
    'real': Filling of real gaps automatically included in 'hhs' scenario.

Gap filling techniques (could also be called interpolation techniques or methods)
    *** Simple interpolation methods
    'IP_lin': Linear interpolation (for 'days' scenario, linear interpolation of daily means)
    'IP_mov': Interpolation with moving (rolling) average
    *** Diurnal interpolation methods (after Falge, Reichstein), programmed same way as look-up-tables:
    'WDM':    Weighted Daylight Mean
    'FDA':    Fixed Diurnal Average
    'MDA':    Moving Diurnal Average (Reichstein)
    'MDC':    Mean Diurnal Course (Falge)
    *** Various look-up tables, LUT (after Falge, Reichstein):
    'LUT_V1':       Look-up tables with 1 dependent variable
    'LUT_V1V2':     Look-up tables with 2 dependent variables
    'LUT_V1V2V3':   Look-up tables with number of dependent variables from 1 to 3
    *** Results from other techniques/methods/models:
    'ANN_XYZ:       Imported gap filling results from my ANNs
    'Model_XYZ':    Additional imported gap filling results


Overview of scripts
***********************
'MGFilling.py':         Multiple gap filling (fill fluxes, analyze, plot) to be used with Jupyter Notebooks 
'Utils.py':             Utilities, i.e. helpful small functions
'GFTechniques.py':      Implementation of gap filling techniques
'FFAnalyze.py':         Analyze filled fluxes (e.g. by bootstrapping)
7x 'Plot....py':        Scripts for plots of data, results, tables

"""

# Description
__name__ = "mgf"
__version__ = "beta_20220812"
__author__ = "Antje Lucas-Moffat"
__email__ = "amm@moffats.de"
# For testing and mansucript, e.g. include extra plot features
__TEST__ = False
__MS__ = True

#### Import own scripts

# Script for filling and analyzing gaps
import mgf.MGFilling as mgf
import mgf.Utils as utl
import mgf.GFTechniques as gft
import mgf.FFAnalyze as ffa

# Scripts for plotting
import mgf.PlotBoxStats as pbs
import mgf.PlotSums as pse
import mgf.PlotScatter as psc
import mgf.PlotGapDistr as pgd
import mgf.PlotTimeSeries as pts
import mgf.PlotDailySeries as pds
import mgf.PlotTable as pta
