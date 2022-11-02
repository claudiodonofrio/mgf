#!/usr/bin/env python
# coding: utf-8

# # Multiple Gap Filling Tool
# 
# 
# ### Original code of 
# 
#     Lucas-Moffat et al., 2022, "Multiple gap-filling for eddy covariance datasets", AgrForMet.
# 
# <font size="3"> This notebook will load the example data, run the gap-filling, and analyze the results. <br> Please cite the paper if using this code. </font>

# Check python environment
# Please make sure the module versions are compatible. ***
# Output of original environment (11 August 2022) ***
# platform os: Darwin 21.4.0
# python: 3.8.8 (incl. modules sys, datetime,shutil)
# numpy: 1.20.3
# pandas: 1.3.5
# matplotlib: 3.5.0
# decimal: 1.70
# scipy: 1.7.3

# Check python libraries
print("*** Check current environment ***")
import platform
print("platform os:", platform.system(), platform.release())
print("python:", platform.python_version(), "(incl. modules sys, datetime,shutil)")
import numpy
print("numpy:", numpy.__version__)
import pandas
print("pandas:", pandas.__version__)
import matplotlib
print("matplotlib:",matplotlib.__version__)
import decimal
print("decimal:",decimal.__version__)
import scipy
print("scipy:",scipy.__version__)


# Import necessary libraries

#Import pandas for reading ini-file
import pandas as pd
#Import multiple gap filling package
import mgf
print("Version mgf:", mgf.__version__)


# Load example data and settings


#Set file path and ini file name
file_path='./example/DE-BaF_tNr'
ini_name='tNr_BaF_2016_v03.ini'
# run_number = '202208112034' #Old run with results saved in example directory

#Print ini
pd.read_csv(file_path+'/'+ini_name, sep=',', encoding='latin1', index_col=['Variable'])

# Run multiple gap-filling techniques (GFT)

#Fill artificial gaps of length 'hhs' and 'days' with multiple techniques for the whole dataset
run_number = mgf.mgf.run_GFT(ini_name, file_path)

# Inspect the filled fluxes

### Load GFT results and (if available) model data, fill real gaps, save and plot results
mgf.mgf.inspect_FF(file_path, run_number)

# Bootstrap scenarios from the filled fluxes

#Bootstrap 999 'hhs' scenarios and 999 'days' scenario to evaluate the performances of the GFT¶
mgf.mgf.bootstrap_FF(file_path, run_number)


# Analyse bootstrap results

# Plot results, calculate sums and errors
mgf.mgf.analyse_BS(file_path, run_number)

# Pick emsemble

#Set GFT to be used for the ensemble results
good_GFTs = 'FDA_hh6|MDA_hh5|LUT_V1V2_d7|LUT_V1V2V3_d7|ANN_CRU|ANN_all'

# Plot ensemble results
mgf.mgf.pick_GE(file_path, run_number, good_GFTs)

