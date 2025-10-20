import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../processing')
from helpers import extend_with_decreasing_growth

new_ssp_extension = True
sy = 1950
ey = 2200
nyears = 251
scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp460', 'ssp585']
n_sce=5


####### These are global population and GDP data as multimodel mean from ssp definition:
df = pd.read_csv('SSP_IAM_V2_201811.csv')
# Indices of gdp in the csv file
# Population is always 5 rows below
# There are data for several IAMs : [AIM/CGE, GCAM4, IMAGE, MESSAGE-GLOBIOM (no SSP5), REMIND-MAGPIE, WITCH-GLOBIOM]
# GDP data is in billion US$2005 PPP

gdp_global = np.zeros((nyears,n_sce))
pop_global = np.zeros((nyears,n_sce))

time_raw = np.linspace(2010,2100,10)
time = np.linspace(sy, ey, nyears)

for si, sce in enumerate(scenarios):
    SSP=int(sce[3])
    gdp_SSP_raw = np.zeros((11))
    pop_SSP_raw = np.zeros((11))
    nsum = 0
    if SSP == 1:
        lidxs = [2944, 34597, 50675, 61647, 73235]
        padds = [5, 5, 4, 5, 1]
    elif SSP == 2:
        lidxs = [6537, 38037, 55185, 65828, 76358]
        padds = [5, 5, 4, 5, 1]
    elif SSP == 3:
        lidxs = [8996, 40789, 58153, 78628]
        padds = [18, 5, 4,1 ]  
    elif SSP == 4:
        lidxs = [10958, 44229, 82036]
        padds = [5, 5, 1]
    elif SSP == 5:
        lidxs = [14446, 46981, 70079, 84306]
        padds = [5, 5, 18, 1]
    
    for i, lidx in enumerate(lidxs):
        # Skip GCAM4 data, because they miss data for 2005
        if i==1: continue

        gdp_SSP_raw = gdp_SSP_raw + np.asarray(df.values[lidx,5:16], dtype=float)
        pop_SSP_raw = pop_SSP_raw + np.asarray(df.values[lidx+padds[i],5:16], dtype=float)
        nsum+=1
    
    gdp_global[60:151,si] = np.interp(time[60:151], time_raw, gdp_SSP_raw[1:] / nsum) * 1.31 # conversion from USD2005 to USD2019 as is SLIIDERS
    pop_global[60:151,si] = np.interp(time[60:151], time_raw, pop_SSP_raw[1:] / nsum)
    
# SSP4-6.0 and SSP1-2.6 do not match the other scenarios and population data, rebase upon the ssp2-4.5 scenario!
for i in [0,3]: pop_global[60:151,i] = pop_global[60:151,i] * (pop_global[60,1] / pop_global[60,i])

if new_ssp_extension:
    constant_pop_in_years = 50
    constant_gdp_in_years = 50
    
    pop_global = extend_with_decreasing_growth(pop_global, constant_pop_in_years, ndim=2, istart=150)
    gdp_global = extend_with_decreasing_growth(gdp_global, constant_gdp_in_years, ndim=2, istart=150)
else:
    # extend population and gdp data with constant values between 2100 and 2200
    pop_global = np.append(pop_global[:151], np.zeros((100,n_sce))+pop_global[np.newaxis,150,:], axis=0)
    gdp_global = np.append(gdp_global[:151], np.zeros((100,n_sce))+gdp_global[np.newaxis,150,:], axis=0)
    

    
##### UN WORLD population data for SLR input (Land water storage component)
# For overwriting pop data from 2000 to 2015 (start of scenarios)
df = pd.read_csv('WPP2022_TotalPopulation.csv')

pop_UN_1950_2020 = np.zeros((71,n_sce))
pop_UN_1950_2020[:,:] = df.values[:,1][:,np.newaxis]

pop_global[:,:] = np.append(pop_UN_1950_2020[:65] * 1e-3, 1.025*pop_global[65:], axis=0) # SSP data scaled to match UN World data in 2015


data = {
    'Time': time
}

for isc, scenario in enumerate(scenarios):
    data['Global pop '+scenario] = pop_global[:,isc]


df = pd.DataFrame(data)
if new_ssp_extension: df.to_csv('SSP_global_pop_decreasing_growth_continuation.csv', index=False)
else: df.to_csv('SSP_global_pop_constant_value_continuation.csv', index=False)

df.head()
