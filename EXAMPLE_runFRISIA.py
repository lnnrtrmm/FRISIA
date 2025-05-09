
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from src.SLRModel import globalSLRModel
from src.SLRImpactsModel import SLRImpactModel

########## FAIR input data for SLR model ########################
scenarios = ['ssp119', 'ssp245', 'ssp585']
colors_sce = ['lightblue', 'tab:green', 'orange']

Temperature = []
OHC = []

# Load only data between 1850 and 2200!
si = 100
ei = 451

for scenario in scenarios:
    df = pd.read_csv('input/FaIR/T_and_OHC_from_FaIRv2.1.1_calibrationv1.1.0_'+scenario+'.csv')

    Temperature.append(df.values[si:ei,1:4] - np.mean(df.values[100:150,np.newaxis,1:4], axis=0)) # take anomaly to 1850-1900
    OHC.append(df.values[si:ei,4:7] - np.mean(df.values[100:150,np.newaxis,4:7], axis=0))
    
time = df.values[:,0]
Temperature = np.asarray(Temperature)
OHC = np.asarray(OHC)
OHC_change = np.append(np.zeros_like(OHC[:,0,np.newaxis,:]), OHC[:,1:,:] - OHC[:,:-1,:], axis=1)



########## Population and GDP data ##############################
### This is for global GDP and pop, as well as coastal information
### The file was prepared using input data from the IIASA SSP database (multimodel mean data)
### and the DIVA spreadsheet taken from the repository of the Wong et al. (2022) CIAM paper.
###
### There are 22 columns per SSP scenario for the three scenarios SSP1-1.9, SSP2-4.5 and SSP5-8.5
ncol = 22
df = pd.read_csv('input/coastal_GDP_and_population.csv')

global_GDP_2010_2200 = np.zeros((191,3))
global_pop_2010_2200 = np.zeros((191,3))


coastal_gdp_global = np.zeros((191,3))
coastal_gdp_bipolar = np.zeros((191,2,3))
coastal_gdp_regional = np.zeros((191,7,3))
coastal_pop_global = np.zeros((191,3))
coastal_pop_bipolar = np.zeros((191,2,3))
coastal_pop_regional = np.zeros((191,7,3))
for i in range(3):
    global_GDP_2010_2200[:,i] = df.values[:,1+i*ncol]
    global_pop_2010_2200[:,i] = df.values[:,12+i*ncol]
    coastal_gdp_global[:,i] = df.values[:,2+i*ncol]
    coastal_pop_global[:,i] = df.values[:,13+i*ncol]
    coastal_gdp_bipolar[:,:,i] = df.values[:,3+i*ncol:5+i*ncol]
    coastal_pop_bipolar[:,:,i] = df.values[:,14+i*ncol:16+i*ncol]
    coastal_gdp_regional[:,:,i] = df.values[:,5+i*ncol:12+i*ncol]
    coastal_pop_regional[:,:,i] = df.values[:,16+i*ncol:23+i*ncol]
    
coastal_gdp_bipolar = coastal_gdp_bipolar.transpose(1, 0, 2)
coastal_pop_bipolar = coastal_pop_bipolar.transpose(1, 0, 2)
coastal_gdp_regional = coastal_gdp_regional.transpose(1, 0, 2)
coastal_pop_regional = coastal_pop_regional.transpose(1, 0, 2)

    
##### UN WORLD population data for SLR input (Land water storage component)
df = pd.read_csv('input/WPP2022_TotalPopulation.csv')

pop_UN_1950_2020 = np.zeros((71,3))
pop_UN_1950_2020[:,:] = df.values[:,1][:,np.newaxis]

global_pop_1950_2200 = np.append(pop_UN_1950_2020[:60]*1.0e-3, global_pop_2010_2200, axis=0)



##############################################################
##
## CO2 emissions from RCMIP 
## -> these are total CO2 emissions as used in 
##    MPI-ESM simulation (1850-2300)
##
##############################################################
all_C_emissions_list = []
for sce in scenarios:
    infile = 'input/CO2_emissions/global_C_emission_'+sce+'.nc'
    ds=Dataset(infile)
    C_emission = ds.variables['carbon_emission_global'][1:,0,0]
    ds.close()
    
    all_C_emissions_list.append(C_emission)
    
all_C_emissions = np.zeros_like(Temperature[:,:,0])
for i in range(3): all_C_emissions[i,:] = np.copy(all_C_emissions_list[i][:-100]) # remove last 100, to only keep data until 2200


############### Calculating SLR Model
SLR_total = np.zeros_like(Temperature)
SLR_thermo = np.zeros_like(Temperature)
SLR_MG = np.zeros_like(Temperature)
SLR_LWS = np.zeros_like(Temperature)
SLR_GIS = np.zeros_like(Temperature)
SLR_AIS = np.zeros_like(Temperature)

population = np.append(np.zeros((100,3)), global_pop_1950_2200, axis=0) # population data only needed after 1950, but array has to cover the full length
SLRModel = globalSLRModel(Temperature[0,:,0], OHC_change[0,:,0], population=population[:,1], sy=1850, ey=2200)


for sidx, scenario in enumerate(scenarios):    
    for pidx in range(3):

        SLRModel.T_anomaly = Temperature[sidx,:,pidx]
        SLRModel.OHC_change = OHC_change[sidx,:,pidx]
        SLRModel.population = population[:,sidx]
        SLRModel.reset_SLR()
        SLRModel.integrate()
        SLRModel.align(2010, silent=True)
        

        SLR_total[sidx,:,pidx] = SLRModel.getSLRTotal()
        SLR_thermo[sidx,:,pidx] = SLRModel.getSLRThermo()
        SLR_MG[sidx,:,pidx] = SLRModel.getSLRMG()
        SLR_LWS[sidx,:,pidx] = SLRModel.getSLRLWS()
        SLR_GIS[sidx,:,pidx] = SLRModel.getSLRGIS()
        SLR_AIS[sidx,:,pidx] = SLRModel.getSLRAIS()
        
# Reduce data arrays even further to period 2010 - 2200!
SLR        = SLR_total[:,160:,:]
SLR_thermo = SLR_thermo[:,160:,:]
SLR_LWS    = SLR_LWS[:,160:,:]
SLR_MG     = SLR_MG[:,160:,:]
SLR_GIS    = SLR_GIS[:,160:,:]
SLR_AIS    = SLR_AIS[:,160:,:]



############### Running FRISIA
T = Temperature[:,160:451,:]
CO2emis = all_C_emissions[:,160:451]

time = np.linspace(2010, 2200, 191)

nyears = 191
nsce = 3
nens = 3
ndim = 1 # global FRISIA setup



##### No protection, no retreat ("No adaptation")

### CIAM comparison
total_relocation_cost_global_no_adapt = np.zeros((nyears, nsce, nens))
total_flood_cost_global_no_adapt = np.zeros((nyears, nsce, nens))
total_storm_cost_global_no_adapt = np.zeros((nyears, nsce, nens))
total_construct_cost_global_no_adapt = np.zeros((nyears, nsce, nens))
total_people_flooded_global_no_adapt = np.zeros((nyears, nsce, nens))
total_flood_fatalities_global_no_adapt = np.zeros((nyears, nsce, nens))

total_relocation_cost_global_protect = np.zeros((nyears, nsce, nens))
total_flood_cost_global_protect = np.zeros((nyears, nsce, nens))
total_storm_cost_global_protect = np.zeros((nyears, nsce, nens))
total_construct_cost_global_protect = np.zeros((nyears, nsce, nens))
total_people_flooded_global_protect = np.zeros((nyears, nsce, nens))
total_flood_fatalities_global_protect = np.zeros((nyears, nsce, nens))


for sidx in range(3):
    print('Start calculating scenario '+scenarios[sidx])
    
    for i in range(nens):
        ### No adaptation strategy
        SLR_components = [SLR_thermo[sidx,:,i], SLR_LWS[sidx,:,i], SLR_MG[sidx,:,i], SLR_GIS[sidx,:,i], SLR_AIS[sidx,:,i]]
        Model = SLRImpactModel(SLR[sidx,:,i], T[sidx,:,i], CO2emis[sidx,:], coastal_pop_global[np.newaxis,:,sidx],
                               coastal_gdp_global[np.newaxis,:,sidx], ey=2200, ndim=ndim, version='DIVA_global',
                               input_path='input/', include_initial_fp=True, randomize=False,
                               include_SLR_components=True, SLR_components=SLR_components)
        
        Model.willingness_to_retreat = 0.0
        Model.willingness_to_invest_in_fp = 0.0
                
        Model.integrate()  
        
        total_relocation_cost_global_no_adapt[:,sidx,i] = Model.getTotalRelocationCost()
        total_flood_cost_global_no_adapt[:,sidx,i] = Model.getTotalFloodCost()
        total_storm_cost_global_no_adapt[:,sidx,i] = Model.getTotalStormCost()
        total_people_flooded_global_no_adapt[:,sidx,i] = Model.getTotalPeopleFlooded()
        total_flood_fatalities_global_no_adapt[:,sidx,i] = Model.getTotalFloodFatalities()
        
        
        ### Protect strategy
        SLR_components = [SLR_thermo[sidx,:,i], SLR_LWS[sidx,:,i], SLR_MG[sidx,:,i], SLR_GIS[sidx,:,i], SLR_AIS[sidx,:,i]]
        Model = SLRImpactModel(SLR[sidx,:,i], T[sidx,:,i], CO2emis[sidx,:], coastal_pop_global[np.newaxis,:,sidx],
                               coastal_gdp_global[np.newaxis,:,sidx], ey=2200, ndim=ndim, version='DIVA_global',
                               input_path='input/', include_initial_fp=True, randomize=False,
                               include_SLR_components=True, SLR_components=SLR_components)
        
        Model.willingness_to_retreat = 0.0
        Model.willingness_to_invest_in_fp = 1.0
                
        Model.integrate()  
        
        total_relocation_cost_global_protect[:,sidx,i] = Model.getTotalRelocationCost()
        total_flood_cost_global_protect[:,sidx,i] = Model.getTotalFloodCost()
        total_storm_cost_global_protect[:,sidx,i] = Model.getTotalStormCost()
        total_people_flooded_global_protect[:,sidx,i] = Model.getTotalPeopleFlooded()
        total_flood_fatalities_global_protect[:,sidx,i] = Model.getTotalFloodFatalities()
