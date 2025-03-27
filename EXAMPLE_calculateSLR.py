import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from SLRModel import globalSLRModel

########## FAIR input data for SLR model ########################
scenarios = ['ssp119', 'ssp245', 'ssp585']
colors_sce = ['lightblue', 'tab:green', 'orange']

Temperature = []
OHC = []

for scenario in scenarios:
    df = pd.read_csv('input/FaIR/T_and_OHC_from_FaIRv2.1.1_calibrationv1.1.0_'+scenario+'.csv')

    Temperature.append(df.values[:,1:4] - np.mean(df.values[100:150,np.newaxis,1:4], axis=0)) # take anomaly to 1850-1900
    OHC.append(df.values[:,4:7] - np.mean(df.values[100:150,np.newaxis,4:7], axis=0))
    
time = df.values[:,0]
Temperature = np.asarray(Temperature)
OHC = np.asarray(OHC)
OHC_change = np.append(np.zeros_like(OHC[:,0,np.newaxis,:]), OHC[:,1:,:] - OHC[:,:-1,:], axis=1)




############### Calculating SLR Model
SLR_total = np.zeros_like(Temperature)
SLR_thermo = np.zeros_like(Temperature)
SLR_MG = np.zeros_like(Temperature)
SLR_LWS = np.zeros_like(Temperature)
SLR_GIS = np.zeros_like(Temperature)
SLR_AIS = np.zeros_like(Temperature)

# Here we do not provide population data, so that the BRICK formulation is used for SLR_LWS.
SLRModel = globalSLRModel(Temperature[0,:,0], OHC_change[0,:,0])

for sidx, scenario in enumerate(scenarios):    
    for pidx in range(3):

        SLRModel.T_anomaly = Temperature[sidx,:,pidx]
        SLRModel.OHC_change = OHC_change[sidx,:,pidx]
        SLRModel.reset_SLR()
        SLRModel.integrate()
        SLRModel.align(1993, silent=True) # align in 1993

        SLR_total[sidx,:,pidx] = SLRModel.getSLRTotal()
        SLR_thermo[sidx,:,pidx] = SLRModel.getSLRThermo()
        SLR_MG[sidx,:,pidx] = SLRModel.getSLRMG()
        SLR_LWS[sidx,:,pidx] = SLRModel.getSLRLWS()
        SLR_GIS[sidx,:,pidx] = SLRModel.getSLRGIS()
        SLR_AIS[sidx,:,pidx] = SLRModel.getSLRAIS()



fig = plt.figure(figsize=(14,6), facecolor='white')
font=15

plt.rcParams.update({'font.size': font})
plt.rc('xtick', labelsize=font)
plt.rc('ytick', labelsize=font)

axes = []


gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1.5, 1, 1])

# Erstellung der Subplots
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
ax3 = fig.add_subplot(gs[0, 2], sharey=ax1)

axes = [ax1, ax2, ax3]
for si in range(3):
        
    ax = axes[si]
        
    d1 = np.median(SLR_thermo[si,:,:], axis=1)*1e3
    d2 = np.median(SLR_LWS[si,:,:], axis=1)*1e3
    d3 = np.median(SLR_MG[si,:,:], axis=1)*1e3
    d4 = np.median(SLR_GIS[si,:,:], axis=1)*1e3
    d5 = np.median(SLR_AIS[si,:,:], axis=1)*1e3


    if si==0: ax.set_ylabel('SLR (mm)')
    
    ax.set_xlabel('')
    if si==0:
        ax.set_xlim([1990,2200])
        ax.set_title(scenarios[si], fontsize=font+2)
    else:
        ax.set_xlim([1990, 2200])
        ax.set_title(scenarios[si], fontsize=font+2)
        
    if si == 1: ax.tick_params(labelleft=False)
    if si==2: ax.yaxis.tick_right()
        
    #ax.plot(time, np.median(data[:,si,:], axis=1), label='FRISIAv1.0 (median)', linewidth=5, color=colors_ref[si])
    ax.stackplot(time, d1, d2, d3, d4, d5, alpha=1.0, edgecolor='k',
                 colors=['tab:red', 'tab:orange', 'tab:green', 'lightblue', 'darkblue'],
                 labels=['thermal expansion', 'land water storage change', 'Mountain glaciers', 'Greenland ice sheet', 'Antarctic ice sheet'])
    
    ax.set_xticks([2000,2050,2100,2150])
        
    ax.grid()

                            
ax1.legend(loc='upper left', fontsize=font-1)
ax1.set_ylim(0,3000)

plt.subplots_adjust(left=0.1, right=0.95, wspace=0.03)
plt.savefig('SLR_components.png')
