import numpy as np
from scipy.stats import gumbel_r
import sys

#### SLR

# This function is taken from the BRICK model:
##==============================================================================
## Copyright 2016 Tony Wong, Alexander Bakker
## This file is part of BRICK (Building blocks for Relevant Ice and Climate
## Knowledge). BRICK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## BRICK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with BRICK.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================
# It was translated into python and simplified by Lennart Ramme, MPI-M, Hamburg 2024
def getLocalSLR(lat, lon, SLRs, weights):
    
    # convert input longitude to degrees east and get indices
    if lon<0: lon += 360
    ilat = int(np.floor(lat)+90)
    ilon = int(np.floor(lon))
    
    wt_thermo = 1.0
    wt_LWS = 1.0
    
    wt_MG = float('NaN')
    shift=0
    while np.isnan(np.sum(wt_MG)):
        shift+=1
        wt_MG = np.nanmean(weights[0][ilat-shift:ilat+shift+1, ilon-shift:ilon+shift+1])
        wt_GIS = np.nanmean(weights[1][ilat-shift:ilat+shift+1, ilon-shift:ilon+shift+1])
        wt_AIS = np.nanmean(weights[2][ilat-shift:ilat+shift+1, ilon-shift:ilon+shift+1])
    
    if np.isnan(np.sum(wt_MG)) or np.isnan(np.sum(wt_GIS)) or np.isnan(np.sum(wt_AIS)):
        sys.exit('Only NaNs found around this location: '+str(lat)+' '+str(lon)+' !')
    
    return wt_thermo*SLRs[0] + wt_LWS*SLRs[1] + wt_MG*SLRs[2] + wt_GIS*SLRs[3] + wt_AIS*SLRs[4]

def getLocalSLR_map(lat, lon, SLRs, weights):
    # implemented to only work without time dimension
    # SLRs should be a list of 5 scalars, one for each SLR compoennt (thermo, LWS, MG, GIS, AIS)

    # lat lon are 2d maps
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()

    # convert input longitude to degrees east
    lon_flat = np.where(lon_flat < 0, lon_flat + 360, lon_flat)

    local_SLR = np.empty_like(lat_flat, dtype=float)
    
    nlat, nlon = weights[0].shape

    for i in range(lat_flat.size):

        ilat = int(np.floor(lat_flat[i]) + 90)
        ilon = int(np.floor(lon_flat[i])) % nlon

        # clamp latitude index
        ilat = np.clip(ilat, 0, nlat - 1)
        
        wt_thermo = 1.0
        wt_LWS = 1.0
    
        wt_MG = float('NaN')
        shift = 0
        while np.isnan(wt_MG):
            shift += 1
            i1 = max(ilat - shift, 0)
            i2 = min(ilat + shift + 1, nlat)
            j1 = max(ilon - shift, 0)
            j2 = min(ilon + shift + 1, nlon)

            wt_MG = np.nanmean(weights[0][i1:i2, j1:j2])
            wt_GIS = np.nanmean(weights[1][i1:i2, j1:j2])
            wt_AIS = np.nanmean(weights[2][i1:i2, j1:j2])

            # if we've exhausted the whole map and still get NaN, raise
            if i1 == 0 and i2 == nlat and j1 == 0 and j2 == nlon and np.isnan(wt_MG):
                raise ValueError(f'Only NaNs found around this location: {lat_flat[i]} {lon_flat[i]} !')

        local_SLR[i] = wt_thermo*SLRs[0] + wt_LWS*SLRs[1] + wt_MG*SLRs[2] + wt_GIS*SLRs[3] + wt_AIS*SLRs[4]
    
    return np.reshape(local_SLR, lat.shape)


def getlSLRweights(lat,lon,weights):
    # convert input longitude to degrees east and get indices
    if lon<0: lon += 360
    ilat = int(np.floor(lat)+90)
    ilon = int(np.floor(lon))
    
    wt_MG = float('NaN')
    shift=0
    while np.isnan(np.sum(wt_MG)):
        shift+=1
        wt_MG = np.nanmean(weights[0][ilat-shift:ilat+shift+1, ilon-shift:ilon+shift+1])
        wt_GIS = np.nanmean(weights[1][ilat-shift:ilat+shift+1, ilon-shift:ilon+shift+1])
        wt_AIS = np.nanmean(weights[2][ilat-shift:ilat+shift+1, ilon-shift:ilon+shift+1])
    
    if np.isnan(np.sum(wt_MG)) or np.isnan(np.sum(wt_GIS)) or np.isnan(np.sum(wt_AIS)):
        sys.exit('Only NaNs found around this location: '+str(lat)+' '+str(lon)+' !')
        
    return wt_MG, wt_GIS, wt_AIS





#### Pre-processing of DIVA/CIAM

#### Functions taken from CIAM
## var is SLR and areaparams are taken from DIVA
## For no adaptation, this should be straightforward in a global model.
## For perfect protection, this should be 0 as long as protection is good enough.
def calcInundatedArea(areaparams, var):
    
    area = (areaparams[0] * np.maximum(0, np.minimum(0.5, var - 0)) \
            + (areaparams[0] + areaparams[1]) / 2 * np.maximum(0, np.minimum(1, var - 0.5)) \
            + areaparams[1] * np.maximum(0, np.minimum(0.5, var - 1.5)) \
            + areaparams[2] * np.maximum(0, np.minimum(1, var - 2)) \
            + areaparams[3] * np.maximum(0, np.minimum(1, var - 3)) \
            + areaparams[4] * np.maximum(0, np.minimum(1, var - 4)) \
            + areaparams[5] * np.maximum(0, np.minimum(1, var - 5)) \
            + areaparams[6] * np.maximum(0, np.minimum(1, var - 6)) \
            + areaparams[7] * np.maximum(0, np.minimum(1, var - 7)) \
            + areaparams[8] * np.maximum(0, np.minimum(1, var - 8)) \
            + areaparams[9] * np.maximum(0, np.minimum(1, var - 9)) \
            + areaparams[10] * np.maximum(0, np.minimum(1, var - 10)) \
            + areaparams[11] * np.maximum(0, np.minimum(1, var - 11)) \
            + areaparams[12] * np.maximum(0, np.minimum(1, var - 12)) \
            + areaparams[13] * np.maximum(0, np.minimum(1, var - 13)) \
            + areaparams[14] * np.maximum(0, var - 14))

    return area


def calc_pSIGMA(pParams, lslr, H):
    # Calculates the expected value of the storm surge exposure area for 
    # Protect scenarios (i.e. with initial flood protection)
    #
    # pParams[0] = psig0 in CIAM
    # pParams[1] = psig0coef in CIAM
    # pParams[2] = psigA in CIAM
    # pParams[3] = psigB in CIAM
    
    return (pParams[0] + pParams[1] * np.maximum(0, lslr)) / (1.0 + pParams[2] * np.exp(pParams[3] * np.maximum(0, H - lslr)))

# Not actually used, but useful for comparing FLOPROS approach with dike heights from DIVA
def create_FLOPROS_DIVA(seg_ypcc, seg_gdp, seg_pop, seg_length, params, popdens, s1, s10, s100, s1000):
    #### Calculating the initial flood protection height using the FLOPROS approach (Scussolini et al., 2016)
    #### and segment information from DIVA database
    
    # WB groups the world by GDP per capita in USD2024 into 4 distinct groups
    # for DIVA data from Wong et al. (2022) is used, these have GDP per capita in USD2010
    inflation_factor = 1.0/1.44      
    
    bound1 = 1135. * inflation_factor
    bound2 = 4495. * inflation_factor
    bound3 = 13935. * inflation_factor
    
    FLOPROS_min = 2.0 # minimum protection level in maximum return period against which there will be protection
    FLOPROS_max = 1000.0 # same but for maximum protection level
    
    mask = np.where(seg_ypcc<=bound1, 1, np.where(seg_ypcc<=bound2, 2, np.where(seg_ypcc<=bound3, 3, 4)))
        
    mean_ypcc = np.zeros(4)
    for i in range(4):
        mean_ypcc[i] = np.sum(np.where(mask==i+1, seg_gdp, 0))/np.sum(np.where(mask==i+1, seg_pop, 0))
    
    grouped_F_mins = (mean_ypcc/mean_ypcc[0]) * FLOPROS_min
    grouped_F_maxs = (mean_ypcc/mean_ypcc[-1]) * FLOPROS_max
    
    # the expected annually exposed assets (using GDP as a proxy for assets) per segment length (EAEA_length)
    # are used here instead of EAD_area in FLOPROS.
    # For this we use the exposure area function from CIAM with 0 as input for lslr and H.
    seg_EAEA_length = np.minimum(seg_pop, calc_pSIGMA(params, 0, 0)*popdens) * seg_ypcc / seg_length
    
    grouped_EAEA_min = np.zeros(4)
    grouped_EAEA_max = np.zeros(4)
    seg_FLOPROS = np.zeros_like(seg_ypcc)
    for i in range(4):
        grouped_EAEA_min = np.nanmin(np.where(mask==i+1,seg_EAEA_length,np.nan))
        grouped_EAEA_max = np.nanmax(np.where(mask==i+1,seg_EAEA_length,np.nan))
        
        seg_FLOPROS = np.where(mask==i+1,
                               interpolate(grouped_F_mins[i],
                                           grouped_F_maxs[i],
                                           seg_EAEA_length,
                                           grouped_EAEA_min,
                                           grouped_EAEA_max),
                               seg_FLOPROS)
        
    # Turning the FLOPROS into a flood protection height (FLOPROH) by linearly interpolating
    # between the enclosing surge levels.
    # For this, we set the initial annual fp height to 0 if there are less than 1 person
    # per km of coastline living in the segment
    seg_FLOPROH = np.where(seg_pop/seg_length < 1, 0, 
                           np.where(seg_FLOPROS <= 10,
                                    interpolate(surge1, surge10, seg_FLOPROS, 1, 10),
                                    np.where(seg_FLOPROS <= 100,
                                             interpolate(surge10, surge100, seg_FLOPROS, 10, 100),
                                             interpolate(surge100, surge1000, seg_FLOPROS, 100, 1000)
                                            )
                                   )
                          )                                            
    return seg_FLOPROH





#### Pre-processing of SLIIDERS

def create_FLOPROH_SLIIDERS(seg_ypcc, seg_elev_cap, seg_elev_pop, seg_length, params, surgeHeights):
    # WB groups by GDP per capita in USD2024
    # SLIIDERS has data in USD 2019 PPP
    inflation_factor = 1.0/1.23
    
    bound1 = 1135. * inflation_factor
    bound2 = 4495. * inflation_factor
    bound3 = 13935. * inflation_factor
    
    # FLOPROS is "maximum return period against which there will be protection"
    FLOPROS_min = 2.0 # minimum protection level 
    FLOPROS_max = 1000.0 # maximum protection level
    
    seg_pop = np.sum(seg_elev_pop, axis=1)
    mask = np.where(seg_ypcc<=bound1, 1, np.where(seg_ypcc<=bound2, 2, np.where(seg_ypcc<=bound3, 3, 4)))

    mean_ypcc = np.zeros(4)
    for i in range(4):
        mean_ypcc[i] = np.sum(np.where(mask==i+1, seg_ypcc*seg_pop, 0))/np.sum(np.where(mask==i+1, seg_pop, 0))
    
    grouped_F_mins = (mean_ypcc/mean_ypcc[0]) * FLOPROS_min
    grouped_F_maxs = (mean_ypcc/mean_ypcc[-1]) * FLOPROS_max
    
    # the expected annually exposed asses (or capital) per segment
    # length (EAEA_length) are used here instead of EAD_area in FLOPROS.
    # For this we calculate the annual exposure 
    seg_EAEA_length = calc_exposure_SLIIDERS(0, 0, surgeHeights[:,-1], seg_elev_cap, params, scalar=True) / seg_length
    
    
    # containers for the *old* minima / maxima 
    
    old_mins = np.empty(4)
    old_maxs = np.empty(4)
    seg_FLOPROS = np.zeros_like(seg_ypcc)
    for i in range(4):
        # values that belong to the current group
        group_vals = np.where(mask == i + 1, seg_EAEA_length, np.nan)

        old_mins[i] = np.nanmin(group_vals)
        old_maxs[i] = np.nanmax(group_vals)


        interpolated = interpolate(grouped_F_mins[i],
                                  grouped_F_maxs[i],
                                  group_vals,          # ← already masked
                                  old_mins[i],
                                  old_maxs[i])


        seg_FLOPROS = np.where(mask == i + 1,
                               interpolated,
                               seg_FLOPROS)

    

    # calculate population per coastline length susceptible
    # to a 1 in 1000 year flood, to serve as a limit
    susceptible_popdens = calc_flooding_SLIIDERS(surgeHeights[:,2], 0, seg_elev_pop) / seg_length

        
    # Turning the FLOPROS into a flood protection height (FLOPROH) 
    # by calculating the respective height of a storm surge with the
    # return period of FLOPROS.
    seg_FLOPROH = np.where(susceptible_popdens < 1, 0, gumbel_r.ppf(1 - 1.0/seg_FLOPROS, params[:,0], params[:,1])) 
    
    
    return seg_FLOPROH    

def calc_flooding_SLIIDERS(sl, H, C, delta=0.1):
    # This calculates the flooded amount of quantity C for a given
    # sea level sl and protection height H using a simply bathtub approach
    
    # Some modifications if sl or H are given as scalars not as segment-level arrays
    N = len(C[:,0])          
    if np.isscalar(sl): sl = np.full(N,sl)
    if np.isscalar(H): H = np.full(N,H)
    
    # The maximum index that could be flooded by a sea level sl
    iMax = np.minimum(199, (sl/delta)).astype(int)
    rows = np.arange(C.shape[0])
    
    flooded_C = np.cumsum(C,axis=1)[rows,iMax]
    
    # if sl is higher than H, add remaining fractional stuff not included above
    # Only add up those places where sl > H.
    flooded_C = np.where(sl>H, flooded_C + C[rows,iMax]*(sl/delta - (sl/delta).astype(int)), 0.0)
    
    return flooded_C


def calc_exposure_SLIIDERS(lsl, H, s10000, C, params, elevs=np.linspace(0.05,19.95,200), delta=0.1,
                           pop=False, scalar=False):
    
    N = len(elevs)
    total_exposure = np.zeros(len(C[:,0]))
    
    for i in range(0,N):
        sh = (i+1)*delta # reference height of the surge 
        sl = sh + lsl # actual height of the surge
        iexp = sl/delta # omit the -1 here, because would then just be used as ":iexp+1"
        iexp = np.asarray(iexp, dtype=int)
                
        probability = np.where(sh > s10000, 0, np.where(sl <= H, 0, gumbel_r.pdf(sh, params[:,0], params[:,1])))
        
        elev_idx = np.arange(len(elevs))
        if not scalar:
            if pop: exposure = np.sum(np.where(elev_idx < iexp[:,np.newaxis], C*ddf_P(sl[:,np.newaxis]-elevs[np.newaxis,:]), 0.0), axis=1)
            else:   exposure = np.sum(np.where(elev_idx < iexp[:,np.newaxis], C*ddf_C(sl[:,np.newaxis]-elevs[np.newaxis,:]), 0.0), axis=1)            
        else:
            if pop: exposure = np.sum(C[:,:iexp]*ddf_P(sl-elevs[np.newaxis,:iexp]), axis=1)
            else:   exposure = np.sum(C[:,:iexp]*ddf_C(sl-elevs[np.newaxis,:iexp]), axis=1)
        
        total_exposure += delta * probability * exposure
    
    return total_exposure






#### Depth damage functions

def ddf_C(d):
    # Depth damage function for capital as in CIAM.
    # Even though we only calculate exposure using this function, we still apply the 
    # depth damage function already now, to have some accounting of the depth of the flooding,
    # before aggregating the information.
    return np.where(d<=0, 0, d/(1+d))

def ddf_P_CIAM(d):
    # Fatalitiy rate in CIAM is just a constant fraction
    # Here, we are only interested whether people are exposed or not.
    return np.where(d<=0, 0, 1)

def ddf_P(d):
    # In contrast to CIAM let us here use a linear increase from 0 to 1 for the first 1 m,
    # to account for reduced damages to people at very shallow floods.
    return np.where(d<=0, 0, np.minimum(d,1))



#### Other functions

def extend_with_decreasing_growth(x, horizon, ndim=3, istart=90):
    # For extending input SSP data beyond 2100
    slope = np.asarray([x[istart] - x[istart-1]])

    years = np.linspace(1,100,100) # from 2101 to 2200
    for i in range(ndim): years = years[:,np.newaxis]        
    
    rslope = slope[np.newaxis,...]*np.maximum(1.0 - years/horizon, 0.0)
    
    for i in range(istart+1,len(x)):
        x[i] = x[i-1] + rslope[i-istart-1]
    return x


def interpolate(New_min, New_max, x, Old_xmin, Old_xmax):
    return New_min + (New_max-New_min)*(x-Old_xmin)/(Old_xmax-Old_xmin)


def logistic_function(x,k=1.0,x0=0.0, A=1.0):
    return A / (1.0 + np.exp(-k*(x-x0)))

def logistic_function2(x,k=1.0,x0=0.0, A=1.0, c=0.0):
    return A / (1.0 + np.exp(-k*(x-x0))) + c

def log_function(x,a,b):
    return a*np.log(x*b+1)

def log_function2(x,a,b,c):
    return a*np.log(x*b+1)+c

