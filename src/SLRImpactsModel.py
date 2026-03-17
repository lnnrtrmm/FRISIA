# © Lennart Ramme, MPI-M, 2024-2026

import numpy as np
import sys
import pandas as pd
import warnings

warnings.filterwarnings('error', category=RuntimeWarning)

SLIIDERS_versions = ['SLIIDERS_global',
                     'SLIIDERS_regional',
                     'SLIIDERS_CapitalDens',
                     'SLIIDERS_PopDens',
                     'SLIIDERS_PopDens100_lim0']

class SLRImpactModel:
    '''
    This is the impacts and adaptation module of FRISIA

    Inputs:
    - global SLR or SLR components -> will be translated into coastal/regional population weighted SLR
    - global mean surface temperature anomaly -> for estimating the expected SLR in 50 years
    - total global CO2 emissions -> for estimating the expected SLR in 50 years
    - asset time series -> for reference growth
    - population time series -> for reference growth
    - GDP time series -> for reference growth
    - Information on aggregation level


    Abbreviatios used:
    - SLR: sea level rise
    - fp: flood protection
    - sy: start year
    - ey: end year
    - dt: length of time step
    - dbg: debugging option
    - nreg: number of regions
    
    '''



    def __init__(self, SLR_total, T, CO2emis, population, gdp, assets=[],
                 nreg=2, sy=2010, ey=2100, dt=1.0, dbg=0, damage=True, 
                 include_SLR_components=False, SLR_components=[], randomize=False,
                 USDyear=2010, version='SLIIDERS_global', input_path='../input/'):


        if version in SLIIDERS_versions: self.database = 'SLIIDERS'
        else: sys.exit('Given version is not defined! '+version)
        self.version = version


        # Timestep:
        self.sy = sy
        self.ey = ey
        self.dt = dt
        # We assume that the final year ey, is included in the data
        self.nyears = int((self.ey - self.sy) / self.dt) + 1
        self.time = np.linspace(self.sy, self.ey, self.nyears)

        self.dbg = dbg
        self.nreg = nreg
        self.randomize = randomize

        #### Set some global variables from input
        # All input has to start in the specified start year!
        self.SLR_total = SLR_total
        self.SLR_components = np.asarray(SLR_components)
        self.include_SLR_components = include_SLR_components
        self.T_anomaly = np.copy(T) # anomaly wrt 1750 in Kelvin
        self.CO2emis = np.copy(CO2emis) # GT C per year


        #### Some checks and initial processing
        # GDP and population time series have to be region-specific!
        # -> First dimension is region (nreg) and second dimension is time (nyears)
        # -> These are only used for getting the (theoretical) internal growth rates
        if population.shape != (self.nreg, self.nyears) or gdp.shape != (self.nreg, self.nyears):
            sys.exit('Population or GDP arrays require ndim=2 (nreg,ntime)!')
        self.population = np.copy(population) # Mio. people
        self.GDP = np.copy(gdp) # billion USD per year

        # Same for asset timeseries, if it is provided
        if len(assets) == 0:
            self.assets = self.GDP * 3.0 # If assets are not given, assume constant asset/GDP ratio of 3 as in CIAM
        else:
            if assets.shape != (self.nreg, self.nyears):
                sys.exit('Asset array requires ndim=2 (nreg,ntime)!')
            self.assets = np.copy(assets)

        if self.database == 'SLIIDERS':
            # SLIIDERS data are in USD2019 PPP, but here we assume that all values are USD2010 initially
            # (prior to conversion with USD_fac). Hence, rescale monetary input arrays with USD inflation rate 
            # between 2010 and 2019 (17.2 %).
            self.assets *= 1./1.172
            self.GDP *= 1./1.172




        ####################################################################################
        ###########          Main structural switches             ##########################
        ####################################################################################

        self.damage = damage

        self.willingness_to_invest_in_fp = 0.0 # [0,1] Can also be time-dependent array
        self.willingness_to_retreat = 0.0 # [0,1]

        self.asset_feedback_switch = 1.0 # [0,1] Do storm damages reduce the value of coastal assets?
        self.people_feedback_switch = 1.0 # [0,1] Do storm fatalities reduce the coastal population?

        # Will protection or retreat be against current SLR or against SLR in 50 years?
        self.include_foresight_in_adaptation = True

        # Include reduced investment in coastal zones in case there will be increased flood heights in the future?
        self.include_reduced_growth = True
        # Include the possibility that assets are moved from unsafe to safe coastal zones (if above param is True)
        self.move_around_growth = True

        # Include maximum available money for flood protection as fraction of GDP?
        self.include_fp_investment_cap = True

        # Include the possibility that raised flood protection is breached (i.e. protection is raised, but SLR is even faster)
        self.include_failing_protection = True

        # Include the effect of asset reduction and changing population on GDP per capita?
        self.include_gdp_effect = True

        # Include the effect of retreat leading to reduced storm surge exposure?
        self.include_retreat_exposure_reduction = True

        # Include productivity losses from annually flooded people in GDP calculation?
        self.include_productivity_feedback = True




        ####################################################################################
        ###########          Uncertainty Parameters          ###############################
        ####################################################################################
        # If self.randomize = True, the model will randomize these parameters within the given 
        # range and using a uniform distribution. Otherwise, the default values are used,
        # which are based on CIAM values (if possible).
        #
        # Watch out: every new uncertainty parameter has to also be added to the randomisation routine
        #            and the input & output routines of the uncertainty parameters

        #### Coastal population and asset evolution parameters
        # Flood fatalities as feedback from storm surges to population stock and for output
        self.flood_event_fatality_rate = 0.01                  # dmnl
        self.flood_event_fatality_rate_range = (0.005, 0.02)

        # Flood damages that feed back to asset stock
        self.fraction_of_storm_damages_that_is_repaired = 0.9                # dmnl
        self.fraction_of_storm_damages_that_is_repaired_range = (0.75, 1.0)

        # The fraction of exposed assets that is going to be damaged (after applying resilience) (=1 in CIAM ?!)
        # This is used to calibrate to CIAM output
        self.flood_event_damage_fraction = 0.3         # dmnl
        self.flood_event_damage_fraction_range = (0.2, 0.4)

        # Determines how fast susceptibility is reduced under retreat
        # A value of 1 means that all generally susceptible assets/people are equally likely to retreat.
        # A higher value means that those that are more likely to be exposed will preferably retreat.
        self.susceptibility_reduction_exponent = 1.5
        self.susceptibility_reduction_exponent_range = (1.0, 2.0)   # dmnl

        # Elasticity of coastal GDP per capita with respect to coastal asset intensity
        self.gdp_asset_elasticity = 0.35
        self.gdp_asset_elasticity_range = (0.2, 0.6)

        # Relative productivity of fully flood-affected workers compared to unaffected workers.
        # 1.0 means no productivity loss; values in [0.9, 0.99] imply a 1-10% reduction.
        self.relative_worker_productivity_if_flooded = 0.96
        self.relative_worker_productivity_if_flooded_range = (0.9, 0.99)

        # Calibration parameter for reduced investment in unprotected coastal zones
        self.coastal_asset_depreciation_time_scale = 30.
        self.coastal_asset_depreciation_time_scale_range = (20.,40.)
        self.effective_flood_height_at_which_investment_is_halved = 1.0                   # m
        self.effective_flood_height_at_which_investment_is_halved_range = (0.5, 3.0)
        self.safe_coastal_zone_likelihood_threshold = 0.95                          # dmnl
        self.safe_coastal_zone_likelihood_threshold_range = (0.9, 1.0)
        self.fraction_of_investments_that_must_be_at_the_coast = 0.5                # dmnl
        self.fraction_of_investments_that_must_be_at_the_coast_range = (0.2, 0.8)

        #### Coastal protection parameters
        self.maximum_gdp_fraction_for_fp_investment = 0.03                  # dmnl
        self.maximum_gdp_fraction_for_fp_investment_range = (0.01, 0.05)
        self.fp_construction_duration = 10.0                                # year
        self.fp_construction_duration_range = (5, 25)
        self.fp_construction_cost_reference_USD2010 = 0.00602                       # b$/m^2/km, from CIAM
        self.fp_construction_cost_reference_USD2010_range = (0.005, 0.007)
        self.maintenance_cost_fraction = 0.02                                             # 1/ year, from CIAM
        self.maintenance_cost_fraction_range = (0.015, 0.03)
        self.coastal_land_value_init_USD2010 = 0.005376                             # b$/km^2, from CIAM
        self.coastal_land_value_init_USD2010_range = (0.005, 0.006)
        self.land_opportunity_cost_rate = 0.04                              # 1/year, from CIAM
        self.land_opportunity_cost_rate_range = (0.03, 0.05)
        self.maximum_fp_deterioration_rate = 0.05                           # m/year
        self.maximum_fp_deterioration_rate_range = (0.02, 0.1)


        #### Retreat parameters; default from CIAM with added uncertainty ranges
        self.asset_relocation_cost_factor = 0.1                                         # dmnl
        self.asset_relocation_cost_factor_range = (0.05, 0.15)
        self.asset_demolition_cost_factor = 0.05                                        # dmnl
        self.asset_demolition_cost_factor_range = (0.025, 0.075)
        self.not_depreciated_fraction_of_assets_at_time_of_retreat = 0.1                # dmnl
        self.not_depreciated_fraction_of_assets_at_time_of_retreat_range = (0.0, 0.2)
        self.people_retreat_cost_factor = 8.0                                           # dmnl
        self.people_retreat_cost_factor_range = (3.0, 10.9)                             # values as in DSCIM-Coastal, but extended to the lower end values of other studies
        self.proactive_retreat_time_scale = 10.                                         # year
        self.proactive_retreat_time_scale_range = (5.0, 25.0)
        self.retreat_sensitivity = 0.5                                                  # dmnl
        self.retreat_sensitivity_range = (0,1)                                          # 0: retreat when reached by mean regional sea level, 1: retreat when reached by 1-in-1 year return period height

        ########## Fixed parameters ######################

        ### Expected SLR in 50 years
        # Sensitivity parameters for calculation of expected SLR in 50 years.
        # Parameter values are from SLR model analysis.
        self.expSLR_sens_CO2emis = 0.0131  # m per Gt C
        self.expSLR_sens_Tano = 0.091      # m per K
        self.expSLR_timeHorizon = 50.0

        # Calculating the construction cost index (cci) as a linear function fitted to CIAM data, bounded between 0.5 and 2.5
        self.fp_cci_fit_slope = 0.036853
        self.fp_cci_fit_intercept = 0.3786
        self.fp_cci_min = 0.5
        self.fp_cci_max = 2.5

        # GDP per capita in the US in 2010 (reference for calculating damage resilience)
        self.ypc_US_2010_USD2010 = 54.41  # b$/Mp/year

        # constant GDP parameters for now
        self.max_gdp_pc_asset_factor = 3.0

        # USD conversion
        # Conversion in case that USD values should be for different year than 2010;
        # all internal, default dollar values are in USD 2010.
        # Provide the conversion factors for individual years. 
        # USD_fac can be overwritten before the integration in case other years are requested.
        if USDyear==2010:
            self.USD_fac = 1.0
        elif USDyear==2005:
            self.USD_fac = 0.8957
        elif USDyear==2017:
            self.USD_fac = 1.1241
        elif USDyear==2019:
            self.USD_fac = 1.172
        elif USDyear==2022:
            self.USD_fac = 1.3421
        else:
            sys.exit(str(USDyear)+' not available as USD value year!')

        ################################################################################
        ##########         Initialize arrays         ###################################
        ################################################################################
        self.SLR                      = np.zeros((self.nreg, self.nyears))
        self.expected_SLR_in_50_years = np.zeros((self.nreg, self.nyears))

        ### Flood protection (fp) variables
        self.average_fp_height                                = np.zeros((self.nreg, self.nyears))
        self.annual_increase_in_fp_height_from_investment     = np.zeros((self.nreg, self.nyears))
        self.annual_reduction_in_fp_height_from_deterioration = np.zeros((self.nreg, self.nyears))
        self.construction_cost                                = np.zeros((self.nreg, self.nyears))
        self.effective_annual_investment_in_fp                = np.zeros((self.nreg, self.nyears))
        self.potential_fp_height_increase_over_50_years       = np.zeros((self.nreg, self.nyears))
        self.fp_land_opportunity_cost                         = np.zeros((self.nreg, self.nyears))
        self.annual_costs_of_fp_maintenance                   = np.zeros((self.nreg, self.nyears))
        self.annual_SLR_costs_of_fp_maintenance               = np.zeros((self.nreg, self.nyears))
        self.mask_annual_fp_investment_limited                = np.zeros((self.nreg, self.nyears))


        ### Population variables
        self.annual_people_flooded       = np.zeros((self.nreg, self.nyears))
        self.annual_flood_fatalities     = np.zeros((self.nreg, self.nyears))
        self.orig_susceptible_people_fraction = np.zeros((self.nreg, self.nyears))
        self.orig_exposed_people_fraction     = np.zeros((self.nreg, self.nyears))
        self.inundated_original_people_fraction   = np.zeros((self.nreg, self.nyears))
        self.retreated_original_people_fraction   = np.zeros((self.nreg, self.nyears))
        # People retreat
        self.annual_reactive_people_retreat  = np.zeros((self.nreg, self.nyears))
        self.annual_proactive_people_retreat = np.zeros((self.nreg, self.nyears))
        self.annual_total_people_retreat     = np.zeros((self.nreg, self.nyears))
        self.people_retreat_cost             = np.zeros((self.nreg, self.nyears))
        
        
        ### Asset variables
        self.annual_storm_damage_to_assets             = np.zeros((self.nreg, self.nyears))
        self.storm_damage_resilience                   = np.zeros((self.nreg, self.nyears))
        self.landvalue_appreciation_factor             = np.ones((self.nreg, self.nyears))
        self.likelihood_of_investment_in_coastal_zones = np.zeros((self.nreg, self.nyears))
        self.orig_susceptible_asset_fraction           = np.zeros((self.nreg, self.nyears))
        self.orig_exposed_asset_fraction               = np.zeros((self.nreg, self.nyears))
        self.inundated_original_asset_fraction         = np.zeros((self.nreg, self.nyears))
        self.retreated_original_asset_fraction         = np.zeros((self.nreg, self.nyears))

        # Asset retreat
        self.annual_reactive_asset_retreat  = np.zeros((self.nreg, self.nyears))
        self.annual_proactive_asset_retreat = np.zeros((self.nreg, self.nyears))
        self.annual_total_asset_retreat     = np.zeros((self.nreg, self.nyears))
        self.asset_demolition_cost          = np.zeros((self.nreg, self.nyears))
        self.asset_relocation_cost          = np.zeros((self.nreg, self.nyears))
        self.assets_lost_during_retreat     = np.zeros((self.nreg, self.nyears))
        self.total_annual_cost_of_retreat   = np.zeros((self.nreg, self.nyears))

        # Loss of area
        self.inundated_area = np.zeros((self.nreg, self.nyears))
        self.abandoned_area = np.zeros((self.nreg, self.nyears))
        self.effective_retreat_height = np.zeros((self.nreg, self.nyears))

        # net flood height
        # This is the net of global SLR and the building of flood protection
        self.effective_flood_height        = np.zeros((self.nreg, self.nyears))


        # Variables for quantities that are flooded on average once per year
        # This is for calculating reactive retreat under no adaptation
        self.surge1_inundated_area                     = np.zeros((self.nreg, self.nyears))
        self.surge1_inundated_original_asset_fraction  = np.zeros((self.nreg, self.nyears))
        self.surge1_inundated_original_people_fraction = np.zeros((self.nreg, self.nyears))

        # GDP productivity feedback diagnostics
        self.annual_GDP_loss_from_productivity_effect = np.zeros((self.nreg, self.nyears))
        ################################################################################
        #### Initial values and fit parameters for pre-defined model versions ##########
        ################################################################################
        
        # These are predefined aggregation levels that load input parameters from files
        # Input parameters are general aggregated information from SLIIDERS database 
        # and fit parameters used to match aggregated infomation for inundation,
        # susceptibility and surge exposure etc.

        ### Loading the general aggregated information
        path = input_path+'aggregated_data/'
        df = pd.read_csv(path+self.version+'_information.csv', delimiter=',')
        
        # The factor multiplied with SLR_total to get regional SLR. Only used in case the individual components are not given
        # so that regional SLR has to be fitted linearly to total global SLR
        self.SLR_factor              = df.SLR_factor.values[:]
        # The weights of the individual SLR components for each region to get regional SLR (in case this is activated).
        # Not included are thermosteric and LWS SLR, because weights are 1 in BRICK.
        self.SLR_weight_MG           = df.SLR_weight_MG.values[:]   
        self.SLR_weight_GIS          = df.SLR_weight_GIS.values[:]  
        self.SLR_weight_AIS          = df.SLR_weight_AIS.values[:]
        
        # Information on the coastal regions
        self.average_fp_height_init  = df.average_fp_height.values[:]
        self.total_fp_length         = df.total_fp_length.values[:]
        self.mobile_asset_fraction   = df.mobcapfrac.values[:,np.newaxis]

        # Information on coastal surge heights
        self.surge_heights = np.ones((self.nreg,4))
        self.surge_heights[:,0] = df.surge1_height.values[:]
        self.surge_heights[:,1] = df.surge10_height.values[:]
        self.surge_heights[:,2] = df.surge100_height.values[:]
        self.surge_heights[:,3] = df.surge1000_height.values[:]

        ### Loading in the fit parameters
        # 2 cases (parameters with and without initial flood protection), 4 parameters and X regions
        self.inund_params_area            = np.zeros((2,4,self.nreg))
        self.inund_params_assets          = np.zeros((2,4,self.nreg))
        self.inund_params_people          = np.zeros((2,4,self.nreg))
        self.surge1_inund_params_area     = np.zeros((2,4,self.nreg))
        self.surge1_inund_params_assets   = np.zeros((2,4,self.nreg))
        self.surge1_inund_params_people   = np.zeros((2,4,self.nreg))
        self.storm_suscept_params_assets  = np.zeros((2,4,self.nreg))
        self.storm_suscept_params_people  = np.zeros((2,4,self.nreg))
        self.storm_exposure_params_assets = np.zeros((2,4,self.nreg))
        self.storm_exposure_params_people = np.zeros((2,4,self.nreg))
        
        self.variables = [self.inundated_area, self.inundated_original_asset_fraction, self.inundated_original_people_fraction,
                          self.surge1_inundated_area, self.surge1_inundated_original_asset_fraction, self.surge1_inundated_original_people_fraction,
                          self.orig_susceptible_asset_fraction, self.orig_susceptible_people_fraction,
                          self.orig_exposed_asset_fraction, self.orig_exposed_people_fraction]
        self.parameters = [self.inund_params_area, self.inund_params_assets, self.inund_params_people,
                           self.surge1_inund_params_area, self.surge1_inund_params_assets, self.surge1_inund_params_people,
                           self.storm_suscept_params_assets, self.storm_suscept_params_people,
                           self.storm_exposure_params_assets, self.storm_exposure_params_people]
        self.limits = [None, 1.0, 1.0, None, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        filenames_extensions = ['inund_params_area', 'inund_params_assets', 'inund_params_people',
                                'ann_inund_params_area', 'ann_inund_params_assets', 'ann_inund_params_people',
                                'storm_suscept_params_assets', 'storm_suscept_params_people',
                                'storm_exposure_params_assets', 'storm_exposure_params_people']

        path = input_path+'fit_function_parameters/'
        for iparam, parameter in enumerate(self.parameters):
            # Fits including the initial flood protection
            # Here, the data were fitted using a logistic function: y = A / (1 + exp(-k*(x-x0)))
            df = pd.read_csv(path+self.version+'_'+filenames_extensions[iparam]+'.csv', delimiter=',')
            parameter[0,0,:] = df.k.values[:]
            parameter[0,1,:] = df.x0.values[:]
            parameter[0,2,:] = df.A.values[:]
            parameter[0,3,:] = df.c.values[:]
            
            # Fits without initial flood protection, same function as above
            df = pd.read_csv(path+self.version+'_'+filenames_extensions[iparam]+'_no_initial_dikes.csv', delimiter=',')
            parameter[1,0,:] = df.k.values[:]
            parameter[1,1,:] = df.x0.values[:]
            parameter[1,2,:] = df.A.values[:]
            parameter[1,3,:] = df.c.values[:]

        ################################################################################
        ################################################################################
        ################################################################################


    def __safe_divide(self, numerator, denominator, default=0.0):
        out = np.full(numerator.shape, default, dtype=float)
        np.divide(numerator, denominator, out=out, where=denominator > 0.0)
        return out

    def __cobb_douglas_output(self, tfp, capital, labor, alpha):
        capital = np.maximum(np.asarray(capital, dtype=float), 0.0)
        labor = np.maximum(np.asarray(labor, dtype=float), 0.0)
        return tfp * np.power(capital, alpha) * np.power(labor, 1.0 - alpha)

    def fit_function(self, x, params, limit=1.0, ifunc=0):
        k = params[ifunc,0]
        x0 = params[ifunc,1]
        A = params[ifunc,2]
        c = params[ifunc,3]
        result = np.maximum(1e-6,A/(1.0 + np.exp(-k*(x-x0))) + c)
        if limit is not None:
            result = np.minimum(limit, result)
        return result

    def __fitted_variable(self, dflood, dh, params, limit=1.0):

        # dh >= 0:
        # This accounts for the possibility that dikes are raised, but SLR was faster (include_failing_protection=True)
        # E.g. for inundated asset fractions:
        # Assume a SLR of 2m, but protection was also increased by dh=1m. Then dflood = SLR - dh = 1m.
        # The fit_function would simply calculate the inundated asset fraction "f0" for 1m in this case.
        # However, in reality the value would be higher, because those segments, where the protection is breached would be
        # inundated by 2m of SLR not 1m. Inundated asset fractions of 2m of SLR are f1 here.
        # The answer is somewhere in the middle, hence the approximated equation below, which was tested on the segment 
        # data to produce reasonable estimates of the inundated asset fractions (and the other parameters) in these scenarios.

        # If this is not activated, then dh = 0 and f = f0, so the option is automatically turned off.
        f0_pos = self.fit_function(dflood, params, limit=limit)
        f1_pos = self.fit_function(dflood+dh, params, limit=limit)

        # dh < 0:
        # If the average dike height has been reduced, then we should not use dflood itself to calculate the damages.
        # Instead we can start by defining the minimum and maximum inundated/susceptible/exposed fractions:
        # The minimum fraction is that of only using SLR (SLR = dflood + dh), assuming dikes would still be as in the initial state
        # The maximum fraction is derived from using the same SLR, but assuming that all dikes are gone
        f0_neg = self.fit_function(dflood+dh, params, limit=limit)
        f1_neg = self.fit_function(dflood+dh, params, limit=limit, ifunc=1)
        

        return np.where(dh >= 0, np.where(dflood < 0, f0_pos, f0_pos + (f1_pos-f0_pos) * (f0_pos/(params[0,2,:] + params[0,3,:]))),
                                 f0_neg + (f1_neg - f0_neg) * (np.abs(dh) / self.average_fp_height[:,0]) )

    def __calc_autonomous_retreat_fraction(self,i, variable):
        if variable == 'area':
            lower_bound = self.inundated_area[:,i]
            upper_bound = self.surge1_inundated_area[:,i]
        elif variable == 'people':
            lower_bound = self.inundated_original_people_fraction[:,i]
            upper_bound = self.surge1_inundated_original_people_fraction[:,i]
        elif variable == 'assets':
            lower_bound = self.inundated_original_asset_fraction[:,i]
            upper_bound = self.surge1_inundated_original_asset_fraction[:,i]
        else: sys.exit('Calling retreat with wrong variable (must be area, people or assets):', variable)
        return lower_bound + self.retreat_sensitivity * (upper_bound - lower_bound)


    def __update_single_param(self, prange, scale_factor):
        return prange[0] + scale_factor * (prange[1] - prange[0])

    def __randomize_parameters(self):
        self.flood_event_fatality_rate                             = self. __update_single_param(self.flood_event_fatality_rate_range, np.random.rand())
        self.fraction_of_storm_damages_that_is_repaired            = self. __update_single_param(self.fraction_of_storm_damages_that_is_repaired_range, np.random.rand())
        self.flood_event_damage_fraction                           = self. __update_single_param(self.flood_event_damage_fraction_range, np.random.rand())
        self.coastal_asset_depreciation_time_scale                 = self. __update_single_param(self.coastal_asset_depreciation_time_scale_range, np.random.rand())
        self.effective_flood_height_at_which_investment_is_halved  = self. __update_single_param(self.effective_flood_height_at_which_investment_is_halved_range, np.random.rand())
        self.safe_coastal_zone_likelihood_threshold                = self. __update_single_param(self.safe_coastal_zone_likelihood_threshold_range, np.random.rand())
        self.fraction_of_investments_that_must_be_at_the_coast     = self. __update_single_param(self.fraction_of_investments_that_must_be_at_the_coast_range, np.random.rand())
        self.maximum_gdp_fraction_for_fp_investment                = self. __update_single_param(self.maximum_gdp_fraction_for_fp_investment_range, np.random.rand())
        self.fp_construction_duration                              = self. __update_single_param(self.fp_construction_duration_range, np.random.rand())
        self.fp_construction_cost_reference_USD2010                = self. __update_single_param(self.fp_construction_cost_reference_USD2010_range, np.random.rand())
        self.maintenance_cost_fraction                             = self. __update_single_param(self.maintenance_cost_fraction_range, np.random.rand())
        self.coastal_land_value_init_USD2010                       = self. __update_single_param(self.coastal_land_value_init_USD2010_range, np.random.rand())
        self.land_opportunity_cost_rate                            = self. __update_single_param(self.land_opportunity_cost_rate_range, np.random.rand())
        self.asset_relocation_cost_factor                          = self. __update_single_param(self.asset_relocation_cost_factor_range, np.random.rand())
        self.asset_demolition_cost_factor                          = self. __update_single_param(self.asset_demolition_cost_factor_range, np.random.rand())
        self.not_depreciated_fraction_of_assets_at_time_of_retreat = self. __update_single_param(self.not_depreciated_fraction_of_assets_at_time_of_retreat_range, np.random.rand())
        self.people_retreat_cost_factor                            = self. __update_single_param(self.people_retreat_cost_factor_range, np.random.rand())
        self.proactive_retreat_time_scale                          = self. __update_single_param(self.proactive_retreat_time_scale_range, np.random.rand())
        self.susceptibility_reduction_exponent                     = self. __update_single_param(self.susceptibility_reduction_exponent_range, np.random.rand())
        self.gdp_asset_elasticity                                  = self. __update_single_param(self.gdp_asset_elasticity_range, np.random.rand())
        self.relative_worker_productivity_if_flooded               = self. __update_single_param(self.relative_worker_productivity_if_flooded_range, np.random.rand())
        self.maximum_fp_deterioration_rate                         = self. __update_single_param(self.maximum_fp_deterioration_rate_range, np.random.rand())
        self.retreat_sensitivity                                   = self. __update_single_param(self.retreat_sensitivity_range, np.random.rand())
        return

    def getUncertaintyParameters(self):
        uncertaintyParameters = [
                np.copy(self.flood_event_fatality_rate),
                np.copy(self.fraction_of_storm_damages_that_is_repaired),
                np.copy(self.flood_event_damage_fraction),
                np.copy(self.coastal_asset_depreciation_time_scale),
                np.copy(self.effective_flood_height_at_which_investment_is_halved),
                np.copy(self.safe_coastal_zone_likelihood_threshold),
                np.copy(self.fraction_of_investments_that_must_be_at_the_coast),
                np.copy(self.maximum_gdp_fraction_for_fp_investment),
                np.copy(self.fp_construction_duration),
                np.copy(self.fp_construction_cost_reference_USD2010),
                np.copy(self.maintenance_cost_fraction),
                np.copy(self.coastal_land_value_init_USD2010),
                np.copy(self.land_opportunity_cost_rate),
                np.copy(self.asset_relocation_cost_factor),
                np.copy(self.asset_demolition_cost_factor),
                np.copy(self.not_depreciated_fraction_of_assets_at_time_of_retreat),
                np.copy(self.people_retreat_cost_factor),
                np.copy(self.proactive_retreat_time_scale),
                np.copy(self.susceptibility_reduction_exponent),
                np.copy(self.gdp_asset_elasticity),
                np.copy(self.relative_worker_productivity_if_flooded),
                np.copy(self.maximum_fp_deterioration_rate),
                np.copy(self.retreat_sensitivity),
            ]
        return uncertaintyParameters

    def setUncertaintyParameters(self, InputParameters):
        self.flood_event_fatality_rate                              = InputParameters[0]
        self.fraction_of_storm_damages_that_is_repaired             = InputParameters[1]
        self.flood_event_damage_fraction                            = InputParameters[2]
        self.coastal_asset_depreciation_time_scale                  = InputParameters[3]
        self.effective_flood_height_at_which_investment_is_halved   = InputParameters[4]
        self.safe_coastal_zone_likelihood_threshold                 = InputParameters[5]
        self.fraction_of_investments_that_must_be_at_the_coast      = InputParameters[6]
        self.maximum_gdp_fraction_for_fp_investment                 = InputParameters[7]
        self.fp_construction_duration                               = InputParameters[8]
        self.fp_construction_cost_reference_USD2010                 = InputParameters[9]
        self.maintenance_cost_fraction                              = InputParameters[10]
        self.coastal_land_value_init_USD2010                        = InputParameters[11]
        self.land_opportunity_cost_rate                             = InputParameters[12]
        self.asset_relocation_cost_factor                           = InputParameters[13]
        self.asset_demolition_cost_factor                           = InputParameters[14]
        self.not_depreciated_fraction_of_assets_at_time_of_retreat  = InputParameters[15]
        self.people_retreat_cost_factor                             = InputParameters[16]
        self.proactive_retreat_time_scale                           = InputParameters[17]
        self.susceptibility_reduction_exponent                      = InputParameters[18]
        self.gdp_asset_elasticity                                   = InputParameters[19]
        self.relative_worker_productivity_if_flooded                = InputParameters[20]
        self.maximum_fp_deterioration_rate                          = InputParameters[21]
        self.retreat_sensitivity                                    = InputParameters[22]

        return




    def __initialise_variables(self):

        if self.dbg == 1: print('Initialising variables!')

        # First, randomize the uncertainty parameters, if this is activated.
        if self.randomize: self.__randomize_parameters()

        ### The (regional/segregated/global but wheighted) SLR is different from the calculated total SLR!
        if self.include_SLR_components:
            self.SLR_thermo = self.SLR_components[0,np.newaxis,:] * np.ones(self.nreg)[:,np.newaxis]
            self.SLR_LWS    = self.SLR_components[1,np.newaxis,:] * np.ones(self.nreg)[:,np.newaxis]
            self.SLR_MG     = self.SLR_components[2,np.newaxis,:] * self.SLR_weight_MG[:,np.newaxis]
            self.SLR_GIS    = self.SLR_components[3,np.newaxis,:] * self.SLR_weight_GIS[:,np.newaxis]
            self.SLR_AIS    = self.SLR_components[4,np.newaxis,:] * self.SLR_weight_AIS[:,np.newaxis]
            self.SLR = self.SLR_thermo + self.SLR_LWS + self.SLR_MG + self.SLR_GIS + self.SLR_AIS
        else:
            self.SLR = self.SLR_total[np.newaxis,:] *  self.SLR_factor[:,np.newaxis]

        self.effective_flood_height[:,0] = 0.0 # Assume that net flood height is always 0 in the first time step.

        
        # Expected SLR can be calculated already now, because T and CO2 emission time series are external anyway.
        self.expected_SLR_in_50_years[:,:] = (self.T_anomaly[np.newaxis,:] * self.expSLR_sens_Tano \
                                             + self.CO2emis[np.newaxis,:] * self.expSLR_sens_CO2emis) * self.SLR_factor[:,np.newaxis]

        # Convert dollar variables to USD of a specific year (default: USD 2010)
        self.ypc_US_2010 = self.ypc_US_2010_USD2010 * self.USD_fac
        self.coastal_land_value_init = self.coastal_land_value_init_USD2010 * self.USD_fac
        self.fp_construction_cost_reference = self.fp_construction_cost_reference_USD2010 * self.USD_fac

        # Copy the input references arrays into the variables that we track
        self.coastal_population = np.copy(self.population)
        self.coastal_assets = np.copy(self.assets) * self.USD_fac
        self.coastal_GDP = np.copy(self.GDP) * self.USD_fac
        self.coastal_GDPperCapita = self.__safe_divide(self.coastal_GDP, self.coastal_population)

        # Calculating flood damage resilience as in CIAM, depending on GDP per capita
        self.storm_damage_resilience[:,0] = self.coastal_GDPperCapita[:,0] / (self.coastal_GDPperCapita[:,0] + self.ypc_US_2010)

        # Initialise main stocks
        self.average_fp_height[:,0]  = self.average_fp_height_init[:]
        

        # Now calculate the initial fractions of assets and people that are exposed to storm surges generally, as well as the
        # initial values for the inundated asset, people and area fractions. The latter should theoretically be 0 with the original data,
        # but the fitted function produces very small values for a SLR of 0 meters. 
        # Also, if SLR is not 0 in the initial year, then this will also be > 0.
    
        for iparam, params in enumerate(self.parameters):
            self.variables[iparam][:,0] = self.__fitted_variable(self.effective_flood_height[:,0], 0, params, limit=self.limits[iparam])

        # Check if given willingness_to_invest variable is time dependent variable with correct length or a scalar
        if not isinstance(self.willingness_to_invest_in_fp, (list, np.ndarray)):
            self.willingness_to_invest_in_fp = np.zeros((self.nreg, self.nyears))+self.willingness_to_invest_in_fp
        elif len(self.willingness_to_invest_in_fp[0,:]) != self.nyears or  len(self.willingness_to_invest_in_fp[:,0]) != self.nreg:
            sys.exit('Length of willingness to invest array does not fit with time or aggregation dimension')
        # Check if given willingness_to_retreat variable is time dependent variable with correct length or a scalar
        if not isinstance(self.willingness_to_retreat, (list, np.ndarray)):
            self.willingness_to_retreat = np.zeros((self.nreg, self.nyears))+self.willingness_to_retreat
        elif len(self.willingness_to_retreat[0,:]) != self.nyears or  len(self.willingness_to_retreat[:,0]) != self.nreg:
            sys.exit('Length of willingness to retreat array does not fit with time or aggregation dimension')


    ################################################
    # MAIN INTEGRATION LOOP
    def integrate(self, silent=True):
        if not silent: print('Start integrating...')

        # INITIALISATION
        self.__initialise_variables()

        # TIME LOOP
        for i in range(0,self.nyears):
            if self.dbg==1: print('Year: ', self.time[i])

            self.__update_flood_protection(i)
            self.__update_population(i)
            self.__update_assets(i)
            if i<self.nyears-1: self.__prepare_next_timestep(i)

        # OUTPUT
        self.__make_outputs()

        return

    ################################################


    def __update_flood_protection(self, i):
        if self.dbg==1: print('    update flood protection')

        # Calculate construction cost
        cci = np.maximum(self.fp_cci_min, np.minimum(self.fp_cci_max , 
                       self.fp_cci_fit_slope * self.coastal_GDPperCapita[:,i] + self.fp_cci_fit_intercept))
        self.construction_cost[:,i] = self.fp_construction_cost_reference * cci


        ### MAINTENANCE COSTS
        self.annual_costs_of_fp_maintenance[:,i] = self.average_fp_height[:,i] * self.total_fp_length \
                                                   * self.construction_cost[:,i] *  self.maintenance_cost_fraction
        ref_maintenance_costs = self.average_fp_height[:,0] * self.total_fp_length \
                                                   * self.construction_cost[:,i] *  self.maintenance_cost_fraction
        
        self.annual_SLR_costs_of_fp_maintenance[:,i] = np.maximum(0, self.annual_costs_of_fp_maintenance[:,i] - ref_maintenance_costs)

        ### FLOOD PROTECTION INVESTMENT
        # Current net flood height + expected SLR is missing protection
        if self.include_foresight_in_adaptation: missing_protection = self.effective_flood_height[:,i] + self.expected_SLR_in_50_years[:,i]
        else: missing_protection = self.effective_flood_height[:,i]
        desired_protection = missing_protection * self.willingness_to_invest_in_fp[:,i]

        cost_of_reaching_desired_protection = np.maximum(0, self.construction_cost[:,i] * self.total_fp_length * \
                    (( desired_protection + self.average_fp_height[:,i])**2 - self.average_fp_height[:,i]**2) )
    
        maintenance_cost = self.annual_costs_of_fp_maintenance[:,i]
        money_available_for_fp = np.maximum(0, self.coastal_GDP[:,i] * self.maximum_gdp_fraction_for_fp_investment \
                                            - maintenance_cost)

            
        if self.include_fp_investment_cap:
            self.effective_annual_investment_in_fp[:,i] = np.minimum(cost_of_reaching_desired_protection / self.fp_construction_duration,
                                                                       money_available_for_fp)
            self.mask_annual_fp_investment_limited[:,i] = np.where(cost_of_reaching_desired_protection / self.fp_construction_duration > money_available_for_fp, 1.0, 0.0)
        else:
            self.effective_annual_investment_in_fp[:,i] = cost_of_reaching_desired_protection / self.fp_construction_duration

        # Calculate potential maximum fp height increase for asset investment decision making
        # Use current maximum spending and the expected willingness to invest in fp, based on previous willingness in the past
        istart = np.maximum(0,i-int(self.proactive_retreat_time_scale))
        expected_willingness_to_invest_in_fp = np.mean(self.willingness_to_invest_in_fp[:,istart:i+1], axis=1)
        potential_investment_over_50_years = money_available_for_fp * expected_willingness_to_invest_in_fp * self.expSLR_timeHorizon
        self.potential_fp_height_increase_over_50_years[:,i] =  np.sqrt( self.average_fp_height[:,i]**2 \
                        + potential_investment_over_50_years / (self.total_fp_length * self.construction_cost[:,i]) ) - self.average_fp_height[:,i]



        # change in fp height from investment
        self.annual_increase_in_fp_height_from_investment[:,i] = np.sqrt( self.average_fp_height[:,i]**2 + self.effective_annual_investment_in_fp[:,i] \
                    / (self.total_fp_length * self.construction_cost[:,i]) ) - self.average_fp_height[:,i]


        ### FLOOD PROTECTION DECLINE
        if self.include_failing_protection:
            # decline in average fp height, if it is below the 1-year return period surge height
            surge1_height = self.surge_heights[:,0] + (self.SLR[:,i] - self.SLR[:,0])
            # maximum deterioration rate if surge1 is more than 1 m higher than seawall, linear increase inbetween 
            annual_deterioration = self.maximum_fp_deterioration_rate * np.minimum(1.0, np.maximum(0.0, surge1_height - self.average_fp_height[:,i]))
            self.annual_reduction_in_fp_height_from_deterioration[:,i] = np.minimum(self.average_fp_height[:,i], annual_deterioration)
                                                                                   

        if i >= self.nyears-1: return
        ##########################################################################################
        ###          
        ### UPDATE THE FLOOD PROTECTION HEIGHT
        ###
        ##########################################################################################
        self.average_fp_height[:,i+1] = self.average_fp_height[:,i] + self.annual_increase_in_fp_height_from_investment[:,i] \
                                                                    - self.annual_reduction_in_fp_height_from_deterioration[:,i]
        
        return


    def __update_population(self, i):
        if self.dbg==1: print('    update population')

        if self.damage:
            ##########################################################################################
            ###
            ### RETREAT of people from the coast
            ###
            ##########################################################################################

            # There are two ways how people retreat from the coast:
            #  1. Reactive retreat in response to flooding.
            #  2. Proactive retreat in response to expected exposure and flooding.
        
            #### Reactive retreat ###################################################################
            # only after first timestep, because it is in response to SLR from one time step to the next
            if i==0: 
                self.annual_reactive_people_retreat[:,i] = 0.0
                total_removed_fraction = self.__calc_autonomous_retreat_fraction(0,'people')
            else:
                # First calculate the previous fraction of original people distribution that was already removed (inundated or retreated)
                total_removed_fraction = np.maximum(self.__calc_autonomous_retreat_fraction(i-1,'people'), self.retreated_original_people_fraction[:,i-1])
                total_removed_fraction = np.clip(total_removed_fraction, 0.0, 1.0 - 1e-12)
                self.annual_reactive_people_retreat[:,i] = self.coastal_population[:,i] \
                                                            * np.maximum(0.0, (self.__calc_autonomous_retreat_fraction(i,'people')\
                                                              - total_removed_fraction) / (1.0 - total_removed_fraction) )
                # update removed fraction for calculation of proactive retreat
                total_removed_fraction = np.maximum(self.__calc_autonomous_retreat_fraction(i,'people'), self.retreated_original_people_fraction[:,i-1])
                total_removed_fraction = np.clip(total_removed_fraction, 0.0, 1.0 - 1e-12)


            #### Proactive retreat #####################################################################
            # If this strategy is chosen, this will lead to additional retreat based on two quantities:
            #  1. If annual exposure will be higher than in the start year, this drives retreat.
            #  2. If people will become inundated, this will drive their retreat.
            #  -> In contrast to protection, which uses expected SLR in 50 years for estimating required investment,
            #     we use a shorter time horizon for retreat (the proactive_retreat_time_scale). We
            #     interpolate linearly between the current and the expected SLR and fp height in 50 years.
            if self.include_foresight_in_adaptation:
                # what is the flood height in the future?
                expected_effective_flood_height = self.effective_flood_height[:,i] + (self.expected_SLR_in_50_years[:,i] \
                                                - self.potential_fp_height_increase_over_50_years[:,i]) \
                                                * (self.proactive_retreat_time_scale / self.expSLR_timeHorizon)
            else:
                expected_effective_flood_height = self.effective_flood_height[:,i]

            # including the possibility that raised dikes are breached
            if self.include_failing_protection: 
                dh = self.average_fp_height[:,i] + self.potential_fp_height_increase_over_50_years[:,i] - self.average_fp_height[:,0]
            else: 
                dh = 0.0

            if self.include_retreat_exposure_reduction:
                # Update the susceptible fraction based on previous retreat
                # Storm damage to people is depending on what is the actually susceptible fraction of people. 
                # This has to account for the fraction of people that is generally susceptible and the fraction
                # of people that has already been removed from the coast (via inundation or retreat)
                expected_orig_sus_people_fraction = self.__fitted_variable(expected_effective_flood_height,
                                                                                   dh, self.storm_suscept_params_people)
        
                expected_actual_sus_people_fraction = np.maximum(0, (expected_orig_sus_people_fraction - total_removed_fraction) \
                                                              / (1.0 - total_removed_fraction))

                exposure_reduction = np.where(expected_orig_sus_people_fraction == 0, 0, 
                        (expected_actual_sus_people_fraction/expected_orig_sus_people_fraction)**self.susceptibility_reduction_exponent)
            else:
                exposure_reduction = 1.0

            # Expected exposure fraction of original distribution times the reduction factor to account for previous retreat
            expected_exposed_fraction = self.__fitted_variable(expected_effective_flood_height, dh, self.storm_exposure_params_people) \
                                        * exposure_reduction

            lower_bound = self.__fitted_variable(expected_effective_flood_height, dh, self.inund_params_people)
            upper_bound = self.__fitted_variable(expected_effective_flood_height, dh, self.surge1_inund_params_people)
            expected_inundated_fraction = np.minimum(1.0, lower_bound + self.retreat_sensitivity * (upper_bound - lower_bound))

            retreating_people_fraction = self.willingness_to_retreat[:,i] \
                    *  np.maximum(0, np.maximum(expected_exposed_fraction - self.orig_exposed_people_fraction[:,0],
                                                expected_inundated_fraction - total_removed_fraction))
            
            # Safeguard against retreating too many people
            retreating_people_fraction = np.minimum(retreating_people_fraction, 1.0 - total_removed_fraction)


            # Here also approximate the abandoned land to estimate opportunity costs later on
            lower_bound = self.__fitted_variable(expected_effective_flood_height, dh, self.inund_params_area, limit=None)
            upper_bound = self.__fitted_variable(expected_effective_flood_height, dh, self.surge1_inund_params_area, limit=None)
            expected_inundated_land = lower_bound + self.retreat_sensitivity * (upper_bound - lower_bound)
            current_inundated_land = self.__calc_autonomous_retreat_fraction(i,'area')
            self.abandoned_area[:,i] = current_inundated_land + self.willingness_to_retreat[:,i] * (expected_inundated_land - current_inundated_land) \

            if i==0: self.retreated_original_people_fraction[:,i] = self.retreated_original_people_fraction[:,0] + retreating_people_fraction
            else: self.retreated_original_people_fraction[:,i] = self.retreated_original_people_fraction[:,i-1] + retreating_people_fraction
    
            self.annual_proactive_people_retreat[:,i] = (self.coastal_population[:,i] - self.annual_reactive_people_retreat[:,i]) \
                                                         * retreating_people_fraction / (1.0 - total_removed_fraction)
    
            # Total asset retreat is sum of reactive and proactive retreat
            self.annual_total_people_retreat[:,i] = self.annual_reactive_people_retreat[:,i] +  self.annual_proactive_people_retreat[:,i]

            ##########################################################################################
            ####
            #### Number of people flooded in storms and flood fatalities
            ####
            ########################################################################################## 
            # Update the susceptible fraction again
            total_removed_fraction = np.maximum(self.__calc_autonomous_retreat_fraction(i,'people'), self.retreated_original_people_fraction[:,i])
            total_removed_fraction = np.clip(total_removed_fraction, 0.0, 1.0 - 1e-12)
            actually_susceptible_people_fraction = np.maximum(0, (self.orig_susceptible_people_fraction[:,i] - total_removed_fraction) \
                                                     / (1.0 - total_removed_fraction))

            if self.include_retreat_exposure_reduction:
                exposure_reduction = np.where(self.orig_susceptible_people_fraction[:,i] == 0, 0, 
                        (actually_susceptible_people_fraction/self.orig_susceptible_people_fraction[:,i])**self.susceptibility_reduction_exponent)
            else: exposure_reduction = 1.0       

            # Subtract here the initial exposure fraction from current exposure fraction to get the SLR driven number
            self.annual_people_flooded[:,i] = (self.coastal_population[:,i] - self.annual_total_people_retreat[:,i]) * np.maximum(0.0, 
                    self.orig_exposed_people_fraction[:,i] * exposure_reduction - self.orig_exposed_people_fraction[:,0])

            self.annual_flood_fatalities[:,i] = self.flood_event_fatality_rate * (1.0 - self.storm_damage_resilience[:,i]) \
                    * self.annual_people_flooded[:,i]


        if i >= self.nyears-1: return


        ##########################################################################################
        ####
        #### Updating the population variables:
        ####
        ##########################################################################################
        self.coastal_population[:,i+1] = self.coastal_population[:,i] * self.population[:,i+1]/self.population[:,i] \
                                         - self.annual_total_people_retreat[:,i] - self.people_feedback_switch*self.annual_flood_fatalities[:,i]

        return



    def __update_assets(self, i):
        if self.dbg==1: print('    update assets')

        # No retreat or storm damage in "No damage" reference case
        if self.damage:
            ##########################################################################################
            ###
            ### RETREAT of assets from the coast
            ###
            ##########################################################################################

            # There are two ways how assets retreat from the coast:
            #  1. Reactive retreat in response to flooding.
            #  2. Proactive retreat in response to expected exposure and flooding.

            #### Reactive retreat ####################################################################
            # only after first timestep, because it is in response to SLR from one time step to the next
            if i==0: 
                self.annual_reactive_asset_retreat[:,i] = 0.0
                total_removed_fraction = self.__calc_autonomous_retreat_fraction(0,'assets')
            else:
                # First calculate the previous fraction of original asset distribution that was already removed (inundated or retreated)
                total_removed_fraction = np.maximum(self.__calc_autonomous_retreat_fraction(i-1,'assets'), self.retreated_original_asset_fraction[:,i-1])
                total_removed_fraction = np.clip(total_removed_fraction, 0.0, 1.0 - 1e-12)
                self.annual_reactive_asset_retreat[:,i] = self.coastal_assets[:,i] \
                                                            * np.maximum(0.0, (self.__calc_autonomous_retreat_fraction(i,'assets')\
                                                              - total_removed_fraction) / (1.0 - total_removed_fraction) )
                # update removed fraction for calculation of proactive retreat
                total_removed_fraction = np.maximum(self.__calc_autonomous_retreat_fraction(i,'assets'), self.retreated_original_asset_fraction[:,i-1])
                total_removed_fraction = np.clip(total_removed_fraction, 0.0, 1.0 - 1e-12)

            
            if self.dbg==1: print('        reactive asset retreat done')

            #### Proactive retreat #####################################################################
            # If this strategy is chosen, this will lead to additional retreat based on two quantities:
            #  1. If annual exposure will be higher than in the start year, this drives retreat.
            #  2. If people will become inundated, this will drive their retreat.
            #  -> In contrast to protection, which uses expected SLR in 50 years for estimating required investment,
            #     we use a shorter time horizon for retreat (the proactive_retreat_time_scale). We
            #     interpolate linearly between the current and the expected SLR and fp height in 50 years.
            if self.include_foresight_in_adaptation:
                # what is the flood height in the future?
                expected_effective_flood_height = self.effective_flood_height[:,i] + (self.expected_SLR_in_50_years[:,i] \
                                                - self.potential_fp_height_increase_over_50_years[:,i]) \
                                                * (self.proactive_retreat_time_scale / self.expSLR_timeHorizon)
            else:
                expected_effective_flood_height = self.effective_flood_height[:,i]

            # including the possibility that raised dikes are breached
            if self.include_failing_protection: 
                dh = self.average_fp_height[:,i] + self.potential_fp_height_increase_over_50_years[:,i] - self.average_fp_height[:,0]
            else: 
                dh = 0.0


            if self.include_retreat_exposure_reduction:
                # Update the susceptible fraction based on previous retreat
                # Storm damage to people is depending on what is the actually susceptible fraction of people. 
                # This has to account for the fraction of people that is generally susceptible and the fraction
                # of people that has already been removed from the coast (via inundation or retreat)
                expected_orig_sus_asset_fraction = self.__fitted_variable(expected_effective_flood_height,
                                                                                   dh, self.storm_suscept_params_assets)
        
                expected_actual_sus_asset_fraction = np.maximum(0, (expected_orig_sus_asset_fraction - total_removed_fraction) \
                                                              / (1.0 - total_removed_fraction))

                exposure_reduction = np.where(expected_orig_sus_asset_fraction == 0, 0, 
                        (expected_actual_sus_asset_fraction/expected_orig_sus_asset_fraction)**self.susceptibility_reduction_exponent)
            else: exposure_reduction = 1.0



            # Expected exposure fraction of original distribution times the reduction factor to account for previous retreat
            expected_exposed_fraction = self.__fitted_variable(expected_effective_flood_height, dh, self.storm_exposure_params_assets) \
                                        * exposure_reduction

            lower_bound = self.__fitted_variable(expected_effective_flood_height, dh, self.inund_params_assets)
            upper_bound = self.__fitted_variable(expected_effective_flood_height, dh, self.surge1_inund_params_assets)
            expected_inundated_fraction = np.minimum(1.0, lower_bound + self.retreat_sensitivity * (upper_bound - lower_bound))

            retreating_asset_fraction = self.willingness_to_retreat[:,i] \
                                          *  np.maximum(0, np.maximum(expected_exposed_fraction - self.orig_exposed_asset_fraction[:,0],
                                                                      expected_inundated_fraction - total_removed_fraction))

            # Safeguard against retreating too many assets
            retreating_asset_fraction = np.minimum(retreating_asset_fraction, 1.0 - total_removed_fraction)

            if i==0: self.retreated_original_asset_fraction[:,i] = self.retreated_original_asset_fraction[:,0] + retreating_asset_fraction
            else: self.retreated_original_asset_fraction[:,i] = self.retreated_original_asset_fraction[:,i-1] + retreating_asset_fraction
    
            self.annual_proactive_asset_retreat[:,i] = (self.coastal_assets[:,i] - self.annual_reactive_asset_retreat[:,i]) * retreating_asset_fraction / (1.0 - total_removed_fraction)
    
            # Total asset retreat is sum of reactive and proactive retreat
            self.annual_total_asset_retreat[:,i] = self.annual_reactive_asset_retreat[:,i] + self.annual_proactive_asset_retreat[:,i]

            if self.dbg==1: print('        proactive retreat done')

            ##########################################################################################
            ###
            ### STORM DAMAGE to assets
            ###
            ##########################################################################################
        
            # Storm damage to assets is depending on what is the actually susceptible fraction of assets. This has to account for the fraction of
            # assets that is generally susceptible and the fraction of assets that has already been removed from the coast
            # (via inundation or retreat)
            total_removed_fraction = np.maximum(self.__calc_autonomous_retreat_fraction(i,'assets'), self.retreated_original_asset_fraction[:,i])
            total_removed_fraction = np.clip(total_removed_fraction, 0.0, 1.0 - 1e-12)
            actually_susceptible_asset_fraction = np.maximum(0, (self.orig_susceptible_asset_fraction[:,i] - total_removed_fraction)\
                                                                  / (1.0 - total_removed_fraction))
            if self.include_retreat_exposure_reduction:
                exposure_reduction = np.where(self.orig_susceptible_asset_fraction[:,i] == 0, 0, 
                        (actually_susceptible_asset_fraction/self.orig_susceptible_asset_fraction[:,i])**self.susceptibility_reduction_exponent)
            else: exposure_reduction = 1.0
            
            # Subtract here the initial exposure fraction from current exposure fraction to get the SLR driven number
            self.annual_storm_damage_to_assets[:,i] = (self.coastal_assets[:,i] - self.annual_total_asset_retreat[:,i]) \
                                                       *self.flood_event_damage_fraction * (1.0 - self.storm_damage_resilience[:,i]) * np.maximum(0.0, 
                                                         self.orig_exposed_asset_fraction[:,i] * exposure_reduction - self.orig_exposed_asset_fraction[:,0])



            if self.dbg==1: print('        damage calculation done')
        else:
            total_removed_fraction = 0.0
        ##########################################################################################
        ###
        ### GROWTH of assets.
        ###
        ##########################################################################################

        # Calculate likelihood of investment in coastal zone
        # -> Less growth if it is expected that net flood height will increase (or is already high)
        expected_effective_flood_height = self.effective_flood_height[:,i] + self.expected_SLR_in_50_years[:,i] \
                                            - self.potential_fp_height_increase_over_50_years[:,i]
        if self.include_failing_protection:  dh = self.average_fp_height[:,i] + self.potential_fp_height_increase_over_50_years[:,i] \
                                                     - self.average_fp_height[:,0]
        else: dh = 0
        expected_orig_susceptible_fraction = self.__fitted_variable(expected_effective_flood_height, dh, self.storm_suscept_params_assets)
        expected_actual_susceptible_fraction = (expected_orig_susceptible_fraction - total_removed_fraction) \
                                                / (1.0 - total_removed_fraction)

        if self.include_reduced_growth and self.damage:
            dflood =  np.maximum(0, self.expected_SLR_in_50_years[:,i] + self.effective_flood_height[:,i]\
                                     - self.potential_fp_height_increase_over_50_years[:,i])
            self.likelihood_of_investment_in_coastal_zones[:,i] = (1.0 - np.maximum(0, dflood / (dflood \
                                            + self.effective_flood_height_at_which_investment_is_halved ))) * expected_actual_susceptible_fraction \
                                            + (1.0 - expected_actual_susceptible_fraction)
        else:
            self.likelihood_of_investment_in_coastal_zones[:,i] = 1.0

        # No need to update next years asset value, if we are already in the last time step...
        if i >= self.nyears-1: return      



        # First calculate reference asset growth and investment in each coastal zone
        reference_asset_growth = self.coastal_assets[:,i] * (self.assets[:,i+1] / self.assets[:,i] - 1.0)
        asset_depreciation = self.coastal_assets[:,i] / self.coastal_asset_depreciation_time_scale
        asset_investment = reference_asset_growth + asset_depreciation

        # initialising some arrays for investment redistribution
        actual_asset_investment = np.copy(asset_investment)
        safe_zone = np.ones_like(asset_investment)
        investment_moved_to_safe_zones = np.zeros(self.nreg)


        # Remove small part of the investment from the insufficiently protected assets
        for icz, likelihood in enumerate(self.likelihood_of_investment_in_coastal_zones[:,i]):
            if likelihood < self.safe_coastal_zone_likelihood_threshold:
                actual_asset_investment[icz] = asset_investment[icz] * likelihood
                investment_moved_to_safe_zones[icz] = asset_investment[icz] * (1.0 - likelihood)
                safe_zone[icz] = 0
        

        # Get the fraction of each coastal zone with likelihood > threshold to all coastal zones with likelihood > threshold
        # avoid div by 0, in case there is no safe coastal zone
        if np.sum(safe_zone) != 0.0:
            safe_asset_fractions = (self.coastal_assets[:,i]*safe_zone) / np.sum(self.coastal_assets[:,i]*safe_zone)
        
        # Add some of the removed investment to well protected assets,
        # but not all of it! Some is moved away from the coast, because it is not coast-specific.
        investment_staying_at_coast = investment_moved_to_safe_zones * self.fraction_of_investments_that_must_be_at_the_coast
        investment_moving_away_from_coast = investment_moved_to_safe_zones - investment_staying_at_coast

        if np.sum(safe_zone) == 0 or not self.move_around_growth:
            actual_asset_investment += investment_staying_at_coast
        else:
            total_investment_staying_at_coast = np.sum(investment_staying_at_coast)
            for icz, likelihood in enumerate(self.likelihood_of_investment_in_coastal_zones[:,i]):
                if likelihood >= self.safe_coastal_zone_likelihood_threshold:
                    actual_asset_investment[icz] += total_investment_staying_at_coast * safe_asset_fractions[icz]

        if self.dbg==1: print('        asset growth calculation done')        
        ##########################################################################################
        ###
        ### UPDATE Coastal asset stock
        ###
        ##########################################################################################
        # The asset feedback switch (default: 0) determines what fraction of storm damages actually reduces coastal asset values.
        # In CIAM this is 0, because it is assumed that storm surge damages are repaired always.
        self.coastal_assets[:,i+1] = self.coastal_assets[:,i] + actual_asset_investment - asset_depreciation \
                                     - self.annual_total_asset_retreat[:,i] - self.asset_feedback_switch \
                                   * (1.0 - self.fraction_of_storm_damages_that_is_repaired) * self.annual_storm_damage_to_assets[:,i]
        return
        


    def __prepare_next_timestep(self,i):
        
        # Calculate GDP(perCapita) internally via Cobb-Douglas formulation or leave as external input
        # Cobb-Douglas formulation: GDP = TFP * K^alpha * L^(1-alpha)
        if self.include_gdp_effect:

            ref_gdp = self.GDP[:,i+1] * self.USD_fac
            ref_capital = self.assets[:,i+1] * self.USD_fac
            ref_labor = self.population[:,i+1]

            ref_effective_inputs = np.power(np.maximum(ref_capital, 0.0), self.gdp_asset_elasticity) \
                                 * np.power(np.maximum(ref_labor, 0.0), 1.0 - self.gdp_asset_elasticity)

            ref_tfp = self.__safe_divide(ref_gdp, ref_effective_inputs, default=0.0)

            # Baseline GDP without productivity feedback (for diagnostics and counterfactual loss)
            gdp_without_productivity_feedback = self.__cobb_douglas_output(
                ref_tfp,
                self.coastal_assets[:,i+1],
                self.coastal_population[:,i+1],
                self.gdp_asset_elasticity
            )

            if self.include_productivity_feedback:
                flooded_people_fraction = self.__safe_divide(
                    self.annual_people_flooded[:,i],
                    self.coastal_population[:,i],
                    default=0.0
                )
                productivity_multiplier = 1.0 - (1.0 - self.relative_worker_productivity_if_flooded) * flooded_people_fraction
            else:
                productivity_multiplier = np.ones(self.nreg)

            self.coastal_GDP[:,i+1] = gdp_without_productivity_feedback * productivity_multiplier
            self.annual_GDP_loss_from_productivity_effect[:,i+1] = np.maximum(
                0.0,
                gdp_without_productivity_feedback - self.coastal_GDP[:,i+1]
            )

            self.coastal_GDPperCapita[:,i+1] = self.__safe_divide(
                self.coastal_GDP[:,i+1],
                self.coastal_population[:,i+1],
                default=0.0
            )
        else:
            self.coastal_GDP[:,i+1] = self.GDP[:,i+1] * self.USD_fac
            self.annual_GDP_loss_from_productivity_effect[:,i+1] = 0.0
            self.coastal_GDPperCapita[:,i+1] = self.__safe_divide(
                self.coastal_GDP[:,i+1],
                self.coastal_population[:,i+1],
                default=0.0
            )

        # Calculating flood damage resilience as in CIAM, depending on GDP per capita
        self.storm_damage_resilience[:,i+1] = self.coastal_GDPperCapita[:,i+1] / (self.coastal_GDPperCapita[:,i+1] + self.ypc_US_2010)
        
        gdp_pc_ratio = self.__safe_divide(
            self.coastal_GDPperCapita[:,i+1],
            self.coastal_GDPperCapita[:,i],
            default=1.0
        )
        pop_ratio = self.__safe_divide(
            self.coastal_population[:,i+1],
            self.coastal_population[:,i],
            default=1.0
        )

        #### Updating the value of land using a formula from CIAM repository
        landvalue_growth = 0.565 * (gdp_pc_ratio - 1.0) + 0.313 * (pop_ratio - 1.0)
        landvalue_growth = np.clip(landvalue_growth, -0.1, 0.1)

        self.landvalue_appreciation_factor[:,i+1] = self.landvalue_appreciation_factor[:,i] * (1.0 + landvalue_growth)


        # The net global flood height 
        self.effective_flood_height[:,i+1] = (self.SLR[:,i+1]-self.SLR[:,0]) - (self.average_fp_height[:,i+1] - self.average_fp_height[:,0])
        
        # Update the susceptible and inundated theoretical fractions for the next time step, based on the updated net flood height
        dflood = self.effective_flood_height[:,i+1]
        if self.include_failing_protection:  dh = self.average_fp_height[:,i+1] - self.average_fp_height[:,0]
        else: dh = 0

        for iparam, param in enumerate(self.parameters):
            self.variables[iparam][:,i+1] = self.__fitted_variable(dflood,dh,param, limit=self.limits[iparam])

        return

    def __make_outputs(self):
        
        # ASSETS

        self.assets_lost_during_reactive_retreat = (1.0 - self.mobile_asset_fraction) * self.annual_reactive_asset_retreat
        self.assets_lost_during_proactive_retreat = (1.0 - self.mobile_asset_fraction) * self.annual_proactive_asset_retreat * self.not_depreciated_fraction_of_assets_at_time_of_retreat
        self.assets_lost_during_retreat = self.assets_lost_during_reactive_retreat + self.assets_lost_during_proactive_retreat

        self.asset_relocation_cost_reactive = self.annual_reactive_asset_retreat * self.mobile_asset_fraction * self.asset_relocation_cost_factor
        self.asset_relocation_cost_proactive = self.annual_proactive_asset_retreat * self.mobile_asset_fraction * self.asset_relocation_cost_factor
        self.asset_relocation_cost = self.annual_total_asset_retreat * self.mobile_asset_fraction * self.asset_relocation_cost_factor

        self.asset_demolition_cost_reactive = self.annual_reactive_asset_retreat * self.asset_demolition_cost_factor * (1.0 - self.mobile_asset_fraction)
        self.asset_demolition_cost_proactive = self.annual_proactive_asset_retreat * self.asset_demolition_cost_factor * (1.0 - self.mobile_asset_fraction)
        self.asset_demolition_cost = self.annual_total_asset_retreat * self.asset_demolition_cost_factor * (1.0 - self.mobile_asset_fraction)

        # PEOPLE
        
        self.people_retreat_cost_reactive = self.annual_reactive_people_retreat * self.coastal_GDPperCapita * self.people_retreat_cost_factor
        self.people_retreat_cost_proactive = self.annual_proactive_people_retreat * self.coastal_GDPperCapita * self.people_retreat_cost_factor
        self.people_retreat_cost = self.people_retreat_cost_reactive + self.people_retreat_cost_proactive


        # CIAM reference: This is for comparison to CIAM output

        # LAND
        # Opportunity costs
        # The opportunity cost is assumed to be 4% of the landvalue annually
        # First calculate the value of land per square kilometer
        land_cost_per_sqkm = self.landvalue_appreciation_factor * self.coastal_land_value_init * self.land_opportunity_cost_rate

        # 1. Opportunity cost: value of land used up for flood protection
        #   - DON'T (subtract the initial fp height, in order to only get the additional opportunity costs), not done in CIAM
        #     -> From now on, do actually only account for the increase in fp height, despite CIAM not doing it
        #   - Factor 1.7 comes from ratio between dike height and width (60° angle, as in CIAM)
        #   - In the end divided by 2, as in CIAM
        self.fp_land_opportunity_cost = self.total_fp_length[:,np.newaxis] * 1.7 * self.average_fp_height * 0.001 * land_cost_per_sqkm / 2.0

        # 2. Opportunity cost: value of land that is inundated or abandoned
        lower_bound = self.inundated_area
        upper_bound = self.surge1_inundated_area
        inundated_area = lower_bound + self.retreat_sensitivity * (upper_bound - lower_bound)
        self.inundated_area_opportunity_cost = inundated_area * land_cost_per_sqkm
        self.abandoned_area_opportunity_cost = (self.abandoned_area - inundated_area) * land_cost_per_sqkm

        # Total opportunity costs of the land that is lost
        self.total_opportunity_cost = self.fp_land_opportunity_cost + self.abandoned_area_opportunity_cost + self.inundated_area_opportunity_cost

        # All land that is lost either for protection or from inundation or from retreat
        self.land_lost = np.maximum(self.abandoned_area, self.total_fp_length[:,np.newaxis] * 1.7 * self.average_fp_height * 0.001 / 2.0)

        #####
        # CIAM reference: These are for comparison to CIAM output
        #####

        ### SLR driven maintenance costs
        ref_maintenance_costs = self.average_fp_height[:,0,None] * self.total_fp_length[:,None] \
                                                   * self.construction_cost[:,:] *  self.maintenance_cost_fraction
        
        annual_SLR_costs_of_fp_maintenance = np.maximum(0, self.annual_costs_of_fp_maintenance - ref_maintenance_costs)


        self.relocation_cost = self.asset_relocation_cost + self.asset_demolition_cost + self.people_retreat_cost
        self.flood_cost = self.assets_lost_during_retreat + self.abandoned_area_opportunity_cost + self.inundated_area_opportunity_cost
        self.construct_cost = self.fp_land_opportunity_cost + self.effective_annual_investment_in_fp + annual_SLR_costs_of_fp_maintenance
        self.storm_cost = self.annual_storm_damage_to_assets


        ### Overall SLR costs, all potential costs added upp
        self.annual_SLR_costs = self.relocation_cost + self.flood_cost +  self.construct_cost + self.storm_cost

        ### SLR costs separated by damages or adaptation_costs
        self.annual_SLR_damages = self.assets_lost_during_reactive_retreat + self.asset_relocation_cost_reactive + self.asset_demolition_cost_reactive\
                                            + self.people_retreat_cost_reactive + self.storm_cost + self.inundated_area_opportunity_cost

        self.annual_SLR_adaptation_costs = self.assets_lost_during_proactive_retreat + self.asset_relocation_cost_proactive + self.asset_demolition_cost_proactive \
                                            + self.people_retreat_cost_proactive + self.construct_cost + self.abandoned_area_opportunity_cost


        return

    def getTotalCoastalGDP(self): return np.sum(self.coastal_GDP,axis=0)
    def getTotalGDPperCapita(self):
        return self.__safe_divide(
            np.sum(self.coastal_GDP, axis=0),
            np.sum(self.coastal_population, axis=0),
            default=0.0
        )
        
    def getTotalCoastalAssets(self): return np.sum(self.coastal_assets, axis=0)
    def getTotalCoastalPopulation(self): return np.sum(self.coastal_population, axis=0)

    def getTotalStormCost(self): return np.sum(self.storm_cost, axis=0)
    def getTotalRelocationCost(self): return np.sum(self.relocation_cost, axis=0)
    def getTotalFloodCost(self): return np.sum(self.flood_cost, axis=0)
    def getTotalConstructCost(self): return np.sum(self.construct_cost, axis=0)


    def getTotalProtectionInvestment(self): return np.sum(self.effective_annual_investment_in_fp, axis=0)
    def getTotalMaintenanceCost(self): return np.sum(self.annual_costs_of_fp_maintenance, axis=0)
    def getTotalOpportunityCost(self): return np.sum(self.total_opportunity_cost, axis=0)

    def getTotalPeopleFlooded(self): return np.sum(self.annual_people_flooded, axis=0)
    def getTotalFloodFatalities(self): return np.sum(self.annual_flood_fatalities, axis=0)
    def getTotalAnnualGDPLossFromProductivityEffect(self): return np.sum(self.annual_GDP_loss_from_productivity_effect, axis=0)
    





        





