# © Lennart Ramme, MPI-M, 2024

import numpy as np
import sys
import pandas as pd

class SLRImpactModel:
    '''
    This is the impacts and adaptation module of FRISIA

    Inputs:
    - global SLR or SLR components -> will be translated into coastal/regional population weighted SLR
    - global mean surface temperature anomaly -> for estimating the expected SLR in 50 years
    - total global CO2 emissions -> for estimating the expected SLR in 50 years
    - GDP time series -> for reference growth
    - population time series -> for reference growth
    - Information on aggregation level (currently available: DIVA_global, DIVA_bipolar, DIVA_regional; see below)


    Abbreviatios used:
    - SLR: sea level rise
    - fp: flood protection
    - sy: start year
    - ey: end year
    - dt: length of time step
    - dbg: debugging option
    - ndim: number of regions, i.e. 1 for global, 2 for bipolar and 7 for DIVA_regional
    
    '''



    def __init__(self, SLR_total, T, CO2emis, population, gdp, ndim=2, sy=2010, ey=2100, dt=1.0, dbg=0, damage=True, 
                 include_SLR_components=False, SLR_components=[], randomize=False,
                 USDyear=2010, version='DIVA_global', include_initial_fp=True, input_path='../input/'):

        # Timestep:
        self.sy = sy
        self.ey = ey
        self.dt = dt
        self.nyears = int((self.ey - self.sy) / self.dt) + 1
        self.time = np.linspace(self.sy+self.dt*0.5, self.ey-self.dt*0.5, self.nyears)

        self.dbg = dbg
        self.ndim = ndim
        self.randomize = randomize

        #### Set some global variables from input
        # All input has to start in the specified start year!
        self.SLR_total = SLR_total
        self.SLR_components = np.asarray(SLR_components)
        self.include_SLR_components = include_SLR_components
        self.T_anomaly = np.copy(T) # anomaly wrt 1750 in Kelvin
        self.CO2emis = np.copy(CO2emis) # GT C per year

        # GDP and population time series have to be region-specific!
        # -> First dimension is region (ndim) and second dimension is time (nyears)
        # -> These are only used for getting the (theoretical) internal growth rates
        if population.ndim != 2 or gdp.ndim != 2:
            sys.exit('Population or GDP arrays require ndim=2 (first dimension is number of coastal regions, second dimension is time)!')
        self.population = np.copy(population) # Mio. people
        self.GDP = np.copy(gdp) # billion USD per year




        ####################################################################################
        ###########          Main structural switches             ##########################
        ####################################################################################

        self.damage = damage
        self.include_initial_fp = include_initial_fp
        # This determines the type of fit that is used
        if self.include_initial_fp: self.iFit = 0
        else: self.iFit = 1

        self.willingness_to_invest_in_fp = 0.0 # [0,1] Can also be time-dependent array
        self.willingness_to_retreat = 0.0 # [0,1]

        self.asset_feedback_switch = 0.0 # [0,1] Do storm damages reduce the value of coastal assets?
        self.people_feedback_switch = 0.0 # [0,1] Do storm fatalities reduce the coastal population?

        ### Switches to make model reproduce CIAM output (if all following switches == False)
        # Include reduced investment in coastal zones in case there will be increased flood heights in the future?
        self.include_reduced_growth = False
        # Include the possibility that assets are moved from unsafe to safe coastal zones (if above param is True)
        self.move_around_growth = True

        # Include maximum available money for flood protection as fraction of GDP?
        self.include_fp_investment_cap = False

        # Include the effect of asset reduction and changing population on GDP per capita?
        self.include_gdp_effect = False

        # Include the effect of retreat leading to reduced storm surge exposure?
        self.include_retreat_exposure_reduction = False


        ####################################################################################
        ###########          Uncertainty Parameters          ###############################
        ####################################################################################
        # If self.randomize = True, the model will randomize these parameters within the given 
        # range and using a uniform distribution. Otherwise, the default values are used,
        # which are based on CIAM values (if possible).
        #
        # Watch out: every new uncertainty parameter has to also be added to the randomisation routine

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

        # Calibration parameter for reduced investment in unprotected coastal zones
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
        self.fp_construction_cost_reference = 0.00602                       # b$/m^2/km, from CIAM
        self.fp_construction_cost_reference_range = (0.005, 0.007)
        self.maintenance_cost_fraction = 0.02                                             # 1/ year, from CIAM
        self.maintenance_cost_fraction_range = (0.015, 0.03)
        self.coastal_land_value_init = 0.005376                             # b$/km^2, from CIAM
        self.coastal_land_value_init_range = (0.005, 0.006)
        self.land_opportunity_cost_rate = 0.04                              # 1/year, from CIAM
        self.land_opportunity_cost_rate_range = (0.03, 0.05)


        #### Retreat parameters; default from CIAM with added uncertainty ranges
        self.mobile_asset_fraction = 0.25                                               # dmnl
        self.mobile_asset_fraction_range = (0.2, 0.3)
        self.asset_relocation_cost_factor = 0.1                                         # dmnl
        self.asset_relocation_cost_factor_range = (0.05, 0.15)
        self.asset_demolition_cost_factor = 0.05                                        # dmnl
        self.asset_demolition_cost_factor_range = (0.025, 0.075)
        self.not_depreciated_fraction_of_assets_at_time_of_retreat = 0.1                # dmnl
        self.not_depreciated_fraction_of_assets_at_time_of_retreat_range = (0.0, 0.2)
        self.increase_factor_for_costs_of_reactive_retreat = 4.0                        # dmnl
        self.increase_factor_for_costs_of_reactive_retreat_range = (3.0, 5.0)
        self.proactive_retreat_time_scale = 10.                                         # year
        self.proactive_retreat_time_scale_range = (5.0, 25.0)       

        ########## Fixed parameters ######################

        ### Expected SLR in 50 years
        # Sensitivity parameters for calculation of expected SLR in 50 years.
        # Parameter values are from SLR model analysis.
        self.expSLR_sens_CO2emis = 0.0131  # m per Gt C
        self.expSLR_sens_Tano = 0.091      # m per K

        # Calculating the construction cost index (cci) as a linear function fitted to CIAM data, bounded between 0.5 and 2.5
        self.fp_cci_fit_slope = 0.036853
        self.fp_cci_fit_intercept = 0.3786
        self.fp_cci_min = 0.5
        self.fp_cci_max = 2.5

        # GDP per capita in the US in 2010 (reference for calculating damage resilience)
        self.ypc_US_2010 = 54.41  # b$/Mp/year


        # USD conversion
        # Conversion in case that USD values should be for different year than 2010;
        # all internal, default dollar values should be USD 2010.
        # Provide the conversion factors for individual years. 
        # USD_fac can be overwritten before the integration in case other years are requested.
        if USDyear==2010:
            self.USD_fac = 1.0
        elif USDyear==2017:
            self.USD_fac = 1.1241
        elif USDyear==2005:
            self.USD_fac = 0.8957
        elif USDyear==2022:
            self.USD_fac = 1.3421
        else:
            sys.exit(str(USDyear)+' not available as USD value year!')

        ################################################################################
        ##########         Initialize arrays         ###################################
        ################################################################################
        self.SLR                      = np.zeros((self.ndim, self.nyears))
        self.expected_SLR_in_50_years = np.zeros((self.ndim, self.nyears))

        self.coastal_population       = np.zeros((self.ndim, self.nyears))
        self.coastal_assets           = np.zeros((self.ndim, self.nyears))
        self.coastal_GDP              = np.zeros((self.ndim, self.nyears))
        self.coastal_GDPperCapita     = np.zeros((self.ndim, self.nyears))
        
        ### Flood protection (fp) variables
        self.average_fp_height                          = np.zeros((self.ndim, self.nyears))
        self.annual_change_in_fp_height                 = np.zeros((self.ndim, self.nyears))
        self.construction_cost                          = np.zeros((self.ndim, self.nyears))
        self.effective_annual_investment_in_fp          = np.zeros((self.ndim, self.nyears))
        self.potential_fp_height_increase_over_50_years = np.zeros((self.ndim, self.nyears))
        self.fp_land_opportunity_cost                   = np.zeros((self.ndim, self.nyears))
        self.annual_costs_of_fp_maintenance             = np.zeros((self.ndim, self.nyears))
        self.annual_costs_of_fp_maintenance_noadapt     = np.zeros((self.ndim, self.nyears))


        ### Population variables
        self.annual_people_flooded       = np.zeros((self.ndim, self.nyears))
        self.annual_flood_fatalities     = np.zeros((self.ndim, self.nyears))
        self.orig_susceptible_people_fraction = np.zeros((self.ndim, self.nyears))
        self.orig_exposed_people_fraction     = np.zeros((self.ndim, self.nyears))
        self.inundated_original_people_fraction   = np.zeros((self.ndim, self.nyears))
        self.retreated_original_people_fraction   = np.zeros((self.ndim, self.nyears))
        # People retreat
        self.annual_reactive_people_retreat  = np.zeros((self.ndim, self.nyears))
        self.annual_proactive_people_retreat = np.zeros((self.ndim, self.nyears))
        self.annual_total_people_retreat     = np.zeros((self.ndim, self.nyears))
        self.people_retreat_cost             = np.zeros((self.ndim, self.nyears))
        
        
        ### Asset variables
        self.annual_storm_damage_to_assets             = np.zeros((self.ndim, self.nyears))
        self.storm_damage_resilience                   = np.zeros((self.ndim, self.nyears))
        self.landvalue_appreciation_factor             = np.ones((self.ndim, self.nyears))
        self.likelihood_of_investment_in_coastal_zones = np.zeros((self.ndim, self.nyears))
        self.orig_susceptible_asset_fraction           = np.zeros((self.ndim, self.nyears))
        self.orig_exposed_asset_fraction               = np.zeros((self.ndim, self.nyears))
        self.inundated_original_asset_fraction         = np.zeros((self.ndim, self.nyears))
        self.retreated_original_asset_fraction         = np.zeros((self.ndim, self.nyears))

        # Asset retreat
        self.annual_reactive_asset_retreat  = np.zeros((self.ndim, self.nyears))
        self.annual_proactive_asset_retreat = np.zeros((self.ndim, self.nyears))
        self.annual_total_asset_retreat     = np.zeros((self.ndim, self.nyears))
        self.asset_demolition_cost          = np.zeros((self.ndim, self.nyears))
        self.asset_relocation_cost          = np.zeros((self.ndim, self.nyears))
        self.assets_lost_during_retreat     = np.zeros((self.ndim, self.nyears))
        self.total_annual_cost_of_retreat   = np.zeros((self.ndim, self.nyears))

        # Loss of area
        self.inundated_area = np.zeros((self.ndim, self.nyears))
        self.abandoned_area = np.zeros((self.ndim, self.nyears))
        self.effective_retreat_height = np.zeros((self.ndim, self.nyears))

        # net flood height
        # This is the net of global SLR and the building of flood protection
        self.effective_flood_height        = np.zeros((self.ndim, self.nyears))

        
        ################################################################################
        #### Initial values and fit parameters for pre-defined model versions ##########
        ################################################################################
        
        if version in ['DIVA_global', 'DIVA_bipolar', 'DIVA_regional']:
            # These are predefined aggregation levels that load input parameters from files
            # Input parameters are general aggregated information from the DIVA dataset 
            # and the fit parameters used to match DIVA and CIAM data of inundation, surge exposure etc.
            # -> We include all segments from DIVA, also those with a population of 0, as there will still be a loss of land.

            ### Loading the general aggregated information
            path = input_path+'aggregated_data/'
            df = pd.read_csv(path+version+'_information.csv', delimiter=',')
            
            # The factor multiplied with SLR_total to get regional SLR. Only used in case the individual components are not given
            # so that regional SLR has to be fitted linearly to total global SLR
            self.SLR_factor              = df.SLR_factor.values[:]
            # The weights of the individual SLR components for each region to get regional SLR (in case this is activated).
            # Not included are thermosteric and LWS SLR, because weights are 1 in BRICK.
            self.SLR_weight_MG           = df.SLR_weight_MG.values[:]   
            self.SLR_weight_GIS          = df.SLR_weight_GIS.values[:]  
            self.SLR_weight_AIS          = df.SLR_weight_AIS.values[:]
            
            # Information on the coastal regions (aggregated segment data from DIVA)
            self.average_fp_height_init  = df.average_fp_height.values[:]
            self.total_fp_length         = df.total_fp_length.values[:]
            self.coastal_population_init = df.population_fraction.values[:] * 575.6 # The total population in the segments from DIVA
            self.coastal_assets_init     = df.asset_fraction.values[:] * 29079.6 # using segment pop, GDPpc and an asset-GDP ratio of 3
            

            # 2 cases (parameters with and without initial flood protection), 4 parameters and X regions
            self.inund_params_area = np.zeros((2,4,len(self.SLR_factor)))
            self.inund_params_assets = np.zeros((2,4,len(self.SLR_factor)))
            self.inund_params_people = np.zeros((2,4,len(self.SLR_factor)))
            self.storm_suscept_params_assets = np.zeros((2,4,len(self.SLR_factor)))
            self.storm_suscept_params_people = np.zeros((2,4,len(self.SLR_factor)))
            self.storm_exposure_params_assets = np.zeros((2,4,len(self.SLR_factor)))
            self.storm_exposure_params_people = np.zeros((2,4,len(self.SLR_factor)))

            # This determines the parameters of the fit to the DIVA dataset. 
            # In CIAM, the "No adaptation case" apparently assumes that there is no initial seawall
            # at the start of the integration. If "include_initial_fp=True" then the fit
            # assumes that the initial seawall has the height given in the DIVA dataset.
            # --> Fitted only to positive values of effective flood height, becaus negative flood height
            # (i.e. better protection than in reference case) is not interesting here.
            
            parameters = [self.inund_params_area, self.inund_params_assets, self.inund_params_people, self.storm_suscept_params_assets, 
                          self.storm_suscept_params_people, self.storm_exposure_params_assets, self.storm_exposure_params_people]
            filenames_extensions = ['inund_params_area', 'inund_params_assets', 'inund_params_people', 'storm_suscept_params_assets',
                                    'storm_suscept_params_people', 'storm_exposure_params_assets', 'storm_exposure_params_people']

            path = input_path+'fit_function_parameters/'
            for iparam, parameter in enumerate(parameters):
                # Fits including the initial flood protection
                # Here, the data were fitted using a logistic function: y = A / (1 + exp(-k*(x-x0)))
                df = pd.read_csv(path+version+'_'+filenames_extensions[iparam]+'.csv', delimiter=',')
                parameter[0,0,:] = df.k.values[:]
                parameter[0,1,:] = df.x0.values[:]
                parameter[0,2,:] = df.A.values[:]
                parameter[0,3,:] = df.c.values[:]
                
                # Fits without initial flood protection
                # -> only need 3 parameters for the fit, leave 4th at 0
                # In this case, fitted using a log function: y = max(0, a*ln(b*x+1)) + c
                df = pd.read_csv(path+version+'_'+filenames_extensions[iparam]+'_no_initial_dikes.csv', delimiter=',')
                parameter[1,0,:] = df.a.values[:]
                parameter[1,1,:] = df.b.values[:]
                parameter[1,2,:] = df.c.values[:]
        else:
            sys.exit('Specified model version "'+version+'" not allowed! Available options are "DIVA_global", "DIVA_bipolar", "DIVA_regional".')

        ################################################################################
        ################################################################################
        ################################################################################


    def fit_function(self, x, params, ifunc):
        if ifunc==0:
            k = params[ifunc,0]
            x0 = params[ifunc,1]
            A = params[ifunc,2]
            c = params[ifunc,3]
            return A/(1.0 + np.exp(-k*(x-x0))) + c
        elif ifunc==1:
            a=params[ifunc,0]
            b=params[ifunc,1]
            c=params[ifunc,2]
            return a*np.log(b*x+1)+c

    def inverse_fit_function(self, F, params, ifunc):
        if ifunc==0:
            k = params[ifunc,0]
            x0 = params[ifunc,1]
            A = params[ifunc,2]
            c = params[ifunc,3]
            return (-1.0/k) * np.log(A / (F - c) - 1.0) + x0
        elif ifunc==1:
            a=params[ifunc,0]
            b=params[ifunc,1]
            c=params[ifunc,2]
            return (np.exp((F - c) / a) - 1.0) / b


    def __update_single_param(self, prange, scale_factor):
        return prange[0] + scale_factor * (prange[1] - prange[0])

    def __randomize_parameters(self):
        self.flood_event_fatality_rate = self. __update_single_param(self.flood_event_fatality_rate_range, np.random.rand())
        self.fraction_of_storm_damages_that_is_repaired = self. __update_single_param(self.fraction_of_storm_damages_that_is_repaired_range, np.random.rand())
        self.flood_event_damage_fraction = self. __update_single_param(self.flood_event_damage_fraction_range, np.random.rand())
        self.effective_flood_height_at_which_investment_is_halved = self. __update_single_param(self.effective_flood_height_at_which_investment_is_halved_range, np.random.rand())
        self.safe_coastal_zone_likelihood_threshold = self. __update_single_param(self.safe_coastal_zone_likelihood_threshold_range, np.random.rand())
        self.fraction_of_investments_that_must_be_at_the_coast = self. __update_single_param(self.fraction_of_investments_that_must_be_at_the_coast_range, np.random.rand())
        self.maximum_gdp_fraction_for_fp_investment = self. __update_single_param(self.maximum_gdp_fraction_for_fp_investment_range, np.random.rand())
        self.fp_construction_duration = self. __update_single_param(self.fp_construction_duration_range, np.random.rand())
        self.fp_construction_cost_reference = self. __update_single_param(self.fp_construction_cost_reference_range, np.random.rand())
        self.maintenance_cost_fraction = self. __update_single_param(self.maintenance_cost_fraction_range, np.random.rand())
        self.coastal_land_value_init = self. __update_single_param(self.coastal_land_value_init_range, np.random.rand())
        self.land_opportunity_cost_rate = self. __update_single_param(self.land_opportunity_cost_rate_range, np.random.rand())
        self.mobile_asset_fraction = self. __update_single_param(self.mobile_asset_fraction_range, np.random.rand())
        self.asset_relocation_cost_factor = self. __update_single_param(self.asset_relocation_cost_factor_range, np.random.rand())
        self.asset_demolition_cost_factor = self. __update_single_param(self.asset_demolition_cost_factor_range, np.random.rand())
        self.not_depreciated_fraction_of_assets_at_time_of_retreat = self. __update_single_param(self.not_depreciated_fraction_of_assets_at_time_of_retreat_range, np.random.rand())
        self.increase_factor_for_costs_of_reactive_retreat = self. __update_single_param(self.increase_factor_for_costs_of_reactive_retreat_range, np.random.rand())
        self.proactive_retreat_time_scale = self. __update_single_param(self.proactive_retreat_time_scale_range, np.random.rand())
        return

    def getUncertaintyParameters(self):
        uncertaintyParameters = [
                np.copy(self.flood_event_fatality_rate),
                np.copy(self.fraction_of_storm_damages_that_is_repaired),
                np.copy(self.flood_event_damage_fraction),
                np.copy(self.effective_flood_height_at_which_investment_is_halved),
                np.copy(self.safe_coastal_zone_likelihood_threshold),
                np.copy(self.fraction_of_investments_that_must_be_at_the_coast),
                np.copy(self.maximum_gdp_fraction_for_fp_investment),
                np.copy(self.fp_construction_duration),
                np.copy(self.fp_construction_cost_reference),
                np.copy(self.maintenance_cost_fraction),
                np.copy(self.coastal_land_value_init),
                np.copy(self.land_opportunity_cost_rate),
                np.copy(self.mobile_asset_fraction),
                np.copy(self.asset_relocation_cost_factor),
                np.copy(self.asset_demolition_cost_factor),
                np.copy(self.not_depreciated_fraction_of_assets_at_time_of_retreat),
                np.copy(self.increase_factor_for_costs_of_reactive_retreat),
                np.copy(self.proactive_retreat_time_scale)
            ]
        return uncertaintyParameters

    def setUncertaintyParameters(self, InputParameters):

        self.flood_event_fatality_rate                              = InputParameters[0]
        self.fraction_of_storm_damages_that_is_repaired             = InputParameters[1]
        self.flood_event_damage_fraction                            = InputParameters[2]
        self.effective_flood_height_at_which_investment_is_halved   = InputParameters[3]
        self.safe_coastal_zone_likelihood_threshold                 = InputParameters[4]
        self.fraction_of_investments_that_must_be_at_the_coast      = InputParameters[5]
        self.maximum_gdp_fraction_for_fp_investment                 = InputParameters[6]
        self.fp_construction_duration                               = InputParameters[7]
        self.fp_construction_cost_reference                         = InputParameters[8]
        self.maintenance_cost_fraction                              = InputParameters[9]
        self.coastal_land_value_init                                = InputParameters[10]
        self.land_opportunity_cost_rate                             = InputParameters[11]
        self.mobile_asset_fraction                                  = InputParameters[12]
        self.asset_relocation_cost_factor                           = InputParameters[13]
        self.asset_demolition_cost_factor                           = InputParameters[14]
        self.not_depreciated_fraction_of_assets_at_time_of_retreat  = InputParameters[15]
        self.increase_factor_for_costs_of_reactive_retreat          = InputParameters[16]
        self.proactive_retreat_time_scale                           = InputParameters[17]

        return




    def __initialise_variables(self):

        if self.dbg == 1: print('Initialising variables!')

        # First, randomize the uncertainty parameters, if this is activated.
        if self.randomize: self.__randomize_parameters()

        ### The (regional/segregated/global but wheighted) SLR is different from the calculated total SLR!
        if self.include_SLR_components:
            self.SLR_thermo = self.SLR_components[0,np.newaxis,:] * np.ones(self.ndim)[:,np.newaxis]
            self.SLR_LWS    = self.SLR_components[1,np.newaxis,:] * np.ones(self.ndim)[:,np.newaxis]
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
        self.ypc_US_2010 *= self.USD_fac
        self.coastal_assets_init *= self.USD_fac
        self.coastal_land_value_init *= self.USD_fac
        self.fp_construction_cost_reference *= self.USD_fac

        # Copy the input references arrays into the variables that we track
        self.coastal_population = np.copy(self.population)
        self.coastal_GDP = np.copy(self.GDP)
        self.coastal_GDPperCapita = self.GDP / self.population

        # Calculating flood damage resilience as in CIAM, depending on GDP per capita
        self.storm_damage_resilience[:,0] = self.coastal_GDPperCapita[:,0] / (self.coastal_GDPperCapita[:,0] + self.ypc_US_2010)

        # Initialise main stocks
        if self.include_initial_fp:
            self.average_fp_height[:,0]  = self.average_fp_height_init[:]
        else:
            self.average_fp_height[:,0]  = 0.0
        self.coastal_assets[:,0]     = self.coastal_assets_init[:]
        self.coastal_population[:,0] = self.coastal_population_init[:]


        # Now calculate the initial fractions of assets and people that are exposed to storm surges generally, as well as the
        # initial values for the inundated asset, people and area fractions. The latter should theoretically be 0 with the original data,
        # but the fitted function produces very small values for a SLR of 0 meters. 
        # Also, if SLR is not 0 in the initial year, then this will also be > 0.
    
        # Calculation uses different fit functions, depending on the above question (logistic and logarithmic functions)
        self.orig_susceptible_asset_fraction[:,0]  = self.fit_function(self.effective_flood_height[:,0], self.storm_suscept_params_assets, self.iFit)
        self.orig_susceptible_people_fraction[:,0] = self.fit_function(self.effective_flood_height[:,0], self.storm_suscept_params_people, self.iFit)
        self.orig_exposed_asset_fraction[:,0]      = self.fit_function(self.effective_flood_height[:,0], self.storm_exposure_params_assets, self.iFit)
        self.orig_exposed_people_fraction[:,0]     = self.fit_function(self.effective_flood_height[:,0], self.storm_exposure_params_people, self.iFit)

        self.inundated_original_asset_fraction[:,0]  = self.fit_function(self.effective_flood_height[:,0], self.inund_params_assets, self.iFit)
        self.inundated_original_people_fraction[:,0] = self.fit_function(self.effective_flood_height[:,0], self.inund_params_people, self.iFit)
        self.inundated_area[:,0]                     = self.fit_function(self.effective_flood_height[:,0], self.inund_params_area, self.iFit)

        self.retreated_original_asset_fraction[:,0]  = self.fit_function(self.effective_flood_height[:,0], self.inund_params_assets, self.iFit)
        self.retreated_original_people_fraction[:,0] = self.fit_function(self.effective_flood_height[:,0], self.inund_params_people, self.iFit)


        # Check if given willingness_to_invest variable is time dependent variable with correct length or a scalar
        if not isinstance(self.willingness_to_invest_in_fp, (list, np.ndarray)):
            self.willingness_to_invest_in_fp = np.zeros((self.ndim, self.nyears))+self.willingness_to_invest_in_fp
        elif len(self.willingness_to_invest_in_fp[0,:]) != self.nyears or  len(self.willingness_to_invest_in_fp[:,0]) != self.ndim:
            sys.exit('Length of willingness to invest array does not fit with time or aggregation dimension')
        # Check if given willingness_to_retreat variable is time dependent variable with correct length or a scalar
        if not isinstance(self.willingness_to_retreat, (list, np.ndarray)):
            self.willingness_to_retreat = np.zeros((self.ndim, self.nyears))+self.willingness_to_retreat
        elif len(self.willingness_to_retreat[0,:]) != self.nyears or  len(self.willingness_to_retreat[:,0]) != self.ndim:
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
        self.annual_costs_of_fp_maintenance_noadapt[:,i] = self.average_fp_height[:,0] * self.total_fp_length \
                                                   * self.construction_cost[:,i] *  self.maintenance_cost_fraction



        ### HOW MUCH INVESTMENT?
        # Current net flood height + expected SLR is missing protection
        missing_protection = self.effective_flood_height[:,i] + self.expected_SLR_in_50_years[:,i]
        desired_protection = missing_protection * self.willingness_to_invest_in_fp[:,i]

        cost_of_reaching_desired_protection = np.maximum(0, self.construction_cost[:,i] * self.total_fp_length * \
                    (( desired_protection + self.average_fp_height[:,i])**2 - self.average_fp_height[:,i]**2) )
    
        maximum_money_available_for_fp = np.maximum(0, self.coastal_GDP[:,i] * self.maximum_gdp_fraction_for_fp_investment \
                                            - self.annual_costs_of_fp_maintenance[:,i] )

            
        if self.include_fp_investment_cap:
            self.effective_annual_investment_in_fp[:,i] = np.minimum(cost_of_reaching_desired_protection / self.fp_construction_duration,
                                                                       maximum_money_available_for_fp)
        else:
            self.effective_annual_investment_in_fp[:,i] = cost_of_reaching_desired_protection / self.fp_construction_duration

        # Calculate potential maximum fp height increase for asset investment decision making
        potential_max_investment_over_50_years = self.effective_annual_investment_in_fp[:,i] * 50.0
        self.potential_fp_height_increase_over_50_years[:,i] =  np.sqrt( self.average_fp_height[:,i]**2 \
                        + potential_max_investment_over_50_years / (self.total_fp_length * self.construction_cost[:,i]) ) - self.average_fp_height[:,i]



        ### CHANGE IN FLOOD PROTECTION HEIGHT
        self.annual_change_in_fp_height[:,i] = np.sqrt( self.average_fp_height[:,i]**2 + self.effective_annual_investment_in_fp[:,i] \
                    / (self.total_fp_length * self.construction_cost[:,i]) ) - self.average_fp_height[:,i]

        if i >= self.nyears-1: return
        ##########################################################################################
        ###          
        ### UPDATE THE FLOOD PROTECTION HEIGHT
        ###
        ##########################################################################################
        self.average_fp_height[:,i+1] = self.average_fp_height[:,i] + self.annual_change_in_fp_height[:,i]            
        
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
            #  1. Reactive retreat in response to inundation, which has very high costs.
            #  2. Proactive retreat in response to expected inundation in 50 years, which has much lower costs.
        
            # Reactive retreat (and hence flood damage) only after first timestep, because it is in response to SLR from one time step to the next
            if i==0: 
                self.annual_reactive_people_retreat[:,i] = 0.0
                total_removed_fraction = self.inundated_original_people_fraction[:,0] # use the initialisation value in first timestep         
            else:
                # First calculate the previous fraction of original people distribution that was already removed (inundated or retreated)
                total_removed_fraction = np.maximum(self.inundated_original_people_fraction[:,i-1], self.retreated_original_people_fraction[:,i-1])
                self.annual_reactive_people_retreat[:,i] = self.coastal_population[:,i] \
                                                            * np.maximum(0.0, (self.inundated_original_people_fraction[:,i]\
                                                              - total_removed_fraction) / (1.0 - total_removed_fraction) )
                # update removed fraction for calculation of proactive retreat
                total_removed_fraction = np.maximum(self.inundated_original_people_fraction[:,i], self.retreated_original_people_fraction[:,i-1])


            # what is the flood height in 50 years?
            expected_effective_flood_height = self.effective_flood_height[:,i] + self.expected_SLR_in_50_years[:,i] \
                                            - self.potential_fp_height_increase_over_50_years[:,i]

            expected_susceptible_fraction = self.fit_function(expected_effective_flood_height, self.storm_suscept_params_people, self.iFit)

            retreating_people_fraction = (self.willingness_to_retreat[:,i] / self.proactive_retreat_time_scale) \
                    *  np.maximum(0, expected_susceptible_fraction - total_removed_fraction)

            if i==0: self.retreated_original_people_fraction[:,i] = self.retreated_original_people_fraction[:,0] + retreating_people_fraction
            else: self.retreated_original_people_fraction[:,i] = self.retreated_original_people_fraction[:,i-1] + retreating_people_fraction
    
            self.annual_proactive_people_retreat[:,i] = self.coastal_population[:,i] * retreating_people_fraction / (1.0 - total_removed_fraction)
    
            # Total asset retreat is sum of reactive and proactive retreat
            self.annual_total_people_retreat[:,i] = self.annual_reactive_people_retreat[:,i] +  self.annual_proactive_people_retreat[:,i]

            ##########################################################################################
            ####
            #### Number of people flooded in storms and flood fatalities
            ####
            ########################################################################################## 

            # Storm damage to people is depending on what is the actually susceptible fraction of people. 
            # This has to account for the fraction of people that is generally susceptible and the fraction
            # of people that has already been removed from the coast (via inundation or retreat)
            total_removed_fraction = np.maximum(self.inundated_original_people_fraction[:,i], self.retreated_original_people_fraction[:,i])
            actually_susceptible_people_fraction = np.maximum(0, (self.orig_susceptible_people_fraction[:,i] - total_removed_fraction) \
                                                     / (1.0 - total_removed_fraction))

            if self.include_retreat_exposure_reduction:
                people_exposure_reduction_because_of_retreat = actually_susceptible_people_fraction/self.orig_susceptible_people_fraction[:,i]
            else:
                people_exposure_reduction_because_of_retreat = 1.0

            # Subtract here the initial exposure fraction from current exposure fraction to get the SLR driven number
            self.annual_people_flooded[:,i] = self.coastal_population[:,i] * people_exposure_reduction_because_of_retreat \
                                              * np.maximum(0.0, self.orig_exposed_people_fraction[:,i] - self.orig_exposed_people_fraction[:,0])

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

            # There are two ways how asset retreat from the coast:
            #  1. Reactive retreat in response to inundation, which has very high costs.
            #  2. Proactive retreat in response to expected inundation in 50 years, which has much lower costs.

            # Reactive retreat (and hence flood damage) only after first timestep, because it is in response to SLR from one time step to the next
            if i==0: 
                self.annual_reactive_asset_retreat[:,i] = 0.0
                total_removed_fraction = self.inundated_original_asset_fraction[:,0] # use the initialisation value in first timestep
            else:
                # First calculate the previous fraction of original asset distribution that was already removed (inundated or retreated)
                total_removed_fraction = np.maximum(self.inundated_original_asset_fraction[:,i-1], self.retreated_original_asset_fraction[:,i-1])
                self.annual_reactive_asset_retreat[:,i] = self.coastal_assets[:,i] * np.maximum(0.0, (self.inundated_original_asset_fraction[:,i]\
                                                                    - total_removed_fraction) / (1.0 - total_removed_fraction) )
                # update removed fraction for calculation of proactive retreat
                total_removed_fraction = np.maximum(self.inundated_original_asset_fraction[:,i], self.retreated_original_asset_fraction[:,i-1])


            # what is the flood height in 50 years?
            expected_effective_flood_height = self.effective_flood_height[:,i] + self.expected_SLR_in_50_years[:,i] \
                                            - self.potential_fp_height_increase_over_50_years[:,i]

            expected_susceptible_fraction = self.fit_function(expected_effective_flood_height, self.storm_suscept_params_assets, self.iFit)

            retreating_asset_fraction = (self.willingness_to_retreat[:,i] / self.proactive_retreat_time_scale) \
                    * np.maximum(0, expected_susceptible_fraction - total_removed_fraction)
                                            


            if i==0: self.retreated_original_asset_fraction[:,i] = self.retreated_original_asset_fraction[:,0] + retreating_asset_fraction
            else: self.retreated_original_asset_fraction[:,i] = self.retreated_original_asset_fraction[:,i-1] + retreating_asset_fraction
    
            self.annual_proactive_asset_retreat[:,i] = self.coastal_assets[:,i] * retreating_asset_fraction / (1.0 - total_removed_fraction)
    
            # Total asset retreat is sum of reactive and proactive retreat
            self.annual_total_asset_retreat[:,i] = self.annual_reactive_asset_retreat[:,i] +  self.annual_proactive_asset_retreat[:,i]



           ##########################################################################################
            ###
            ### STORM DAMAGE to assets
            ###
            ##########################################################################################
        
            # Storm damage to assets is depending on what is the actually susceptible fraction of assets. This has to account for the fraction of
            # assets that is generally susceptible and the fraction of assets that has already been removed from the coast
            # (via inundation or retreat)
            total_removed_fraction = np.maximum(self.inundated_original_asset_fraction[:,i], self.retreated_original_asset_fraction[:,i])
            actually_susceptible_asset_fraction = np.maximum(0, (self.orig_susceptible_asset_fraction[:,i] - total_removed_fraction)\
                                                                  / (1.0 - total_removed_fraction))
            if self.include_retreat_exposure_reduction:
                asset_exposure_reduction_because_of_retreat = actually_susceptible_asset_fraction/self.orig_susceptible_asset_fraction[:,i]
            else:
                asset_exposure_reduction_because_of_retreat = 1.0
            
            # Subtract here the initial exposure fraction from current exposure fraction to get the SLR driven number
            self.annual_storm_damage_to_assets[:,i] = self.coastal_assets[:,i] * asset_exposure_reduction_because_of_retreat \
                                              * self.flood_event_damage_fraction * (1.0 - self.storm_damage_resilience[:,i]) \
                                              * np.maximum(0.0, self.orig_exposed_asset_fraction[:,i] - self.orig_exposed_asset_fraction[:,0])


        ##########################################################################################
        ###
        ### GROWTH of assets as GDP grows.
        ###
        ##########################################################################################


        # TODO: Account for urbanization

        # Calculate likelihood of investment in coastal zone
        # -> Less growth if it is expected that net flood height will increase (or is already high)
        if self.include_reduced_growth and self.damage:
            dflood =  np.maximum(0, self.expected_SLR_in_50_years[:,i] + self.effective_flood_height[:,i]\
                                     - self.potential_fp_height_increase_over_50_years[:,i])
            self.likelihood_of_investment_in_coastal_zones[:,i] = (1.0 - np.maximum(0, dflood / (dflood \
                                            + self.effective_flood_height_at_which_investment_is_halved ))) * expected_susceptible_fraction \
                                            + (1.0 - expected_susceptible_fraction)
        else:
            self.likelihood_of_investment_in_coastal_zones[:,i] = 1.0

        # No need to update next years asset value, if we are already in the last time step...
        if i >= self.nyears-1: return      


        # Assets growth is being redistributed between coastal zones and inland, in case that there are expected damages
        # First calculate theoretical asset growth in each coastal zone
        theoretical_asset_growth = self.coastal_assets[:,i] * (self.GDP[:,i+1] / self.GDP[:,i] - 1.0)

        # initialising some arrays for asset growth redistribution
        actual_asset_growth = np.copy(theoretical_asset_growth)
        safe_zone = np.ones_like(theoretical_asset_growth)
        growth_moved_to_safe_zones = 0.0

        # Remove small part of the growth from the insufficiently protected assets
        for icz, likelihood in enumerate(self.likelihood_of_investment_in_coastal_zones[:,i]):
            if likelihood < self.safe_coastal_zone_likelihood_threshold:
                actual_asset_growth[icz] = theoretical_asset_growth[icz] * likelihood
                growth_moved_to_safe_zones += theoretical_asset_growth[icz] * (1.0 - likelihood)
                safe_zone[icz] = 0
        

        # Get the fraction of each coastal zone with likelihood > threshold to all coastal zones with likelihood > threshold
        # avoid div by 0, in case there is no safe coastal zone
        if np.sum(safe_zone) != 0.0:
            safe_asset_fractions = (self.coastal_assets[:,i]*safe_zone) / np.sum(self.coastal_assets[:,i]*safe_zone)
        else:
            safe_asset_fractions = np.zeros_like(safe_zone)
        asset_fractions = (self.coastal_assets[:,i]) / np.sum(self.coastal_assets[:,i])
        
        # Add some of the removed growth to well protected assets,
        # but not all of it! Some is moved away from the coast, because it is not coast-specific.
        growth_staying_at_coast = growth_moved_to_safe_zones * self.fraction_of_investments_that_must_be_at_the_coast
        growth_moving_away_from_coast = growth_moved_to_safe_zones - growth_staying_at_coast

        for icz, likelihood in enumerate(self.likelihood_of_investment_in_coastal_zones[:,i]):
            # if there is no safe coastal zone option, then the investment that has to stay in the coastal zone is 
            # made in the original coastal zone without moving it around.
            if sum(safe_zone) == 0 or not self.move_around_growth:
                actual_asset_growth[icz] += growth_staying_at_coast * asset_fractions[icz]
            # distribute the removed investment over the safe coastal zones relative to their fraction of all protected assets
            elif likelihood >= self.safe_coastal_zone_likelihood_threshold:
                actual_asset_growth[icz] += growth_staying_at_coast * safe_asset_fractions[icz]

        ##########################################################################################
        ###
        ### UPDATE Coastal asset stock
        ###
        ##########################################################################################
        # The asset feedback switch (default: 0) determines what fraction of storm damages actually reduces coastal asset values.
        # In CIAM this is 0, because it is assumed that storm surge damages are repaired always.
        self.coastal_assets[:,i+1] = self.coastal_assets[:,i] + actual_asset_growth[:] - self.annual_total_asset_retreat[:,i] \
                                   - self.asset_feedback_switch \
                                   * (1.0 - self.fraction_of_storm_damages_that_is_repaired) * self.annual_storm_damage_to_assets[:,i]
        return
        


    def __prepare_next_timestep(self,i):
        
        # Calculate GDP(perCapita) internally or leave as external input?
        if self.include_gdp_effect:
            self.coastal_GDP[:,i+1] = self.coastal_GDPperCapita[:,0] * self.coastal_population[:,0] \
                                       * self.coastal_assets[:,i+1] / self.coastal_assets[:,0]
            self.coastal_GDPperCapita[:,i+1] = self.coastal_GDP[:,i+1] / self.coastal_population[:,i+1]
           
        # Calculating flood damage resilience as in CIAM, depending on GDP per capita
        self.storm_damage_resilience[:,i+1] = self.coastal_GDPperCapita[:,i+1] / (self.coastal_GDPperCapita[:,i+1] + self.ypc_US_2010)

        
        #### Updating the value of land using a formula from CIAM repository
        self.landvalue_appreciation_factor[:,i+1] = self.landvalue_appreciation_factor[:,i] \
                    * np.exp(0.565 * (self.coastal_GDPperCapita[:,i+1]/self.coastal_GDPperCapita[:,i] - 1.0) \
                            + 0.313 * (self.coastal_population[:,i+1]/self.coastal_population[:,i] - 1.0) )

        # The net global flood height 
        self.effective_flood_height[:,i+1] = (self.SLR[:,i+1]-self.SLR[:,0]) - (self.average_fp_height[:,i+1] - self.average_fp_height[:,0])
        
        # Update the susceptible and inundated theoretical fractions for the next time step, based on the updated net flood height
        dflood = self.effective_flood_height

        self.orig_susceptible_asset_fraction[:,i+1]  = self.fit_function(dflood[:,i+1], self.storm_suscept_params_assets, self.iFit)
        self.orig_susceptible_people_fraction[:,i+1] = self.fit_function(dflood[:,i+1], self.storm_suscept_params_people, self.iFit)
        self.orig_exposed_asset_fraction[:,i+1]      = self.fit_function(dflood[:,i+1], self.storm_exposure_params_assets, self.iFit)
        self.orig_exposed_people_fraction[:,i+1]     = self.fit_function(dflood[:,i+1], self.storm_exposure_params_people, self.iFit)

        self.inundated_original_asset_fraction[:,i+1] = np.maximum(self.inundated_original_asset_fraction[:,i], 
                                                                   self.fit_function(dflood[:,i+1], self.inund_params_assets, self.iFit))
        self.inundated_original_people_fraction[:,i+1] = np.maximum(self.inundated_original_people_fraction[:,i],
                                                                    self.fit_function(dflood[:,i+1], self.inund_params_people, self.iFit))
        self.inundated_area[:,i+1] = np.maximum(self.inundated_area[:,i], self.fit_function(dflood[:,i+1], self.inund_params_area, self.iFit))

        # To get the abandoned area during retreat, we first calculate the effective retreat height using the asset inundation without
        # initial flood protection (the inverse of the log function used to calculate inundated assets in this case)
        self.effective_retreat_height[:,i] = self.inverse_fit_function(self.retreated_original_asset_fraction[:,i], self.inund_params_assets, 1)
        self.abandoned_area[:,i+1] = np.maximum(self.abandoned_area[:,i], 
                                                self.fit_function(self.effective_retreat_height[:,i],self.inund_params_area, 1))


        return

    def __make_outputs(self):
        
        # ASSETS
        self.assets_lost_during_retreat = (1.0 - self.mobile_asset_fraction) * (self.annual_reactive_asset_retreat \
                                               + self.annual_proactive_asset_retreat * self.not_depreciated_fraction_of_assets_at_time_of_retreat)

        self.asset_relocation_cost = self.annual_total_asset_retreat * self.mobile_asset_fraction * self.asset_relocation_cost_factor
        self.asset_demolition_cost = self.annual_total_asset_retreat * self.asset_demolition_cost_factor * (1.0 - self.mobile_asset_fraction)

        
        # PEOPLE
        # Reactive retreat is five times as costly as proactive retreat, just like in CIAM
        self.people_retreat_cost = (self.annual_proactive_people_retreat + self.increase_factor_for_costs_of_reactive_retreat * self.annual_reactive_people_retreat) * self.coastal_GDPperCapita

        # CIAM reference: This is for comparison to CIAM output

        # LAND
        # Opportunity costs
        # The opportunity cost is assumed to be 4% of the landvalue annually
        # First calculate the value of land per square kilometer
        land_cost_per_sqkm = self.landvalue_appreciation_factor * self.coastal_land_value_init * self.land_opportunity_cost_rate

        # 1. Opportunity cost: value of land used up for flood protection
        #   - DON'T (subtract the initial fp height, in order to only get the additional opportunity costs), not done in CIAM
        #   - Factor 1.7 comes from ratio between dike height and width (60° angle, as in CIAM)
        #   - In the end divided by 2, as in CIAM
        self.fp_land_opportunity_cost = self.total_fp_length[:,np.newaxis] * 1.7 * self.average_fp_height * 0.001 * land_cost_per_sqkm / 2.0

        # 2. Opportunity cost: value of land that is inundated or abandoned
        lost_area = np.maximum(self.inundated_area, self.abandoned_area)
        self.lost_land_opportunity_cost =  lost_area * land_cost_per_sqkm

        # Total opportunity costs of the land that is lost
        self.total_opportunity_cost = self.fp_land_opportunity_cost + self.lost_land_opportunity_cost


        #####
        # CIAM reference: These are for comparison to CIAM output
        #####
        self.relocation_cost = self.asset_relocation_cost + self.asset_demolition_cost + self.people_retreat_cost
        self.flood_cost = self.assets_lost_during_retreat + self.lost_land_opportunity_cost
        self.construct_cost = self.fp_land_opportunity_cost + self.effective_annual_investment_in_fp + self.annual_costs_of_fp_maintenance
        self.net_construct_cost = self.fp_land_opportunity_cost + self.effective_annual_investment_in_fp \
                + self.annual_costs_of_fp_maintenance - self.annual_costs_of_fp_maintenance_noadapt
        self.storm_cost = self.annual_storm_damage_to_assets


        return

    def getTotalCoastalGDP(self): return np.sum(self.coastal_GDP,axis=0)
    def getTotalGDPperCapita(self): return np.sum(self.coastal_GDP,axis=0) / np.sum(self.coastal_population,axis=0)
    def getTotalCoastalAssets(self): return np.sum(self.coastal_assets, axis=0)
    def getTotalCoastalPopulation(self): return np.sum(self.coastal_population, axis=0)

    def getTotalStormCost(self): return np.sum(self.storm_cost, axis=0)
    def getTotalRelocationCost(self): return np.sum(self.relocation_cost, axis=0)
    def getTotalFloodCost(self): return np.sum(self.flood_cost, axis=0)
    def getTotalConstructCost(self): return np.sum(self.construct_cost, axis=0)
    def getTotalNetConstructCost(self): return np.sum(self.net_construct_cost, axis=0)


    def getTotalProtectionInvestment(self): return np.sum(self.effective_annual_investment_in_fp, axis=0)
    def getTotalMaintenanceCost(self): return np.sum(self.annual_costs_of_fp_maintenance, axis=0)
    def getNetMaintenanceCost(self): return self.annual_costs_of_fp_maintenance - self.annual_costs_of_fp_maintenance_noadapt
    def getTotalNetMaintenanceCost(self): return np.sum(self.annual_costs_of_fp_maintenance - self.annual_costs_of_fp_maintenance_noadapt, axis=0)

    def getTotalOpportunityCost(self): return np.sum(self.total_opportunity_cost, axis=0)

    def getTotalPeopleFlooded(self): return np.sum(self.annual_people_flooded, axis=0)
    def getTotalFloodFatalities(self): return np.sum(self.annual_flood_fatalities, axis=0)
    





        





