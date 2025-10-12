# © Lennart Ramme, MPI-M, 2024

import numpy as np
import sys

class globalSLRModel:
    '''This is a model to calculate global SLR based on the contributions from five components.
    Required input to this model are only the global mean surface air temperature anomaly and ocean heat content changes.
    Global population is an optional input parameter.
    '''

    def __init__(self, T, OHC_change, population=[-999], sy=1750, ey=2300, dt=1.0, dbg=0, include_highImpact=False, randomize=False):

        ########## Set some global variables ################
        self.T_anomaly = T
        self.OHC_change = OHC_change
        self.include_highImpact = include_highImpact
        self.randomize = randomize
        
        if population[0] == -999:
            self.with_population = False
            self.population = np.nan
        else:
            self.with_population = True
            self.population = population
        
        self.dbg = dbg
        
        # Timestep:
        self.sy = sy
        self.ey = ey
        self.dt = dt
        self.nyears = int((self.ey - self.sy) / self.dt) + 1
        self.time = np.linspace(self.sy, self.ey, self.nyears)

        ########## Parameter settings ###################

        ####################################################################################
        ###########          Calibrated Parameters          ################################
        ####################################################################################

        ### These are the scale factors that can be varied within [0,1] to cover the full uncertainty range
        self.thermo_scale_factor = 0.5
        self.LWS_scale_factor = 0.5
        self.MG_scale_factor = 0.5
        self.GIS_scale_factor = 0.5
        self.AIS_scale_factor = 0.5

        ### These are the ranges within which the uncertainty parameters will be scaled
        ### Default value is the mean of the two.

        ### Tuning for thermosteric component
        self.thermo_m_per_YJ_range = (0.1, 0.12)

        ### Tuning for Land water storage component
        self.LWS_rate_range = (0.0002, 0.0004)   # m per yr ; constant rate in case no population data are given
        self.LWS_pop_sens_range = (1.0e-8, 6.0e-8) # m per Mp per yr

        ### Tuning for Mountain glacier component
        # Varying both the sensitivity and elasticity of the T-dep. of SLR_MG
        self.MG_beta_0_range = (0.0004, 0.001)

        ### Tuning for Greenland Ice Sheet component
        # Vary all parameters that are differing between the low and high discharge versions of MAGICC
        self.GIS_rho_range = (1.0e-4, 5.0e-4)
        self.GIS_v_range = (0.5e-4, 2.0e-4)
        # implement high impact behaviour also for GrIS (like for AntIS, but with GSAT anomaly not T_Ant as reference)
        self.GIS_disintegration_rate_range = (0.001, 0.01)
        self.GIS_T_crit_range = (4,2.5) # Note lower number at 2nd position, as lower T_crit induces more melting

        ### Tuning for Antarctic Ice Sheet component
        self.AIS_a_anto_range = (0.15, 0.3)
        self.AIS_h0_range = (1600, 1900)
        self.AIS_disintegration_rate_range = (0.000, 0.008)
        self.AIS_T_crit_range = (-14, -15) # Note lower number at 2nd position, as lower T_crit induces more melting



        ####################################################################################
        ####################################################################################

        ## Other parameters
        # thermosteric component
        self.thermo_startyear = self.sy

        # Land water storage component
        self.LWS_startyear = 1960 # Assume that global groundwater levels were in balance in 1960

        # Mountain glacier component
        self.MG_startyear = 1850
        self.MG_T_eq = 0.0
        self.MG_V_0 = 0.41
        self.MG_n = 1.646
        self.MG_T_exp = 1.5

        # GrIS component
        self.GIS_startyear = self.sy
        self.GIS_p = 2.0
        self.GIS_SMB_max = 7.36    # m
        self.GIS_eps = 0.39
        self.GIS_outlet_max = 0.42


        # AntIS component
        self.AIS_b_anto = 0.5
        self.AIS_alpha = 0.5
        self.AIS_gamma = 2.
        self.AIS_mu = 8.7
        self.AIS_nu = 0.012
        self.AIS_P0 = 0.35
        self.AIS_kappa = 0.04
        self.AIS_f0 = 1.2
        self.AIS_c = 95.
        self.AIS_b0 = 775.
        self.AIS_slope = 0.0006

        self.AIS_startyear = self.sy+1
        self.AIS_rho_w = 1030.
        self.AIS_rho_i = 917.
        self.AIS_rho_m = 4000.
        self.AIS_Rad0 = 1.864e6
        self.AIS_dell = self.AIS_rho_w / self.AIS_rho_i
        self.AIS_eps1 = self.AIS_rho_i /(self.AIS_rho_m - self.AIS_rho_i)
        self.AIS_eps2 = self.AIS_rho_w /(self.AIS_rho_m - self.AIS_rho_i)

        ########## Initialize arrays ################
        self.SLR_thermo = np.zeros((self.nyears))
        self.SLR_LWS    = np.zeros((self.nyears))
        self.SLR_MG     = np.zeros((self.nyears))
        self.SLR_GIS    = np.zeros((self.nyears))
        self.SLR_AIS    = np.zeros((self.nyears))
        self.SLR_total  = np.zeros((self.nyears))

        self.SLR_GIS_SMB = np.zeros((self.nyears))
        self.SLR_GIS_DIS = np.zeros((self.nyears))
        self.GIS_outlet_vdis = np.zeros((self.nyears))


        self.AIS_Volume = np.zeros((self.nyears))
        self.AIS_Radius = np.zeros((self.nyears))


    def reset_SLR(self):
        ########## Initialize arrays again ################
        self.SLR_thermo = np.zeros((self.nyears))
        self.SLR_LWS    = np.zeros((self.nyears))
        self.SLR_MG     = np.zeros((self.nyears))
        self.SLR_GIS    = np.zeros((self.nyears))
        self.SLR_AIS    = np.zeros((self.nyears))
        self.SLR_total  = np.zeros((self.nyears))

        self.SLR_GIS_SMB = np.zeros((self.nyears))
        self.SLR_GIS_DIS = np.zeros((self.nyears))
        self.GIS_outlet_vdis = np.zeros((self.nyears))

        self.AIS_Volume = np.zeros((self.nyears))
        self.AIS_Radius = np.zeros((self.nyears))


    def align(self, year, silent=False):
        index = int(year - self.sy)
        if index < 0 or index >= self.nyears: sys.exit('Cannot align to year outside of simulated range!')
        if not silent: print('Aligning SLR to be 0 in year '+str(year)+' (index: '+str(index)+')')
        self.SLR_thermo = self.SLR_thermo - self.SLR_thermo[index]
        self.SLR_LWS = self.SLR_LWS - self.SLR_LWS[index]
        self.SLR_MG = self.SLR_MG - self.SLR_MG[index]
        self.SLR_GIS_SMB = self.SLR_GIS_SMB - self.SLR_GIS_SMB[index]
        self.SLR_GIS_DIS = self.SLR_GIS_DIS - self.SLR_GIS_DIS[index]
        self.SLR_GIS = self.SLR_GIS - self.SLR_GIS[index]
        self.SLR_AIS = self.SLR_AIS - self.SLR_AIS[index]
        self.SLR_total = self.SLR_total - self.SLR_total[index]


    def integrate(self, silent=True):
        if not silent: print('Start integrating...')

        self.__init_params()

        self.__AIS_init()
        self.__GIS_init()
        for i in range(0, self.nyears-1):
            if self.dbg==1: print('Year: ', i)
            if self.time[i] >= self.thermo_startyear: self.__update_SLR_thermo(i)
            if self.time[i] >= self.LWS_startyear:    self.__update_SLR_LWS(i)
            if self.time[i] >= self.MG_startyear:     self.__update_SLR_MG(i)
            if self.time[i] >= self.GIS_startyear:    self.__update_SLR_GIS(i)
            if self.time[i] >= self.AIS_startyear:    self.__update_SLR_AIS(i)

            self.SLR_total[i+1] = self.SLR_thermo[i+1] + self.SLR_LWS[i+1] + self.SLR_MG[i+1] + self.SLR_GIS[i+1] + self.SLR_AIS[i+1]

        if not silent: print('   ...finished')

    def __init_params(self):

        if self.randomize:
            # Randomize all scale factors to arbitrarily choose values from the uncertainty range for each contribution individually
            all_factors = [self.thermo_scale_factor, self.LWS_scale_factor, self.MG_scale_factor, self.GIS_scale_factor, self.AIS_scale_factor]
            for factor in all_factors: factor = np.random.rand()


        ### Adjusting the parameters using the scale factor
        # Thermo
        self.thermo_m_per_YJ = self.__scale_single_param(self.thermo_m_per_YJ_range, self.thermo_scale_factor)

        # LWS
        self.LWS_rate = self.__scale_single_param(self.LWS_rate_range, self.LWS_scale_factor)
        self.LWS_pop_sens = self.__scale_single_param(self.LWS_pop_sens_range, self.LWS_scale_factor)

        # MG
        self.MG_beta_0 = self.__scale_single_param(self.MG_beta_0_range, self.MG_scale_factor)

        # GIS
        self.GIS_rho = self.__scale_single_param(self.GIS_rho_range, self.GIS_scale_factor)
        self.GIS_v = self.__scale_single_param(self.GIS_v_range, self.GIS_scale_factor)
        self.GIS_T_crit = self.__scale_single_param(self.GIS_T_crit_range, self.GIS_scale_factor)
        self.GIS_disintegration_rate = self.__scale_single_param(self.GIS_disintegration_rate_range, self.GIS_scale_factor)

        # AIS
        self.AIS_a_anto = self.__scale_single_param(self.AIS_a_anto_range, self.AIS_scale_factor)
        self.AIS_h0 = self.__scale_single_param(self.AIS_h0_range, self.AIS_scale_factor)
        self.AIS_T_crit = self.__scale_single_param(self.AIS_T_crit_range, self.AIS_scale_factor)
        self.AIS_disintegration_rate = self.__scale_single_param(self.AIS_disintegration_rate_range, self.AIS_scale_factor)

        if not self.include_highImpact:
            self.AIS_T_crit = 999.
            self.GIS_T_crit = 999.
        return

    def __scale_single_param(self, prange, scale_factor):
        return prange[0] + scale_factor * (prange[1] - prange[0])
        

    def __update_SLR_thermo(self,i):
        if self.dbg==1: print('   SLR thermo: ', i)
        self.SLR_thermo[i+1] = self.SLR_thermo[i] + self.thermo_m_per_YJ * 1.e-24 * self.OHC_change[i]
        return


    def __update_SLR_LWS(self, i):
        if self.dbg==1: print('   SLR LWS: ', i)

        if self.with_population:
            self.SLR_LWS[i+1] = self.SLR_LWS[i] + self.LWS_pop_sens * self.population[i]
        else:
            self.SLR_LWS[i+1] = self.SLR_LWS[i] + self.LWS_rate
        return


    def __update_SLR_MG(self, i):
        if self.dbg==1: print('   SLR MG: ', i)
        change = 0.0
        if self.SLR_MG[i] < self.MG_V_0 and self.T_anomaly[i] > 0.0: 
            change = self.MG_beta_0 * (self.T_anomaly[i] - self.MG_T_eq)**self.MG_T_exp * (1.0 - self.SLR_MG[i]/self.MG_V_0)**self.MG_n
        self.SLR_MG[i+1] = self.SLR_MG[i] + change
        return


    def __GIS_init(self):
        self.GIS_outlet_vdis[:] = self.GIS_outlet_max
        return

    def __update_SLR_GIS(self, i):
        if self.dbg==1: print('   SLR GIS: ', i)
        T = self.T_anomaly[i]

        if T >= 0.0: term1 = T**self.GIS_p
        else: term1 = 0.0
        term2 = np.sqrt(1.0 - self.SLR_GIS_SMB[i]/self.GIS_SMB_max)
    

        if  self.dbg == 'GIS': print(i, 'GIS_v:', self.GIS_v)
        if  self.dbg == 'GIS': print(i, 'GIS_p:', self.GIS_p)

        if self.dbg == 'GIS': print(i, 'Tanomaly:', T)
        if self.dbg == 'GIS': print(i, 'term1:', term1)
        if self.dbg == 'GIS': print(i, 'term2:', term2)

        

        self.SLR_GIS_SMB[i+1] = self.SLR_GIS_SMB[i] + self.GIS_v * term1 * term2
        
        if self.dbg == 'GIS': print(i, 'SLR_GIS_SMB:',  self.SLR_GIS_SMB[i+1])

        # Marine ice cliff implementation
        if T > self.GIS_T_crit:
            disintegration_discharge = np.maximum(0.0, self.GIS_disintegration_rate * (self.GIS_outlet_max - self.SLR_GIS_DIS[i]))
        else:
            disintegration_discharge = 0.0


        if  self.dbg == 'GIS': print(i, 'GIS_eps:', self.GIS_eps)
        if  self.dbg == 'GIS': print(i, 'GIS_rho:', self.GIS_rho)        
        if  self.dbg == 'GIS': print(i, 'GIS R:', disintegration_discharge)        
        if  self.dbg == 'GIS': print(i, 'GIS_outlet_max:', self.GIS_outlet_max)
        if  self.dbg == 'GIS': print(i, 'SLR_GIS_DIS[i]:', self.SLR_GIS_DIS[i])
        

        SID = np.maximum(0, self.GIS_rho * (self.GIS_outlet_max - self.SLR_GIS_DIS[i]) * np.exp(self.GIS_eps*T))
        self.SLR_GIS_DIS[i+1] = self.SLR_GIS_DIS[i] + SID + disintegration_discharge

        self.SLR_GIS[i+1] = self.SLR_GIS_DIS[i+1] + self.SLR_GIS_SMB[i+1]

        if self.dbg == 'GIS':
            print(i, 'SLR_GIS_DIS:', self.SLR_GIS_DIS[i+1])
            print(i, 'SLR_GIS:', self.SLR_GIS[i+1])
            print()
        return


    def __AIS_init(self):
        R = self.AIS_Rad0
        rc = self.AIS_b0/self.AIS_slope
        V = np.pi * (1+self.AIS_eps1) * ( (8./15.) * self.AIS_mu**0.5 * R**2.5 - (1./3.)*self.AIS_slope*R**3)
        if R>rc: V = V - np.pi*self.AIS_eps2 * ( (2./3.)  * self.AIS_slope*(R**3-rc**3)-self.AIS_b0*(R**2-rc**2) )
        self.AIS_Volume[:] = V
        self.AIS_Radius[:] = R

        if self.dbg == 'AIS':
            print('rc init: ', rc)
            print('V init: ', V)
            print('V init: ', V)
            print('R init: ', R)

        return

    def __update_SLR_AIS(self, i):
        if self.dbg==1: print('   SLR AIS: ', i)
        Tg = self.T_anomaly[i]
        R = self.AIS_Radius[i]
        V = self.AIS_Volume[i]
        SL = self.SLR_total[i]
        dSL = SL - self.SLR_total[i-1]
        AIS_includes_dSLais = 1.0

        # define some constants
        Pi = np.pi
        Tf = -1.8
        Toc_0 = 0.72
        Aoc = 3.619e14
        lf    = -1.18

        # From fitting Antarctic air temperature to global mean temperature anomaly
        c1 = 14.28
        c2 = 0.77
        Ta = (Tg - c1) / c2
    
        # Connecting Antarctic ocean temperature to Antarctic air temperature
        Toc = Tf + (self.AIS_a_anto*Tg + self.AIS_b_anto-Tf) / (1. + np.exp(-Tg+(Tf-self.AIS_b_anto)/self.AIS_a_anto))

        if self.dbg == 'AIS': print(i, 'SL:', SL)
        if self.dbg == 'AIS': print(i, 'dSL:', dSL)
        if self.dbg == 'AIS': print(i, 'Ta:', Ta)
        if self.dbg == 'AIS': print(i, 'Toc:', Toc)

        # Start model
        hr   = self.AIS_h0 + self.AIS_c * Ta        # equation 5
        rc   = (self.AIS_b0 - SL)/self.AIS_slope    # application of equation 1 (paragraph after eq3)
        P    = self.AIS_P0 * np.exp(self.AIS_kappa*Ta) # equation 6
        beta = self.AIS_nu * P**(0.5)      # equation 7 (corrected with respect to text)

             
        rR = R - ((hr - self.AIS_b0 + self.AIS_slope*R)**2) / self.AIS_mu
        Btot_1 = P * Pi * R**2 - \
        Pi * beta * (hr - self.AIS_b0 + self.AIS_slope*R) * (R*R - rR*rR) - \
        (4. * Pi * beta * self.AIS_mu**0.5 *   (R-rR)**2.5) / 5.  + \
        (4. * Pi * beta * self.AIS_mu**0.5 * R*(R-rR)**1.5) / 3.
        Btot_2 = P * Pi*R**2
             
        if hr > 0: Btot = Btot_1
        else: Btot = Btot_2

        if self.dbg == 'AIS': print(i, 'Btot:', Btot)
                                  
        F_1   = 0.   # no ice flux
        ISO_1 = 0.   # (third term equation 14) NAME?
        fac_1 = Pi * (1.+self.AIS_eps1) * (4./3. * self.AIS_mu**0.5 * R**1.5 - self.AIS_slope*R**2) # ratio dV/dR (eq 14) 

        ## if R > rc
        # In case there is a marine ice sheet / grounding line
        fac_2   = fac_1 - 2.*Pi*self.AIS_eps2 * (self.AIS_slope*R**2 - self.AIS_b0*R) # correction fac (eq 14)
        Hw = self.AIS_slope*R - self.AIS_b0 + SL  # equation 10
        
        # Ice speed at grounding line (equation 11)
        Speed = self.AIS_f0 * ((1.-self.AIS_alpha) + self.AIS_alpha * ((Toc - Tf)/(Toc_0 - Tf))**2) * \
                (Hw**self.AIS_gamma) / ( (self.AIS_slope*self.AIS_Rad0 - self.AIS_b0)**(self.AIS_gamma-1.) )
        
        F_2   = 2.*Pi*R * self.AIS_dell * Hw * Speed   # equation 9

        # ISO term depends on dSL_tot (third term equation 14 !! NAME)
        c_iso = 2.*Pi*self.AIS_eps2* (self.AIS_slope*rc**2 - (self.AIS_b0/self.AIS_slope)*rc)

        # first term is zero if dSL represents only non-AIS components (dSLais=0)
        # second term is zero if dSL represents all components (dSLais=1)
        ISO_2 = AIS_includes_dSLais * c_iso *  dSL + (1.-AIS_includes_dSLais) * \
                    ((1.-c_iso)/c_iso) * (dSL - lf * (Btot - F_2) / Aoc)
                 
        if R>rc:
            F = F_2
            ISO = ISO_2
            fac = fac_2
        else:
            F = F_1
            ISO = ISO_1
            fac = fac_1


        if Ta > self.AIS_T_crit: disintegration_volume = self.AIS_disintegration_rate * self.AIS_Volume[0] / 57.0
        else: disintegration_volume = 0.0


        self.AIS_Radius[i+1] = R + (Btot-F+ISO-disintegration_volume)/fac
        self.AIS_Volume[i+1] = V + (Btot-F+ISO-disintegration_volume)
        self.SLR_AIS[i+1] = 57. * (1. - self.AIS_Volume[i+1]/self.AIS_Volume[0])

        if self.dbg == 'AIS':
            print(i, 'F:', F)
            print(i, 'ISO:', ISO)
            print(i, 'fac:', fac)

            print(i, 'R: ', self.AIS_Radius[i+1])
            print(i, 'V: ', self.AIS_Volume[i+1])

            print(i, 'SLR_AIS: ', self.SLR_AIS[i+1])

            print()



        return


    def getSLRTotal(self):     return np.copy(self.SLR_total)
    def getSLRThermo(self):    return np.copy(self.SLR_thermo)
    def getSLRLWS(self):       return np.copy(self.SLR_LWS)
    def getSLRMG(self):        return np.copy(self.SLR_MG)
    def getSLRGIS(self):       return np.copy(self.SLR_GIS)
    def getSLRAIS(self):       return np.copy(self.SLR_AIS)
    def getTime(self):         return np.copy(self.time)
    def getTanomaly(self):     return np.copy(self.T_anomaly)
    def getOHCchange(self):    return np.copy(self.OHC_change)
