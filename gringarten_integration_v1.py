# -*- coding: utf-8 -*-
"""
Created on Sep 29 00:00:00 2024

@author: Luke P Frash
"""

#imports
import numpy as np

# # from plugins.solvers import geodt as gt
from plugins.solvers.libs import stats
from plugins.solvers.libs import csv
from plugins.solvers.libs import iogt
from plugins.solvers.libs import cashflow
from plugins.solvers import gringarten as gg
from plugins.solvers.libs import properties
water = properties.water()
import matplotlib.pyplot as plt
import pylab

from plugins.solvers.libs.units import *
yr=yr

"""
### gringarten solver for P10/P50/P90 with heterogeneity
# - assume non uniform flow
# - include effects of fracture spacing
# - no natural fractures
# - consider optimization of power production
# - plot results
if False:
    #dependencies
    import tkinter as tk
    from tkinter import filedialog
    import os
    
    #select settings file
    root = tk.Tk()
    root.title('GeoDT')
    root.geometry('250x50')
    print("select a settings file (*.gds)")
    filename = filedialog.askopenfilename(initialdir=os.getcwd())
    if filename == '':
        print("no file selected")
    root.destroy()
    
    #load settings
    names, data = csv.load_csv(filename)
    setup = iogt.setup()
    
    #randomize inputs
    def sample():
        for name in names:
            if name in ['Strategy_SourceDirectory_path','Strategy_TargetDirectory_path','Strategy_DiscreteFractures_path',
                        'Strategy_Solver_type','Strategy_Design_type','Strategy_Target_type','Strategy_Stim_type']:
                val = data[name][1]
                setattr(setup,name,val)
            else:
                low = data[name][0]
                if low == '': low = '0.0'
                low = float(low)
                nom = data[name][1]
                if nom == '': nom = '0.0'
                nom = float(nom)
                hig = data[name][2]
                if hig == '': hig = '0.0'
                hig = float(hig)
                dis = data[name][3]
                if dis == 'none':
                    val = nom
                elif dis == 'uniform':
                    val = np.random.uniform()*(hig-low)+low
                elif dis == 'normal':
                    val = stats.norm_trunc(1,mu=0.5*(hig+low),
                                        dev=0.25*(hig-low),
                                        lo=low,hi=hig)[0]
                elif dis == 'loguniform':
                    val = 10.0**(np.random.uniform()*(np.log10(hig)-np.log10(low))+np.log10(low))
            setattr(setup,name,val)
            # if name in ['Strategy_TargetDirectory_path']:
            #     print('inital working directory %s' %(os.getcwd()))
            #     os.chdir(data[name][1])
            #     print('modified working directory %s' %(os.getcwd()))
    
else:
    setup = iogt.setup()
    
if True:
    #iterate and solve
    iters = int(float(setup.Strategy_Iterations_units))        
    for i in range(0,iters): #range(0,1): #
        #sample parameters
        sample()"""

if True:
    if True:
        setup = iogt.setup()
        filename='test.gts'
        i = 0
        
        #intake key values
        qinj = setup.Circulation_InjectionRate_m3ps
        tinj = setup.Circulation_InjectionTemperature_K
        w_s = setup.Well_Spacing_m
        w_n = int(setup.Well_ProducerCount_wells)
        if setup.Strategy_Design_type == "O&G":
            w_n = np.min([2.0,w_n])
        f_n = setup.Fracture_TortuosityMin_ratio
        f_s = setup.Stimulation_Clusters_clusters
        
        # if setup.Strategy_Target_type in ['DEPTH','DEEP']:
        #     rocktemp = setup.Well_TargetDepth_m*setup.Rock_ThermalGradient_Kpm + setup.Air_Temperature_K
        #     rockdeep = setup.Well_TargetDepth_m
        # else:
        #     rocktemp = setup.Well_TargetTemp_K
        #     rockdeep = (setup.Well_TargetTemp_K - setup.Air_Temperature_K)/setup.Rock_ThermalGradient_Kpm

        rockdeep = setup.Well_TargetDepth_m
        rocktemp = setup.Well_TargetTemp_K

        g = 9.81 #m/s2
        w_R = setup.Well_HydraulicRadius_m #m, well radius (hydraulic)
        f_rho = setup.Fluid_NominalDensity_kgpm3 #kg/m3, fluid density
        f_mu = setup.Fluid_Viscosity_Pas #Pa-s, water viscosity
        w_f = setup.Well_Roughness_metric #metric, Hazen-Williams pipe roughness - varies from 150 (smooth) to 80 (rough)
        p_whp = setup.Circulation_Backpressure_Pa #Pa, production pressure at wellhead
        r_E = setup.Rock_YoungsModulus_Pa #Pa, Young's modulus
        r_v = setup.Rock_PoissonRatio_ratio #Poisson's ratio
        r_rho = setup.Rock_Density_kgpm3 #kg/m3, rock density
        r_Kn = setup.StressMin_Coefficient_ratio #Pa/Pa, normal stress to vertical stress ratio
        f_k = setup.Proppant_InitialPermeabilityNom_m2 #m2, proppant permeability
        f_e = setup.Proppant_CompressibilityNom_1pPa #1/Pa, proppant compressibility
        f_u = setup.Stimulation_ProppantConcentration_m3pm3 #m3/m3, proppant concentration
        f_i = setup.Domain_PressureIncrement_Pa #Pa, stimulation solver pressure increment
        
        #fracture geometry parameters #TODO - ensure that the Pc calculation is sensible
        Pc = f_i*20.0 #fractures start at 0.2*spacing and grow by 1.2x per step, meaning 10 growth steps
        f_R = 0.0 #m, fracture radius
        f_V = 0.0 #m3, stimulation volume for all stages
        f_w = 0.0 #m, fracture aperture
        vol = 0.0 #m3, hydropropped volume
        stimtype = setup.Strategy_Stim_type.upper()
        if stimtype in ['RAD','R','RADI','RADIUS']:
            f_R = setup.Stimulation_TargetRadius_m #m
            mult = np.max([0.9,f_R/(0.2*w_s)])
            step = np.log(mult/0.5787)/0.1823
            #Pc = f_i*step
            f_w = 8.0*Pc*(1.0-r_v**2.0)*f_R/(np.pi*r_E)
            vol = (4.0/3.0)*np.pi*f_R**2.0*0.5*f_w
            f_V = vol*f_s
        else: # ['V','VOL','VOLUME']
            f_V = setup.Stimulation_TargetVolume_m3
            vol = np.max([1e-9,f_V/f_s])
            vol0 = (4.0/3.0)*np.pi*(0.2*w_s)**2.0*0.5*8.0*2.0*f_i*(1.0-r_v**2.0)*(0.2*w_s)/(np.pi*r_E)
            step = 1.6384*np.log(vol/vol0)+1.4738
            #Pc = f_i*step
            f_R = ((3.0/16.0)*vol*r_E/(Pc*(1.0-r_v**2.0)))**(1.0/3.0)
            f_w = 8.0*Pc*(1.0-r_v**2.0)*f_R/(np.pi*r_E)
        
        #compensate for leakoff #TODO make this less crappy
        Stimulation_Leakoff_ratio = 4.55
        vol = vol*Stimulation_Leakoff_ratio
        
        #don't let well spacing be larger than fracture radius #TODO - recompute fracture dimensions?
        if w_s > (f_R/1.1):
            w_s = f_R/1.1
            
        #solve gringarten
        fs = setup.Well_DeviatedLength_m*setup.Well_ProducerProportion_ratio/f_s
        water.set_Pref(P=70e6)
        v5 = water.v_from_PT(T=tinj,P=70e6)
        fv = (water.h_from_T(rocktemp)-water.h_from_T(tinj))/((rocktemp-tinj)*v5)
        if i == 0:
            tp = np.linspace(0,setup.Power_LifeSpan_s,31)
        grin = gg.gringarten()
        grin.basic(rocktemp=rocktemp,
               rockKt=setup.Rock_ThermalConductivity_WpmK,
               rockSv=setup.Rock_HeatCapacity_kJpm3K,
               numwells=w_n+1,
               wellspacing=w_s,
               numfracs=f_s,
               fracspacing=fs,
               tinj=tinj,
               bulkflow=qinj,
               tortuosity=f_n,
               poreSv=fv,
               timepassed=tp)
        
        #compressible propped fracture flow (taken from circuit.py)
        prop = vol*f_u #m3, proppant volume per fracture
        mass = prop*1538.0 #kg, proppant mass per fracture
        bd0p = (vol*f_u/0.64)/(np.pi*f_R**2.0) #m, propped zero-stress aperture assuming close packing with porosity of 1-0.64
        f_c = f_k*bd0p
        w_l = rockdeep + setup.Well_DeviatedLength_m*(1.0-0.5*np.sin(setup.Well_Dip_rad)) #m, well length
        bhp = rockdeep*f_rho*g #Pa, bottomhole pressure
        s3 = rockdeep*r_rho*g*r_Kn #Pa, minimum principal stress
        p_pro = p_whp - s3 + bhp #Pa, outlet pressure adjusted for hydrostatic effects
        
        def propped(qi):
            #setup binary search
            pn_guess = np.zeros(3)
            pn_guess[0] = p_pro #maximum compressive stress (compression negative)
            pn_guess[2] = Pc #maximum positive stress (tensile opening)
            pn_guess[1] = 0.5*(pn_guess[2]+pn_guess[0]) #binary search first guess
            pn_calc = 0.0
            for k in range(0,22):
                #calculate hydraulic aperture
                f_p = (f_s*12.0*f_k*bd0p*np.exp(f_e*pn_guess[1]))**(1.0/3.0)
                f_c = f_k*bd0p*np.exp(f_e*pn_guess[1])
                #pressure drops
                dPf = 12.0*qi*f_mu/f_p**3.0
                dPwi = (10.7/0.9e-3)*(w_l*f_mu*f_rho*g*qi**1.852)/(w_f**1.852*(2*w_R)**4.87)
                dPwp = (10.7/0.9e-3)*(w_l*f_mu*f_rho*g*(qi/w_n)**1.852)/(w_f**1.852*(2*w_R)**4.87)
                #head pressures
                P4 = p_pro #production surface
                P3 = P4 + dPwp #production bottomhole
                P2 = P3 + (2.0+1.0/w_n)*dPf #injection bottomhole
                P1 = P2 + dPwi #injection wellhead
                pn_calc = P2 #0.5*(P2+P3) #net pressure, use higher to account for minor hydropropping
                head = np.asarray([P1,P2,P3,P4])
                #next guess
                if pn_calc > pn_guess[1]:
                    pn_guess[0] = pn_guess[1]
                else:
                    pn_guess[2] = pn_guess[1]
                pn_guess[1] = 0.5*(pn_guess[2]+pn_guess[0])
            #absolute pressures
            real = head + s3
            real[0] = real[0] - bhp
            real[-1] = real[-1] - bhp
            return head, real, f_c
        
        #flow and pressures
        head, real, f_c = propped(qinj)
        p_fric = head[0]-head[-1]
        rho_i = 1.0/water.v_from_PT(T=tinj,P=real[1])
        rho_p = 1.0/water.v_from_PT(T=rocktemp,P=real[2])
        p_siph = (rho_i-rho_p)*g*rockdeep #TODO: add fluid density calculations
        p_loss = np.max([0.0,p_fric-p_siph])
        real[0] = real[0]-p_siph 
        real[1] = real[1]-p_siph
        real[2] = real[2]-p_siph #density difference from injector to producer

        #**********************************************************************************************************
        ### estimate power
        #**********************************************************************************************************
        mt = qinj/v5 #kg/s
        injectors = 1
        producers = int(setup.Well_ProducerCount_wells)
        if setup.Strategy_Design_type == "O&G":
            injectors = np.max([producers-1,1])
            mt = mt*injectors
            qinj = qinj*injectors
        h2 = np.zeros(len(tp))
        h5 = water.h_from_PT(real[0],tinj)
        for t in range(0,len(h2)):
            h2[t] = water.h_from_PT(real[-1],grin.Tout[t])
        effy = 0.095639*np.log(h2) - 0.538659 #(Zarrouk and Moon, 2014: Geothermics + Thermoflex, 2025)
        Bulk = effy*h2*mt #kJ/s
        Bulk[Bulk<0] = 0.0
        # Pump = -1e-3*np.ones(len(Bulk))*qinj*p_loss/setup.Power_GeneralEfficiency_ratio
        rho_i = 1.0/water.v_from_PT(T=tinj,P=real[1])
        rho_p = np.zeros(len(Bulk))
        for t in range(0,len(h2)):
            rho_p[t] = 1.0/water.v_from_PT(T=grin.Tout[t],P=real[2]) #TODO: improve calculation to account for thermal viscosity effects
        Pump = -1e-3*qinj*(p_fric-(rho_i-rho_p)*g*rockdeep)/setup.Power_GeneralEfficiency_ratio
        Pump[Pump>0] = 0.0
        Ther = mt*(h2 - h5)
        Net = Bulk + Pump #kJ/s
        # if len(Net[Net>0]) > 0:
        #     Pos = np.asarray(Net)[Net>0]
        #     Gross = np.mean(np.asarray(Bulk)[Net>0])
        # else:
        #     Pos = np.asarray([0.0])
        #     Gross = np.asarray([0.0])
        
        #**********************************************************************************************************
        ### calcuate drilling costs
        #**********************************************************************************************************
        
        #estimate drilling cost
        Economics_DrillingHourly_USDph = 5082.1 #$/hr
        Economics_CasingHourly_USDph = 4001.3 #$/hr
        Economics_DrillingCutVolume_USDpm4 = 4.5e-5 #$/m4
        Economics_DrillingCutOffset_USDpm3 = 1000.0 #$/m3
        Economics_CasingCutVolume_USDpm4 = 1.4 #$/m4
        Economics_CasingCutOffset_USDpm3 = 1460.0 #$/m3
        Economics_MobilizationFixed_USD = 5130000.0 #USD
        Economics_MobilizationTime_h = 336.0 #hr
        Economics_RigWalkFixed_USD = 250000.0 #USD
        Economics_RigWalkTime_h = 48.0 #hr
        Economics_WellheadTime_h = 46.0 #hr
        Economics_DrillingTimeVolume_hpm4 = 2.3e-11 #h/m4
        Economics_DrillingTimeOffset_hpm3 = 0.5 #h/m3
        Economics_CasingTimeVolume_hpm4 = 2.0e-12 #h/m4
        Economics_CasingTimeOffset_hpm3 = 0.15 #h/m3
        Economics_LearningRate_ratio = 0.20 #learning rate
        Strategy_Completion_type = 'DC' #type of completion
        Special_ConductorRad_m = 0.190 #0.39/2 #m, swap to compensate for simple model
        
        #for number of wells
        # Learning_Step_num += 1
        
        #running totals
        Drilling_Time_h = Economics_MobilizationTime_h
        Drilling_Cost_USD = Economics_MobilizationFixed_USD
        Learning_Factor_hph = 1.0
        Cost_Tracker = []
        
        
        #estimated cost for drilling campaign
        for w in range(0,injectors+producers):
            #learning factor
            Learning_Factor_hph = (w+1)**(np.log(1-Economics_LearningRate_ratio)/np.log(2))
            if w > 0:
               Drilling_Time_h += Economics_RigWalkTime_h
               Drilling_Cost_USD += Economics_RigWalkFixed_USD 
                    
            #vertical segment cost
            Vertical_Depth = rockdeep - setup.Well_DeviatedLength_m*0.5*np.sin(setup.Well_Dip_rad)
            Vertical_MidDepth = 0.5*Vertical_Depth
            Vertical_Volume = Vertical_Depth*np.pi*setup.Well_BoreRadius_m**2 #TODO - have Rodney make better
            Vertical_Volume = Vertical_Depth*np.pi*Special_ConductorRad_m**2
            Vertical_DCV = (Economics_DrillingCutVolume_USDpm4*Vertical_MidDepth**2.0 + Economics_DrillingCutOffset_USDpm3) * Vertical_Volume
            Drilling_Cost_USD += Vertical_DCV
            Vertical_DTV = (Economics_DrillingTimeVolume_hpm4*Vertical_MidDepth**3.0 + Economics_DrillingTimeOffset_hpm3) * Vertical_Volume
            Drilling_Time_h += Vertical_DTV*Learning_Factor_hph
            Drilling_Cost_USD += Vertical_DTV*Learning_Factor_hph*Economics_DrillingHourly_USDph
            if True: #cased
                Vertical_CCV = (Economics_CasingCutVolume_USDpm4*Vertical_MidDepth**1.0 + Economics_CasingCutOffset_USDpm3) * Vertical_Volume
                Drilling_Cost_USD += Vertical_CCV
                Vertical_CTV = ((Economics_CasingTimeVolume_hpm4*Vertical_MidDepth**3.0 + Economics_CasingTimeOffset_hpm3) * Vertical_Volume)
                Drilling_Time_h += Vertical_CTV*Learning_Factor_hph
                Drilling_Cost_USD += Vertical_CTV*Learning_Factor_hph*Economics_CasingHourly_USDph
            
            #deviated segement cost
            Deviated_MidDepth = Vertical_Depth + 0.5*setup.Well_DeviatedLength_m #measured depth
            Deviated_Volume = setup.Well_DeviatedLength_m*np.pi*setup.Well_BoreRadius_m**2
            Deviated_DCV = (Economics_DrillingCutVolume_USDpm4*Deviated_MidDepth**2.0 + Economics_DrillingCutOffset_USDpm3) * Deviated_Volume
            Drilling_Cost_USD += Deviated_DCV
            Deviated_DTV = (Economics_DrillingTimeVolume_hpm4*Deviated_MidDepth**3.0 + Economics_DrillingTimeOffset_hpm3) * Deviated_Volume
            Drilling_Time_h += Deviated_DTV*Learning_Factor_hph
            Drilling_Cost_USD += Deviated_DTV*Learning_Factor_hph*Economics_DrillingHourly_USDph
            if (Strategy_Completion_type in ['FEM','FRACTURE EXCHANGE']) or (w < injectors): #cased
                Deviated_CCV = (Economics_CasingCutVolume_USDpm4*Deviated_MidDepth**1.0 + Economics_CasingCutOffset_USDpm3) * Deviated_Volume
                Drilling_Cost_USD += Deviated_CCV
                Deviated_CTV = (Economics_CasingTimeVolume_hpm4*Deviated_MidDepth**3.0 + Economics_CasingTimeOffset_hpm3) * Deviated_Volume
                Drilling_Time_h += Deviated_CTV*Learning_Factor_hph
                Drilling_Cost_USD += Deviated_CTV*Learning_Factor_hph*Economics_CasingHourly_USDph
            
            #non batch-drilling costs
            if w_n == 0:
                Drilling_Time_h +=   6*Economics_WellheadTime_h
                Drilling_Cost_USD += 6*Economics_WellheadTime_h*Economics_DrillingHourly_USDph
            Cost_Tracker += [[Drilling_Time_h,Drilling_Cost_USD]]
        
        if False:
            Cost_Tracker = np.asarray(Cost_Tracker)
            fig = pylab.figure(figsize=(10.0,5.0),dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(np.linspace(1,len(Cost_Tracker),len(Cost_Tracker)),Cost_Tracker[:,0]/24,color='blue')
            # plt.colorbar(ax=ax1,mappable=sc,label=name[2],orientation='vertical')
            ax.set_xlabel('Number of Wells (wells)',fontsize=10)
            ax.set_ylabel('Drilling Time (d)',fontsize=10)
            plt.tight_layout()
            
            fig = pylab.figure(figsize=(10.0,5.0),dpi=100)
            ax = fig.add_subplot(111)
            ax.plot(np.linspace(1,len(Cost_Tracker),len(Cost_Tracker)),Cost_Tracker[:,1]*1e-6,color='green')
            # plt.colorbar(ax=ax1,mappable=sc,label=name[2],orientation='vertical')
            ax.set_xlabel('Number of Wells (wells)',fontsize=10)
            ax.set_ylabel('Drilling Cost ($M USD)',fontsize=10)
            plt.tight_layout()
            
            pylab.show()
        
        #**********************************************************************************************************
        ### calcuate stimulation costs
        #**********************************************************************************************************
                
        #estimate stimulation cost
        Economics_FracFixed_USD = 1233149.0 #$
        Economics_FracDepth_USD = 638.0 #$/m
        Economics_FracHourly_USDph = 5973.3 #$/hr
        Economics_FracSand_USDpkg = 0.662 #$/kg #ISP
        # Economics_FracSand_USDpkg = 0.248 #$/kg #Sand
        Economics_FracFuel_USDpkWh = 0.451 #$/kWh
        Economics_FracWater_USDpm3 = 22.2 #$/m3
        Economics_FracPressureFactor_ratio = 0.8
        Economics_FracPressureFactor_scale = 1e-8 #Pa
        Economics_FracEquipment_USDpm3 = 9.85 #$/m3
        Economics_FracTimeFactor_hph = 1.95
        Well_FracStringRadius_m = 0.075 #m
        Well_FracStringRadius_m = w_R #m
        
        #information
        Pstim = (10.7/0.9e-3)*(w_l*f_mu*f_rho*g*setup.Stimulation_InjectionRate_m3ps**1.852)/(w_f**1.852*(2*Well_FracStringRadius_m)**4.87) + s3+setup.Tension_Cohesion_Pa-bhp
        Stimulation_Cost_USD = Economics_FracFixed_USD
        Stimulation_Cost_USD += Economics_FracDepth_USD*rockdeep
        Stimulation_Time_h = Economics_FracTimeFactor_hph*((vol*f_s)/setup.Stimulation_InjectionRate_m3ps)/(60*60)*injectors
        Stimulation_Cost_USD += Economics_FracEquipment_USDpm3*vol*f_s*(Economics_FracPressureFactor_ratio+(Economics_FracPressureFactor_scale*Pstim)**3)*injectors
        Stimulation_Cost_USD += Economics_FracFuel_USDpkWh*Stimulation_Time_h*Pstim*setup.Stimulation_InjectionRate_m3ps*1e-3
        Stimulation_Cost_USD += Economics_FracWater_USDpm3*vol*f_s*injectors
        Stimulation_Cost_USD += Economics_FracSand_USDpkg*mass*f_s*injectors
        Stimulation_Time_h += Economics_WellheadTime_h
        Stimulation_Cost_USD += Stimulation_Time_h*Economics_FracHourly_USDph
        
        #**********************************************************************************************************
        ### Cost of capital economics
        #**********************************************************************************************************
        #power plant cost
        Power_PlantCost_USD = np.max(Bulk)*setup.Economics_EquipmentCost_USDpkW
        
        #economics solver configuration
        econ = cashflow.cashflow()
        
        if True: #Durham #TODO note that this can be changed
            if True: #Ingore stimulation cost for fault based systems #TODO find a better way to embed this type of thing
                if i==0: print('notice: cost of stimulation ignored')
                Stimulation_Cost_USD = 0.0
                Stimulation_Time_h = 0.0
            if True: #Ignore injection well cost for production-only systems
                if i==0: print('notice: cost of injection well ignored')
                Drilling_Cost_USD = Drilling_Cost_USD - Cost_Tracker[0][1]
                Drilling_Time_h = Drilling_Time_h - Cost_Tracker[0][0]
                p_fric = head[-2]-head[-1]
                Pump = -1e-3*qinj*(p_fric-(rho_i-rho_p)*g*rockdeep)/setup.Power_GeneralEfficiency_ratio
                Pump[Pump>0] = 0.0
        
            econ.Currency_Exchange_USDpGBP = 1.3
            econ.Time_Start_yr = 2025 #yr, calendar year for project initation
            econ.Time_Planning_yr = 1 #yr, final investment decision
            econ.Time_Construction_yr = np.max([2,1+int((Drilling_Time_h+Stimulation_Time_h)/(364.75*24))]) #yr, construction time
            econ.Time_Lifespan_yr = int(econ.Time_Construction_yr+econ.Time_Planning_yr+tp[-1]/(364.75*24*60*60)) #yr, total project life
            econ.Capex_Reservoir_USD = -1.0*(Drilling_Cost_USD+Stimulation_Cost_USD)/econ.Currency_Exchange_USDpGBP #USD, drilling and stimulation costs
            econ.Capex_Facilities_USD = -Power_PlantCost_USD/econ.Currency_Exchange_USDpGBP #USD, power plant cost
            
            #performance, noting that length of these series is shorter than Capex array: Ther, Net, Bulk, Pump #kJ/s
            gap = list(np.zeros(econ.Time_Planning_yr+econ.Time_Construction_yr))
            dt_yr = (tp[1]-tp[0])/(365.25*24*60*60)
            if (dt_yr > 1.05) or (dt_yr < 0.95):
                print("unhandled error: timesteps for power analysis are incompatible with economics - must be 1 yr timesteps")
            heat = econ.Power_Thermal_kWh[2:-2]/(3.0*(tp[1]-tp[0])/(60*60)) #kW, default assumes 3 injectors, 6 producers, so take 1/3
            mh = heat/(setup.Power_GeneralEfficiency_ratio*h2)
            mh[mh<0] = 0.0
            me = mt-mh
            me[me<0] = 0.0
            heat = np.min([heat,Ther],axis=0)
            heat[heat<0] = 0.0
            econ.Power_Gross_kWh =   np.asarray(gap + list(0.5*(Bulk[1:]+Bulk[:-1])*me[1:]/mt))*(tp[1]-tp[0])/(60*60)
            econ.Power_Thermal_kWh = np.asarray(gap + list(0.5*(heat[1:]+heat[:-1])))*(tp[1]-tp[0])/(60*60)
            econ.Power_Pumping_kWh = np.asarray(gap + list(0.5*(Pump[1:]+Pump[:-1])))*(tp[1]-tp[0])/(60*60)
                
        else: #United States
            econ.Currency_Exchange_USDpGBP = 1.0 #per USD, currency exchange rate
            econ.Time_Start_yr = 2025 #yr, calendar year for project initation
            econ.Time_Planning_yr = 1 #yr, final investment decision
            econ.Time_Construction_yr = np.max([2,1+int((Drilling_Time_h+Stimulation_Time_h)/(364.75*24))]) #yr, construction time
            econ.Time_Lifespan_yr = int(econ.Time_Construction_yr+econ.Time_Planning_yr+tp[-1]/(364.75*24*60*60)) #yr, total project life
            econ.Capex_Planning_USD = -3.5e6 #USD, survey & planning costs
            econ.Capex_Reservoir_USD = -1.0*(Drilling_Cost_USD+Stimulation_Cost_USD)/econ.Currency_Exchange_USDpGBP #USD, drilling and stimulation costs
            econ.Capex_Gridconnect_USD = -1.0*(10*1e6)+3e6 #USD, grid connection costs
            econ.Capex_Facilities_USD = -Power_PlantCost_USD/econ.Currency_Exchange_USDpGBP #USD, power plant cost
            econ.Grant_Government_USD = 3.0e6 #USD, government grants
            econ.Opex_Rate_ratio = 0.02 #ratio, opex vs capex
            econ.Sales_ElectricWholesale_USDpkWh = 0.0967 #USD, market rate sales for electricity
            econ.Sales_ElectricContract_USDpkWh = 0.120 #USD, contract sales rate for electricity
            econ.Sales_ElectricPurchase_USDpkWh = 0.0967 #USD, typical purchase price for customers
            econ.Sales_HeatWholesale_USDpkWh = 6.75/293.1 #USD, market rate sales for heat
            econ.Sales_HeatContract_USDpkWh = 6.75/293.1 #USD, contract rate sales for heat
            econ.Sales_CapacityFactor_ratio = 0.95 #ratio, capacity factor for plant uptime
            econ.Sales_ContractTime_yr = 35 #yr, length of contract for heat and electricity
            econ.Equity_CostOfCapital_ratio = 0.12 #ratio, time value of money versus alternative investments
            econ.Equity_Contribution_ratio = 0.40 #ratio, equity contribution assuming debt service reserve account funded by debt 
            econ.Equity_CompanyShare_ratio = 1.00 #ratio, share of equity held by primary company
            econ.Equity_SeriesAShare_ratio = 0.08 #ratio, share of equity held by series A investors
            econ.Equity_SeriesBShare_ratio = 0.90 #ratio, share of equity held by series B investors
            econ.Equity_BrokerFee_ratio = 0.050 #ratio, fee to setup project funds
            econ.Debt_InterestRate_ratio = 0.080 #ratio, interest rate on debt and loans 
            econ.Debt_ArrangementFee_ratio = 0.015 #ratio, fee to setup debt accounts
            econ.Debt_ServiceAccount_ratio = 0.50 #principal and interest fraction of a year
            econ.Debt_Tenor_yr = np.min([30,int(tp[-1]/(364.75*24*60*60))]) #years, loan duration until maturity
            econ.Inflation_CPI_ratio = 0.02 #inflation rate
            econ.Tax_Corporate_ratio = 0.25-0.0341 #ratio, tax rate on profits after losses
            econ.Tax_ProductionCredit_USDpkWh = 0.0035 #0.015 #USD, tax credit for energy production
            econ.Tax_CapitalCredit_ratio = 0.50  #0.40 #USD, tax credit for capital construction
            econ.Tax_DirectPay_ratio = 1.0 #1.0 #ratio, ratio of tax credits that are elligble for direct pay (e.g., Inflation Reduction Act)
            econ.Tax_DepreciationRate_ratio = 0.20
            econ.CO2_Electric_kgpkWh = 0.22499 #kgCO2/kWhe, CO2 emitted from 'grid' power
            econ.CO2_Heat_kgpkWh = 0.18293 #kgCO2/kWht, CO2 emitted from natural gas for heat 
            
            #performance, noting that length of these series is shorter than Capex array: Ther, Net, Bulk, Pump #kJ/s
            gap = list(np.zeros(econ.Time_Planning_yr+econ.Time_Construction_yr))
            dt_yr = (tp[1]-tp[0])/(365.25*24*60*60)
            if (dt_yr > 1.05) or (dt_yr < 0.95):
                print("unhandled error: timesteps for power analysis are incompatible with economics - must be 1 yr timesteps")
            heat = np.zeros(len(Bulk))
            if tinj < 353.15: #optimal binary outlet temeperature is around 80 C, if injection is cooler, sell the heat
                hi = water.h_from_PT(P=real[0],T=353.15)
                lo = water.h_from_PT(P=real[0],T=tinj)
                heat = heat+h2
                heat[heat>hi] = hi
                heat = (heat-lo)*mt*setup.Power_GeneralEfficiency_ratio
                heat[heat<0] = 0.0
            econ.Power_Gross_kWh =   np.asarray(gap + list(0.5*(Bulk[1:]+Bulk[:-1])))*(tp[1]-tp[0])/(60*60)
            econ.Power_Thermal_kWh = np.asarray(gap + list(0.5*(heat[1:]+heat[:-1])))*(tp[1]-tp[0])/(60*60)
            econ.Power_Pumping_kWh = np.asarray(gap + list(0.5*(Pump[1:]+Pump[:-1])))*(tp[1]-tp[0])/(60*60)
        
        #calculated parameters
        econ.re_init()

        #**********************************************************************************************************
        ### Save Results
        #**********************************************************************************************************
        out=[]
        out+=[['pin',i]]
        out+=[['Strategy_Target_type',setup.Strategy_Target_type]]
        out+=[['Strategy_Stim_type',stimtype]]
        out+=[['Strategy_Completion_type',Strategy_Completion_type]]
        out+=[['Well_TargetDepth_m',rockdeep]]
        out+=[['Well_DeviatedLength_m',setup.Well_DeviatedLength_m]]
        out+=[['Well_Dip_rad',setup.Well_Dip_rad]]
        out+=[['Well_Roughness_metric',w_f]]
        out+=[['Well_HydraulicRadius_m',w_R]]
        out+=[['Well_BoreRadius_m',setup.Well_BoreRadius_m]]
        out+=[['Well_FracstringRadius_m',Well_FracStringRadius_m]]
        out+=[['Drilling_BottomTemperature_K',rocktemp]]
        out+=[['Drilling_BottomPressure_MPa',bhp*1e-6]]
        out+=[['Rock_ThermalGradient_Kpm',setup.Rock_ThermalGradient_Kpm]]
        out+=[['Rock_Density_kgpm3',r_rho]]
        out+=[['Rock_ThermalConductivity_WpmK',setup.Rock_ThermalConductivity_WpmK]]
        out+=[['Rock_HeatCapacity_kJpm3K',setup.Rock_HeatCapacity_kJpm3K]]
        out+=[['Rock_PoissonRatio_ratio',r_v]]
        out+=[['Rock_YoungsModulus_Pa',r_E]]
        out+=[['StressMin_Coefficient_ratio',r_Kn]]
        out+=[['Stress_s3_MPa',s3*1e-6]]
        out+=[['Well_ProducerCount_wells',setup.Well_ProducerCount_wells]]
        out+=[['Well_Spacing_m',w_s]]
        out+=[['Stimulation_InjectionRate_m3ps',setup.Stimulation_InjectionRate_m3ps]]
        out+=[['Stimulation_InjectionPressure_MPa',Pstim*1e-6]]
        out+=[['Stimulation_Clusters_clusters',f_s]]
        out+=[['Stimulation_TargetRadius_m',f_R]]
        out+=[['Stimulation_TargetVolume_m3',f_V]]
        out+=[['Stimulation_FractureSpacing_m',fs]]
        out+=[['Stimulation_ProppantConcentration_m3pm3',f_u]]
        out+=[['Stimulation_TotalSandMass_kg',mass*f_s]]
        out+=[['Stimulation_TotalWaterVolume_m3',vol*f_s]]
        out+=[['Stimulation_NetPressure_Pa',Pc]]
        out+=[['Stimulation_Leakoff_ratio',Stimulation_Leakoff_ratio]]
        out+=[['Domain_PressureIncrement_Pa',f_i]]
        out+=[['Proppant_InitialPermeabilityNom_m2',f_k]]
        out+=[['Proppant_CompressibilityNom_1pPa',f_e]]
        out+=[['Circulation_InjectionTemperature_K',tinj]]
        out+=[['Circulation_Backpressure_Pa',p_whp]]
        out+=[['Injection_Rate_m3ps',qinj]]
        out+=[['Injection_Pressure_MPa',real[0]*1e-6]]
        out+=[['Injection_BottomPressure_MPa',real[1]*1e-6]]
        out+=[['Production_BottomPressure_MPa',real[2]*1e-6]]
        out+=[['Production_Pressure_MPa',real[3]*1e-6]]
        out+=[['Fracture_TortuosityMin_ratio',f_n]]
        out+=[['Fracture_Conductivity_m2m',f_c]]
        out+=[['Fluid_HeatCapacity_kJpm3K',fv]]
        out+=[['Fluid_NominalDensity_kgpm3',f_rho]]
        out+=[['Fluid_Viscosity_Pas',f_mu]]
        out+=[['Fluid_Friction_Pa',p_fric]]
        out+=[['Fluid_Siphon_Pa',p_siph]]
        out+=[['Fluid_PumpLoss_Pa',p_loss]]
        out+=[['Power_GeneralEfficiency_ratio',setup.Power_GeneralEfficiency_ratio]]
        out+=[['Stimulation_Time_h',Stimulation_Time_h]]
        out+=[['Drilling_Time_h',Drilling_Time_h]]
        out+=[['Stimulation_Cost_USD',Stimulation_Cost_USD]]
        out+=[['Drilling_Cost_USD',Drilling_Cost_USD]]
        out+=[['Power_Net_kW',np.mean(Net)]]
        out+=[['Power_Bulk_kW',np.mean(Bulk)]]
        out+=[['Power_Pump_kW',np.mean(Pump)]]
        out+=[['Power_Thermal_kW',np.mean(Ther)]]
        for var in vars(econ):
            if np.isscalar(getattr(econ,var)):
                out += [[var,getattr(econ,var)]]
        for var in vars(econ):  
            if var == 'Profit_CashFlowTaxed_USD':
                for t in range(0,len(tp)):
                    out+=[['CFAT %.2e' %(t+1), econ.Profit_CashFlowTaxed_USD[t]]]
        for t in range(0,len(tp)):
            out+=[['Temp %.2e' %(tp[t]/yr), grin.Tout[t]]]
        for t in range(0,len(tp)):
            out+=[['Pout %.2e' %(tp[t]/yr), Net[t]]]
        # gt.save_csv(out=out,filename='gringartenTaiwan1.csv',append=True) #'gringartenl12.csv'
        csv.save_csv(out=out,filename=(filename[:-4] + '.csv'),append=True) #'gringartenl12.csv'

