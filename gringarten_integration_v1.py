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
from plugins.solvers.libs import wells
from plugins.solvers.libs import fractures
from plugins.solvers.libs import power
from plugins.solvers.libs import economics
from plugins.solvers.libs import gringarten as gg
from plugins.solvers.libs import properties
water = properties.water()
import matplotlib.pyplot as plt
import pylab

from plugins.solvers.libs.units import *
yr=yr
g=g
mDft=mDft

warn = 0

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
        s = iogt.setup()
        filename='test.gts'
        i = 0
        
        #**********************************************************************************************************
        ### overrides
        #**********************************************************************************************************
        s.Strategy_SecondOrderNoise_ratio = 0.0 #gringarten does not support 3D heterogeneity
        
        #**********************************************************************************************************
        ### input values
        #**********************************************************************************************************
        qinj = s.Circulation_InjectionRate_m3ps
        tinj = s.Circulation_InjectionTemperature_K
        w_s = s.Well_Spacing_m
        w_l = s.Well_DrilledLength_m
        w_n = int(s.Well_ProducerCount_wells)
        w_i = int(s.Well_InjectorCount_wells)
        f_n = s.Fracture_Tortuosity_ratio
        f_s = s.Stimulation_Clusters_clusters
        rockdeep = s.Well_TargetDepth_m
        rocktemp = s.Well_TargetTemp_K
        w_D = s.Well_ProductionDiameter_m - 2.0*s.Well_WallThickness_m
        f_rho = s.Fluid_InjectionDensity_kgpm3
        f_mu = s.Fluid_InjectionViscosity_Pas
        w_f = s.Well_Roughness_metric
        p_whp = s.Circulation_Backpressure_Pa
        r_E = s.Rock_YoungsModulus_Pa
        r_v = s.Rock_PoissonRatio_ratio
        r_rho = s.Rock_Density_kgpm3
        s3 = s.Stress_S3_Pa
        f_k = s.Stimulation_ProppantPermeability_m2
        f_e = s.Stimulation_ProppantCompressibility_1pPa
        f_u = s.Stimulation_ProppantConcentration_m3pm3
        Pc = s.Stimulation_CriticalPressure_Pa
        f_R = s.Stimulation_TargetRadius_m
        vol = s.Stimulation_TargetVolume_m3pfrac
        f_w = s.Stimulation_ApertureCritical_m
        fs = s.Stimulation_Spacing_m
        r_Kt = s.Rock_ThermalConductivity_WpmK
        r_Sv = s.Rock_HeatCapacity_kJpm3K
        f_Sv = s.Fluid_HeatCapacity_kJpm3K
        pp = s.Fluid_ReservoirStatic_Pa
        s_rho = s.Stimulation_ProppantDensity_kgpm3
        leak = s.Stimulation_Leakoff_ratio
        life = s.Power_LifeSpan_s
        steps = int(s.Domain_TimeSteps_steps)
        v5 = 1.0/f_rho
        
        #**********************************************************************************************************
        ### generate system geometry
        #**********************************************************************************************************
        #well geometry
        w = wells.gen_wells(s)
        wells.vtk(w,'%i_wells' %(s.PIN_PIN_PIN),0.005*s.Domain_Size_m)
        
        #fracture geometry
        f = fractures.gen_hydrofracs(s,w)
        fractures.vtk(f,'%i_hydr' %(s.PIN_PIN_PIN),0.005*s.Domain_Size_m)
        
        #**********************************************************************************************************
        ### solve gringarten
        #**********************************************************************************************************
        #number of directions of flow from the injector to producers
        if s.Strategy_Design_type.upper() in ['O&G']:
            if w_n > 1:
                split = 2
            else:
                split = 1
        elif s.Strategy_Design_type.upper() in ['EGS','FGS','KGS','CGS']:
            split = w_n
        else:
            split = 1
            print('WARNING: Unhandled design type for Gringarten solver')
        #gringarten with log-distributed flow heterogeneity (i.e., relatively few fast paths defined by tortuosity)
        if i == 0:
            ts = np.linspace(0,life,1+steps)
        grin = gg.gringarten()
        grin.basic(rocktemp=rocktemp,
               rockKt=r_Kt,
               rockSv=r_Sv,
               flowsplit=split,
               wellspacing=w_s,
               numfracs=f_s,
               fracspacing=fs,
               tinj=tinj,
               bulkflow=qinj,
               tortuosity=f_n,
               poreSv=f_Sv,
               timepassed=ts)
        #TODO: add well thermal losses
        #get produced enthalpy, density, and viscosity
        h2 = np.zeros(len(ts))
        for t in range(0,len(h2)):
            h2[t] = water.h_from_PT(P=pp,T=grin.Tout[t])

        #**********************************************************************************************************
        ### pressure and flow with thermal siphon
        #**********************************************************************************************************
        #compressible propped fracture flow
        prop = vol*f_u #m3, proppant volume per fracture
        mass = prop*s_rho #kg, proppant mass per fracture
        bd0p = (vol*f_u/0.64)/(np.pi*f_R**2.0) #m, propped zero-stress aperture assuming close packing with porosity of 1-0.64
        f_c = f_k*bd0p
        
        #pressure head from frictional losses (temperature effect omitted)
        Ps = np.zeros((4,len(h2)),dtype=float)
        Fs = np.zeros((4,len(h2)),dtype=float)
        
        for t in range(0,len(h2)):
            #static thermal-pressure effect
            p_mu = water.mu_from_T(grin.Tout[t])
            ps_rho = 1.0/water.v_from_PT(P=p_whp,T=grin.Tout[t])
            pb_rho = 1.0/water.v_from_PT(P=pp,T=grin.Tout[t])
            p_bhp = 0.5*(ps_rho+pb_rho)*g*rockdeep
            i_mu = water.mu_from_T(tinj)
            is_rho = 1.0/water.v_from_PT(P=p_whp,T=tinj)
            ib_rho = 1.0/water.v_from_PT(P=pp,T=tinj)
            i_bhp = 0.5*(is_rho+ib_rho)*g*rockdeep
            
            #flow-friction effect
            pn_guess = np.zeros(3)
            pn_guess[0] = p_whp - s3 + p_bhp  #maximum compressive net-stress (compression negative)
            pn_guess[2] = Pc #maximum positive stress (tensile opening)
            pn_guess[1] = 0.5*(pn_guess[2]+pn_guess[0]) #binary search first guess
            for k in range(0,22):
                #calculate hydraulic aperture
                f_p = (f_s*12.0*f_k*bd0p*np.exp(f_e*pn_guess[1]))**(1.0/3.0)
                f_c = f_k*bd0p*np.exp(f_e*pn_guess[1])
                #pressure drops in fracture(stress), injection well, and production well
                dPf = (2.0+1.0/(w_n/w_i))*12.0*qinj*(0.5*(i_mu+p_mu))/f_p**3.0
                dPwi = (10.7/0.9e-3)*(w_l*i_mu*ib_rho*g*qinj**1.852)/(w_f**1.852*w_D**4.87)
                dPwp = (10.7/0.9e-3)*(w_l*p_mu*pb_rho*g*(qinj/(w_n/w_i))**1.852)/(w_f**1.852*w_D**4.87)
                #head pressures
                Fps = p_whp - s3 + p_bhp #production surface
                Fpb = Fps + dPwp #production bottomhole
                Fib = Fpb + dPf #injection bottomhole
                Fis = Fib + dPwi #injection wellhead
                #next guess with net pressure from injector because of hydropropping
                if Fib > pn_guess[1]:
                    pn_guess[0] = pn_guess[1]
                else:
                    pn_guess[2] = pn_guess[1]
                pn_guess[1] = 0.5*(pn_guess[2]+pn_guess[0])
            
            #absolute pressures
            Pps = Fps + s3 - p_bhp
            Ppb = Fpb + s3
            Pib = Fib + s3
            Pis = Fis + s3 - i_bhp
            
            #no thermal blowout
            if (Pis < 1.02e5):
                if warn == 0:
                    print('Warning: thermal blowout scenario detected (i.e., flashing in the well)')
                    warn = 1
                Pis = 1.02e5
            
            #record result
            Ps[0,t] = Pps
            Ps[1,t] = Ppb
            Ps[2,t] = Pib
            Ps[3,t] = Pis
            Fs[0,t] = Fps
            Fs[1,t] = Fpb
            Fs[2,t] = Fib
            Fs[3,t] = Fis
            
        #information
        p_fric = Fs[3,0] - Fs[0,0]
        p_loss = Ps[3,0] - Ps[0,0]
        
        #**********************************************************************************************************
        ### power
        #**********************************************************************************************************
        #calculkate power production from bulk efficiency (based on inlet enthlpy) TODO: better solver with air temp
        mt = qinj*w_i/v5
        mi = mt*(1.0+leak)
        P5 = Ps[3]
        ther,bulk,pump,net,heat = power.zarrouk_moon(s,h2,mi,mt,P5)
        
        #record key outputs
        cutoff = len(net)
        k = np.where(net<0)[0]
        if len(k)>0:
            cutoff = k[0]
        s.Power_PeakGross_kW = np.max(bulk[0:cutoff])
        s.Power_PeakPump_kW = np.min(pump[0:cutoff])
        s.Power_PeakNet_kW = np.max(net[0:cutoff])
        s.Power_MinHeat_kW = np.min(heat[0:cutoff])
        running = np.zeros(len(bulk)) + bulk
        running = np.mean(running[0:cutoff])
        s.Power_Running_kW = running
        output = np.zeros(len(net)) + net
        output = np.mean(output[0:cutoff])
        s.Power_NetOutput_kW = output
        
        #**********************************************************************************************************
        ### economics using inflation, debt, equity, and fees to get LCOE and IRR
        #**********************************************************************************************************
        #costs
        econ = economics.evaluate(s,w,ts,ther,bulk,pump,heat,visuals=False)
        
        #**********************************************************************************************************
        ### construct output file
        #**********************************************************************************************************
        names = ['Time_Elapsed_yr']
        values = [ts/yr]
        names = ['Production_Temperature_C']
        values = [grin.Tout-273.15]
        names = ['Production_Enthalpy_kJpkg']
        values = [h2]
        names += ['Energy_Thermal_kW']
        values += [ther]
        names += ['Energy_Gross_kW']
        values += [bulk]
        names += ['Energy_Pump_kW']
        values += [pump]
        names += ['Energy_Heat_kW']
        values += [heat]
        names += ['Economics_CapexOpex_USD']
        values += [econ.Cost_Capex_USD+econ.Cost_Opex_USD]
        names += ['Economics_EnergyCredit_USD']
        values += [econ.Sales_Energy_USD+econ.Cost_Energy_USD+econ.Tax_DirectPay_USD]
        names += ['Economics_CashFlowTaxed_USD']
        values += [econ.Profit_CashFlowTaxed_USD]
        names += ['Economics_Equity_USD']
        values += [econ.Profit_Equity_USD]
        iogt.save(s,values,names)


