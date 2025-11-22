# -*- coding: utf-8 -*-
""" 
Tools for calculating the power production from geothermal energy
Stations:
    1 = production bottomhole (reservoir)
    2 = production wellhead (surface)
    2s = flash steam (surface)
    2l = flash liquid (surface)
    3 = turbine exit/chiller inlet (surface)
    4 = chiller exit/pump inlet (surface)
    5 = injection wellhead (surface)
    6 = injection bottomhole (reservoir)
"""
#initializations
import numpy as np
if __package__ is None or __package__ == '':
    import iogt
    import properties
    from units import *
else:
    from . import iogt
    from . import properties
    from .units import *
water = properties.water()
g=g

#**********************************************************************************************************
### cap heat production based on demand
#**********************************************************************************************************
def waste_heat(s=iogt.setup(),h2=[],mt=100.0):
    heat = np.zeros(len(h2))
    if s.Circulation_InjectionTemperature_K < 353.15: #optimal binary outlet temperature is around 80 C, if injection is cooler, sell the heat
        hi = water.h_from_PT(P=s.Circulation_Backpressure_Pa,T=353.15)
        lo = water.h_from_PT(P=s.Circulation_Backpressure_Pa,T=s.Circulation_InjectionTemperature_K)
        heat = heat+h2
        heat[heat>hi] = hi
        heat = (heat-lo)*mt*s.Power_GeneralEfficiency_ratio
        heat[heat<0] = 0.0
    return heat

#**********************************************************************************************************
### estimate power --- ref: Zarrouk and Moon, 2014 + Thermoflex, 2025
#**********************************************************************************************************
def zarrouk_moon(s=iogt.setup(),h2=[],mi=100.0,mt=100.0,P5=[]):
    #key parameters
    P2 = s.Circulation_Backpressure_Pa
    i_rho = s.Fluid_InjectionDensity_kgpm3
    pump_eff = s.Power_GeneralEfficiency_ratio
    h5 = water.h_from_PT(P=P2,T=s.Circulation_InjectionTemperature_K)
    
    #gross power production from plant inlet enthalpy (no correction for air temperature nor use of only a fraction of the heat's potential)
    effy = 0.095639*np.log(h2) - 0.538659
    bulk = effy*h2*mt #kJ/s
    bulk[bulk<0] = 0.0
    
    #pump losses accounting for thermospihon
    pump = -1e-3*(mi/i_rho)*(P5-P2)/pump_eff
    pump[pump>0] = 0.0
    
    #thermal power
    ther = mt*(h2 - h5)
    
    #net power
    net = bulk + pump #kJ/s
    
    #heat
    heat = waste_heat(s,h2,mt)
    
    #results
    return ther, bulk, pump, net, heat

def zarrouk_moon_old(s=iogt.setup(),h2=[],mi=100.0,mt=100.0,P5=30.0e6):
    #key parameters
    P2 = s.Circulation_Backpressure_Pa
    rho_i = s.Fluid_InjectionDensity_kgpm3
    pump_eff = s.Power_GeneralEfficiency_ratio
    rockdeep = s.Well_TargetDepth_m
    h5 = water.h_from_PT(P=P5,T=s.Circulation_InjectionTemperature_K)
    
    #gross power production from plant inlet enthalpy (no correction for air temperature nor use of only a fraction of the heat's potential)
    effy = 0.095639*np.log(h2) - 0.538659
    bulk = effy*h2*mt #kJ/s
    bulk[bulk<0] = 0.0
    
    #pump losses accounting for thermospihon
    rho_p = np.zeros(len(h2))
    visc_p = np.zeros(len(h2))
    pump = np.zeros(len(h2))
    for t in range(0,len(h2)):
        T = water.T_from_Ph(P=P2,h=h2[t])
        rho_p[t] = 1.0/water.v_from_PT(T=T,P=P2)
        visc_p[t] = water.mu_from_T(T)
        pump[t] = -1e-3*(mi/rho_i)*(0.5+0.5*visc_p[t]/visc_p[0])*((P5-P2)-(rho_i-rho_p[t])*g*rockdeep)/pump_eff
    pump[pump>0] = 0.0
    
    #thermal power
    ther = mt*(h2 - h5)
    
    #net power
    net = bulk + pump #kJ/s
    
    #results
    return rho_p, ther, bulk, pump, net

#testing
if False:
    #Zarrouk and Moon, 2014
    s = iogt.setup()
    s.Circulation_Backpressure_Pa = 10e6
    s.Circulation_InjectionTemperature_K = 60.0+273.15
    s.Power_GeneralEfficiency_ratio = 0.85
    s.Strategy_Target_type = 'DEEP'
    s.re_init()
    h2 = np.asarray([1600.0,1400.0,1200.0,1000.0,800.0,600.0,400.0,200.0])
    mi = 100.0
    mt = 100.0
    P5 = np.zeros(len(h2))
    for t in range(0,len(h2)):
        T = water.T_from_Ph(P=s.Fluid_InjectionStatic_Pa,h=h2[t])
        P5[t] = 50e6-(s.Fluid_InjectionDensity_kgpm3-1/water.v_from_PT(P=s.Circulation_Backpressure_Pa,T=T))*g*s.Well_TargetDepth_m
    ther,bulk,pump,net,heat = zarrouk_moon(s,h2,mi,mt,P5)
    print('enthal (kJ/kg)\tbackpress (kPa)\tthermal (kW)\tbulk (kW)\tpump (kW)\tnet (kW)\theat (kW)')
    for i in range(0,len(h2)):
        print('%.6f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e' %(h2[i],s.Circulation_Backpressure_Pa,ther[i],bulk[i],pump[i],net[i],heat[i]))
        
    # #Zarrouk and Moon, 2014
    # s = iogt.setup()
    # s.Circulation_Backpressure_Pa = 10e6
    # s.Circulation_InjectionTemperature_K = 60.0+273.15
    # s.Power_GeneralEfficiency_ratio = 0.85
    # s.re_init()
    # h2 = np.asarray([1800.0,1600.0,1400.0,1200.0,1000.0,800.0,600.0,400.0])
    # mi = 100.0
    # mt = 100.0
    # P5 = 50e6
    # rho,ther,bulk,pump,net = zarrouk_moon_old(s,h2,mi,mt,P5)
    # print('enthal (kJ/kg)\tbackpress (kPa)\tdensity (kg/m3)\tthermal (kW)\tbulk (kW)\tpump (kW)\tnet (kW)')
    # for i in range(0,len(h2)):
    #     print('%.6f\t%.2e\t%.6f\t%.2e\t%.2e\t%.2e\t%.2e' %(h2[i],s.Circulation_Backpressure_Pa,rho[i],ther[i],bulk[i],pump[i],net[i]))

    