# -*- coding: utf-8 -*-
"""
Created on Sep 29 00:00:00 2024

@author: Luke P Frash
"""

#imports
import numpy as np
from scipy.stats import linregress as linreg

if __package__ is None or __package__ == '':
    from libs import iogt
    from libs import wells
    from libs import fractures
    from libs import gringarten as gg
    from libs import power
    from libs import economics
    from libs import pickler
    from libs import properties
    from libs.typ import typ
    from libs.units import *
else:
    from .libs import iogt
    from .libs import wells
    from .libs import fractures
    from .libs import gringarten as gg
    from .libs import power
    from .libs import economics
    from .libs import pickler
    from .libs import properties
    from .libs.typ import typ
    from .libs.units import *
water = properties.water()
yr=yr
g=g
mDft=mDft

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

#**********************************************************************************************************
### primary class structure
#**********************************************************************************************************
class core:
    def __init__(self,s,run=True,visuals=False,save=True):
        #key values
        self.pin = s.PIN_PIN_PIN
        self.setup = s
        self.bound = [] #boundary fractures
        self.fracs = [] #shear fractures
        self.activ = [] #flowing fractures
        self.pipes = [] #pipe network
        self.nodes = [] #node network
        self.wells = [] #well network
        self.time = [] #timestamps (ts)
        self.temp = [] #temperatures of all wells (w_T)
        self.enth = [] #enthalpy of all wells (w_h, h2)
        self.ther = [] #thermal energy production
        self.heat = [] #low-temp heat energy production
        self.bulk = [] #gross electricity production
        self.pump = [] #pumping energy demand
        self.enet = [] #net electricity production
        self.warn = [] #warning and error log
        
        #continue
        if run: self.re_init(s,visuals,save)
    
    def re_init(self,s,visuals=True,save=True):
        #patterned development
        self.pad(s,visuals,save)
        
        #save result
        if save: self.save()
    
    def bench(self,s,visuals=True,save=True):
        #**********************************************************************************************************
        ### overrides
        #**********************************************************************************************************
        s.Strategy_SecondOrderNoise_ratio = 0.0 #gringarten does not support 3D heterogeneity
        
        #**********************************************************************************************************
        ### generate system geometry #TODO: develop a tool to parameterize this model using input well geometry
        #**********************************************************************************************************
        #try loading from file
        w = []
        if (s.Strategy_DiscreteWells_path != ''):
            try:
                w = wells.intake(s.Strategy_DiscreteWells_path)
            except:
                w = []
        #recompute key parameters
        if len(w) > 0:
            s.Strategy_Target_type = 'DEEP'
            depth = 0.0
            inj = 0
            pro = 0
            dev = 0.0
            imid = np.zeros(3)
            pmid = np.zeros(3)
            dip = 0.0
            azn = 0.0
            for i in w:
                depth = np.max([depth,i.c0[2],i.c1[2]])
                if i.type == 'injector':
                    inj += 1
                elif i.type == 'producer':
                    pro += 1
                elif i.type == 'perfcluster':
                    imid = 0.5*(i.c0 + i.c1)
                    dip = i.dip
                    azn = i.azn
                elif i.type in ['screen','procluster']:
                    pmid = 0.5*(i.c0 + i.c1)
                    dev = ((i.c0[0]-i.c1[0])**2+(i.c0[1]-i.c1[1])**2+(i.c0[2]-i.c1[2])**2)**0.5
            s.Well_Spacing_m = ((imid[0]-pmid[0])**2+(imid[1]-pmid[1])**2+(imid[2]-pmid[2])**2)**0.5
            s.Well_Dip_rad = dip
            s.Well_Azimuth_rad = azn
            s.Well_DeviatedLength_m = dev
            s.Well_TargetDepth_m = depth
            s.Well_ProducerCount_wells = int(pro+0.001)
            s.re_init()
            s.Well_InjectorCount_wells = int(inj+0.001)
            self.pin = s.PIN_PIN_PIN
            self.setup = s
        #otherwise generate well geometry
        if len(w) < 1:
            w = wells.gen_wells(s)
        if visuals: wells.vtk(w,'%i_wells' %(s.PIN_PIN_PIN),0.005*s.Domain_Size_m)
        #extract injection and production inner diameters (smallest only)
        w_Di = 2.0*(s.Well_ProductionRadius_m - s.Well_WallThickness_m)
        w_Dp = w_Di
        for i in w:
            if i.type == 'perfcluster':
                w_Di = 2.0*i.ra
            if i.type in ['screen','procluster']:
                w_Dp = 2.0*i.ra
                
        #**********************************************************************************************************
        ### input values
        #**********************************************************************************************************
        qinj = s.Circulation_TargetRate_m3ps
        tinj = s.Circulation_InjectionTemperature_K
        w_s = s.Well_Spacing_m
        w_l = s.Well_DrilledLength_m
        w_n = int(s.Well_ProducerCount_wells)
        w_i = int(s.Well_InjectorCount_wells)
        f_n = s.Fracture_Tortuosity_ratio
        f_s = s.Well_Clusters_count
        rockdeep = s.Well_TargetDepth_m
        rocktemp = s.Well_TargetTemp_K
        #w_D = s.Well_ProductionDiameter_m - 2.0*s.Well_WallThickness_m
        #f_rho = s.Fluid_InjectionDensity_kgpm3
        w_f = s.Well_Roughness_metric
        p_whp = s.Circulation_Backpressure_Pa
        s3 = s.Stress_S3_Pa
        f_k = s.Stimulation_ProppantPermeability_m2
        f_e = s.Stimulation_ProppantCompressibility_1pPa
        f_u = s.Stimulation_ProppantConcentration_m3pm3
        Pc = s.Stimulation_CriticalPressure_Pa
        f_R = s.Stimulation_Radius_m
        vol = s.Stimulation_Volume_m3pfrac
        fs = s.Stimulation_Spacing_m
        r_Kt = s.Rock_ThermalConductivity_WpmK
        r_Sv = s.Rock_HeatCapacity_kJpm3K
        f_Sv = s.Fluid_HeatCapacity_kJpm3K
        pp = s.Fluid_ReservoirStatic_Pa
        #leak = s.Stimulation_Leakoff_ratio
        bleed = s.Circulation_Leakoff_ratio
        life = s.Power_LifeSpan_s
        steps = int(s.Domain_TimeSteps_steps)
        #v5 = 1.0/f_rho
        ts = np.linspace(0,life,1+steps)
        warn = 0
        
        #**********************************************************************************************************
        ### solve gringarten
        #**********************************************************************************************************
        #number of directions of flow from the injector to producers
        if s.Strategy_Design_type.upper() in ['O&G','DCM']:
            if w_n > 1:
                split = 2
            else:
                split = 1
            if s.Strategy_Design_type.upper() in ['O&G']:
                f_s = int(2.0*f_s)
                fs = 0.5*fs
        elif s.Strategy_Design_type.upper() in ['EGS','FGS','KGS','CGS']:
            split = w_n
        else:
            split = 1
            print('WARNING: Unhandled design type for Gringarten solver')
        #gringarten with log-distributed flow heterogeneity (i.e., relatively few fast paths defined by tortuosity)
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
        ### producer temperature drop (temperature rise with flow down injector ignored)
        #**********************************************************************************************************
        #additional parameters
        H = s.Power_ConvectionCoefficient_kWpm2K
        segments = s.Well_SegmentCount_units
        dt = ts[1]-ts[0]
        
        #rock conditions and mesh creation
        node_Dz = np.linspace(1,0,segments+1)*s.Well_TargetDepth_m
        node_Tr = node_Dz*s.Rock_ThermalGradient_Kpm + s.Air_Temperature_K
        T_res = grin.Tout[:]
        h_res = h2
        water.set_Pref(P=pp)
        is_rho = 1.0/water.v_from_PT(P=pp,T=tinj)
        mi = qinj*is_rho
        ms = mi*(w_i/w_n)*(1.0-bleed)
        ResSv = s.Rock_HeatCapacity_kJpm3K
        ResKt = s.Rock_ThermalConductivity_WpmK
        CemKt = s.Cement_ThermalConductivity_WpmK
        
        #paramterization
        node_ht = np.zeros((len(ts),segments+1),dtype=float)
        node_Tt = np.zeros((len(ts),segments+1),dtype=float)
        pipe_Lg = node_Dz[:-1]-node_Dz[1:]
        pipe_Rt = np.zeros((len(ts),segments),dtype=float) #R0
        pipe_Et = np.zeros((len(ts),segments),dtype=float)
        pipe_Qt = np.zeros((len(ts),segments),dtype=float)
        
        #well geometry (effect of radius is negligible, so only smallest radius is used)
        Rir = s.Well_ProductionRadius_m+s.Well_AnnularThickness_m
        Ric = s.Well_ProductionRadius_m
        Ris = s.Well_ProductionRadius_m-s.Well_WallThickness_m
        for i in w:
            if i.type in ['screen','procluster']:
                Ris = i.ra
                Ric = i.rb
                Rir = i.rc
        
        #energy vs thermal radius for specified borehole diameter (per unit dT and dL)
        Ror = np.linspace(1,2000,2000)
        ERor = 2*np.pi*ResSv*1.0*((1-np.log(Rir)/np.log(Rir/Ror))*(Ror**2-Rir**2)/2
            + (1/(2*np.log(Rir/Ror)))*(Ror**2*(np.log(Ror)-0.5)-Rir**2*(np.log(Rir)-0.5)))
        lnRor = np.log(Ror[:])
        lnERor = np.log(ERor[:])
        ERm, ERb, r_value, p_value, std_err = linreg(lnERor,lnRor)
        
        #convergence criteria
        goal = 0.05
        goalE = 0.003
        dT0 = 10.0
        
        #iterate over time
        for t in range(0,len(ts)):
            #calculate nodal temperatures by pipe flow and nodal mixing
            iters = 0
            err = np.ones(segments+1)*goal*10
            E0_update_intervals = 5
            max_iters = 20*E0_update_intervals #gives 100 max iterations
            E0e = np.zeros(segments)
            dE0e = np.zeros(segments)
            
            #initial conditions
            node_Tt[t,:] = node_Tr[:]
            
            #boundary conditions
            node_Tt[t,0] = T_res[t]
            node_ht[t,0] = h_res[t]
            
            #iterate for temperature stability
            while 1:
                #calculate nodal fluid enthalpy
                for n in range(0,segments+1):
                    node_ht[t,n] = water.h_from_T(T=node_Tt[t,n])
                    
                #loop breaker
                kerr = np.max(np.abs(err))
                eerr = np.max(np.abs(dE0e)/(np.abs(E0e)+1.0e-9))
                iters += 1
                if iters > max_iters:
                    # print( 'Heat solver converged to %e K after %i iterations' %(kerr,iters-1))
                    break
                elif (kerr < goal) and (eerr < goalE):
                    # print( 'Heat solver converged to %e K after %i iterations' %(kerr,iters-1))
                    break
                
                #follow the flow to estimate heating/cooling of fluid per pipe
                Ap = np.zeros(segments,dtype=float)
                Bp = np.zeros(segments,dtype=float)
                Cp = np.zeros(segments,dtype=float)
                Ep = np.zeros(segments,dtype=float)
                # hm = np.zeros(segments+1,dtype=float)
                # mi = np.zeros(segments+1,dtype=float)
                for p in range(0,segments):
                    #working variables
                    n0 = p
                    n1 = p+1
                    Tr0 = node_Tr[n0] 
                    Tr1 = node_Tr[n1]
                    
                    #dE0: overshoot protection
                    dT = 0.5*(node_Tt[t,n1]+node_Tt[t,n0])-0.5*(node_Tr[n1]+node_Tr[n0])
                    dT = np.max([np.abs(dT0),np.abs(dT)])
                    a = 1.0/(ResSv*np.abs(dT)*ResKt*10**-3)
                    b = 1.0/H
                    c = -2.0*np.abs(dT)*dt
                    dE0 = (-b + (b**2.0 - 4.0*a*c)**0.5)/(2.0*a)
                    
                    #Et: extracted energy, R0: thermal radius, Qt: heat flow
                    E0p = dE0*2.0*np.pi*Rir*pipe_Lg[p]
                    Esp = pipe_Et[t,p]
                    Etp = np.max([E0p,Esp])
                    dE0e[p] = Etp - E0e[p]
                    E0e[p] = Etp
                    pipe_Rt[t,p] = np.exp(ERm*(np.log(np.abs(Etp/(pipe_Lg[p]*dT))))+ERb) + Rir
                    pipe_Qt[t,p] = pipe_Lg[p]/(1.0/(2.0*np.pi*Ris*H) + np.log(pipe_Rt[t,p]/Rir)/(2.0*np.pi*ResKt*10**-3) + np.log(Rir/Ric)/(2.0*np.pi*CemKt*10**-3))
                    
                    #equilibrium enthalpy
                    heq1 = water.h_from_T(T=node_Tr[n1])
                    #non-equilibrium conduction limited heating
                    Ap[p] = pipe_Qt[t,p]*(0.5*(Tr1 + Tr0) - 0.5*(node_Tt[t,n1] + node_Tt[t,n0]))
                    #equilibrium conduction limited heating
                    Bp[p] = pipe_Qt[t,p]*(Tr1 - 0.5*(Tr1 + node_Tt[t,n0]))
                    #flow limited cooling to equilibrium with positive flow (ms > 0)
                    Cp[p] = ms*(heq1 - node_ht[t,p])

                    #take maximum of conduction terms because this will drive conduction heat deliverability
                    Kp = np.asarray([Ap[p],Bp[p]])[np.argmax(np.abs([Ap[p],Bp[p]]))]
                    #if flow limits heat extraction from rock, it will go to equilibrium
                    Qp = Cp[p]
                    #get the limiting term
                    KorQ = np.argmin(np.abs([Kp,Qp]))
                    Ep[p] = np.asarray([Kp,Qp])[KorQ]
                    
                    #conduction limited
                    if KorQ == 0: 
                        hm = Ep[p]/ms + node_ht[t,n0]
                    else:
                        hm = heq1
                    
                    #calculate next nodal temperature
                    z = water.T_from_h(h=hm)
                    err = node_Tt[t,n1] - z
                    node_ht[t,n1] = hm
                    node_Tt[t,n1] = water.T_from_h(h=hm)
                    
            #remember heat exchange
            if t < (len(ts)-1):
                pipe_Et[t+1,:] = pipe_Et[t,:] + np.abs(Ep[:])*dt
        
        #update surface temperature
        is_temp = tinj
        ib_temp = tinj
        pb_temp = grin.Tout
        ps_temp = node_Tt[:,-1]
        for t in range(0,len(h2)):
            h2[t] = water.h_from_PT(P=pp,T=ps_temp[t])
            #error checking
            # print('reservoir temp: %.1f, surface temp: %.1f' %(pb_temp[t], ps_temp[t]))
        
        #**********************************************************************************************************
        ### pressure with thermal siphon #TODO: this module is slow, can it be faster?
        #**********************************************************************************************************
        #compressible propped fracture flow
        bd0p = (vol*f_u/0.64)/(np.pi*f_R**2.0) #m, propped zero-stress aperture assuming close packing with porosity of 1-0.64
        # f_c = f_k*bd0p
        
        #pressure head from frictional losses (temperature effect omitted)
        Ps = np.zeros((4,len(h2)),dtype=float)
        Fs = np.zeros((4,len(h2)),dtype=float)
        
        for t in range(0,len(h2)):
            #static thermal-pressure effect
            p_mu = water.mu_from_T(0.5*(ps_temp[t] + ps_temp[t]))
            ps_rho = 1.0/water.v_from_PT(P=p_whp,T=ps_temp[t])
            pb_rho = 1.0/water.v_from_PT(P=pp,T=pb_temp[t])
            p_bhp = 0.5*(ps_rho+pb_rho)*g*rockdeep
            i_mu = water.mu_from_T(tinj)
            is_rho = 1.0/water.v_from_PT(P=p_whp,T=is_temp)
            ib_rho = 1.0/water.v_from_PT(P=pp,T=ib_temp)
            i_bhp = 0.5*(is_rho+ib_rho)*g*rockdeep
            
            #flow-friction effect
            pn_guess = np.zeros(3)
            pn_guess[0] = p_whp - s3 + p_bhp  #maximum compressive net-stress (compression negative)
            pn_guess[2] = Pc #maximum positive stress (tensile opening)
            pn_guess[1] = 0.5*(pn_guess[2]+pn_guess[0]) #binary search first guess
            for k in range(0,22):
                #calculate hydraulic aperture
                f_p = (f_s*12.0*f_k*bd0p*np.exp(f_e*pn_guess[1]))**(1.0/3.0)
                # f_c = f_k*bd0p*np.exp(f_e*pn_guess[1])
                #pressure drops in fracture(stress), injection well, and production well
                dPf = (2.0+1.0/(w_n/w_i))*12.0*qinj*(0.5*(i_mu+p_mu))/f_p**3.0
                dPwi = (10.7/0.9e-3)*(w_l*i_mu*ib_rho*g*qinj**1.852)/(w_f**1.852*w_Di**4.87)
                dPwp = (10.7/0.9e-3)*(w_l*p_mu*pb_rho*g*(qinj/(w_n/w_i))**1.852)/(w_f**1.852*w_Dp**4.87)
                
                # #error checking
                # if (t == 0) and (k == 21):
                #     print(w_Di, dPwi, w_Dp, dPwp, w_n, w_i)
                
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
        
        #record information into fractures
        f = fractures.gen_hydrofracs(s,w)
        #if visuals: fractures.vtk(f,'%i_hydr' %(s.PIN_PIN_PIN),0.005*s.Domain_Size_m)
        heterogeneity = np.linspace(0,1,int(s.Well_Clusters_count))**3*(1-s.Fracture_Tortuosity_ratio)+s.Fracture_Tortuosity_ratio
        q_fracs = (qinj*heterogeneity**3.0)/np.sum(heterogeneity**3.0)
        roll = -1
        fdex = 0
        bh_var = (12.0*q_fracs*0.5*(i_mu+p_mu)/(dPf/(2.0+1.0/(w_n/w_i))))**(1.0/3.0)
        bd_var = q_fracs*0.5*(i_mu+p_mu)/(f_k*dPf/(2.0+1.0/(w_n/w_i)))
        for a in range(0,len(f)):
            if typ(f[a].typ) == 'propped':
                if (fdex == 0) or (fdex == s.Well_Clusters_count):
                    roll = a
                    fdex = 0
                fdex += 1
                f[a].bh = bh_var[s.Well_Clusters_count-1-(a-roll)]
                f[a].bd = bd_var[s.Well_Clusters_count-1-(a-roll)] #f[a].bd0p*np.exp(f_e*pn_guess[1])
                f[a].bd0p = f[a].bd/np.exp(f_e*pn_guess[1])
                f[a].bd0 = f[a].bd0p
                f[a].Pcen = 0.5*(Ps[2,-1]+Ps[1,-1])
                f[a].Pmax = Ps[1,-1]
                f[a].prop_load = 0.64*f[a].bd0p*np.pi*f_R**2.0
                f[a].stim = s.Stimulation_Steps_count
                f[a].fc = 3.28084*1.01324e15*(f[a].bh**3.0)/12.0
        if visuals: fractures.vtk(f,'%i_hydr' %(s.PIN_PIN_PIN),0.005*s.Domain_Size_m)
        
        #**********************************************************************************************************
        ### seimsicity information
        #**********************************************************************************************************
        #get shear stress at maximum variance
        critf = fractures.gen_critfrac(s)
        M0max = critf[0].tau*critf[0].arup**(3.0/2.0)
        Mwmax = (np.log10(M0max)-9.1)/1.5
        s.Seismicity_MaxQuake_Mw = Mwmax
        s.Seismicity_NumQuake_ea = s.Stimulation_Steps_count*s.Well_Clusters_count
        s.Fractures_ShearStim_ea = 0.0
        s.Fractures_HydroStim_ea = s.Stimulation_Steps_count*s.Well_Clusters_count
        if Fs[2,-1] > 0.0:
            s.Fractures_Hydroprop_ea = s.Stimulation_Steps_count*s.Well_Clusters_count
        else:
            s.Fractures_Hydroprop_ea = 0.0
        s.Cost_Seismic_USD = s.Cost_Seismic_USDpMw*np.exp(Mwmax*s.Cost_Seismic_exp)
        
        #**********************************************************************************************************
        ### power
        #**********************************************************************************************************
        #calculkate power production from bulk efficiency (based on inlet enthlpy) TODO: better solver with air temp
        mi = qinj*w_i*is_rho
        mt = mi*(1.0-bleed)
        P5 = Ps[3]
        ther,bulk,pump,net,heat = power.zarrouk_moon(s,h2,mi,mt,P5)
        
        #record key outputs
        s.Circulation_InjectionPressure_Pa = Ps[3][0]
        s.Circulation_ProductionPressure_Pa = Ps[0][0]
        s.Circulation_FrictionWell_Pa = (Fs[3][0]-Fs[2][0])+(Fs[1][0]-Fs[0][0])
        s.Circulation_FrictionFracture_Pa = Fs[2][0]-Fs[1][0]
        s.Circulation_FractureConductivity_m2m = f_k*bd0p*np.exp(f_e*pn_guess[1])
        s.Circulation_InjectionRate_kgps = mi
        s.Circulation_ProductionRate_kgps = mt
        s.Circulation_InjectionRate_m3ps = mi/is_rho
        s.Circulation_ProductionRate_m3ps = mt/ps_rho
        s.Circulation_Recovery_ratio = mt/mi
        cutoff = len(net)
        k = np.where(net<0)[0]
        if len(k)>0:
            cutoff = k[0] + 1
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
        ### heat recovery factor
        #**********************************************************************************************************
        heat_volume = s.Well_InjectorCount_wells*s.Well_DeviatedLength_m*np.pi*s.Well_Spacing_m**2.0
        heat_stored = heat_volume*(s.Well_TargetTemp_K - s.Circulation_InjectionTemperature_K)*s.Rock_HeatCapacity_kJpm3K
        heat_extrac = np.sum(ther[:cutoff])*dt
        s.Power_RecoveryFactor_ratio = heat_extrac/heat_stored
        
        #**********************************************************************************************************
        ### economics using inflation, debt, equity, and fees to get LCOE and IRR (setup inherits economics' scalar values)
        #**********************************************************************************************************
        #costs
        econ = economics.evaluate(s,w,ts,ther,bulk,pump,heat,visuals=False)
        s.Land_PowerDensity_Wpm2 = 1e3*s.Power_RunningNet_kW/(np.max([s.Well_DeviatedLength_m*np.cos(s.Well_Dip_rad),2.0*s.Stimulation_Radius_m])
                                    *np.max([s.Well_Spacing_m*np.cos(s.Well_RotationPhase_rad)*(s.Well_ProducerCount_wells+s.Well_InjectorCount_wells),2.0*s.Stimulation_Radius_m]))
        
        # #**********************************************************************************************************
        # ### construct output file
        # #**********************************************************************************************************
        # if save:
        #     names = ['Time_Elapsed_yr']
        #     values = [ts/yr]
        #     names += ['Production_Temperature_C']
        #     values += [ps_temp-273.15]
        #     names += ['Production_Enthalpy_kJpkg']
        #     values += [h2]
        #     names += ['Energy_Thermal_kW']
        #     values += [ther]
        #     names += ['Energy_Gross_kW']
        #     values += [bulk]
        #     names += ['Energy_Pump_kW']
        #     values += [pump]
        #     names += ['Energy_Heat_kW']
        #     values += [heat]
        #     names += ['Economics_CapexOpex_USD']
        #     values += [econ.Cost_Capex_USD+econ.Cost_Opex_USD]
        #     names += ['Economics_EnergyCredit_USD']
        #     values += [econ.Sales_Energy_USD+econ.Cost_Energy_USD+econ.Tax_DirectPay_USD]
        #     names += ['Economics_CashFlowTaxed_USD']
        #     values += [econ.Profit_CashFlowTaxed_USD]
        #     names += ['Economics_Equity_USD']
        #     values += [econ.Profit_Equity_USD]
        #     iogt.save(s,values,names)
        
        #**********************************************************************************************************
        ### extra information
        #**********************************************************************************************************
        self.setup = s
        self.econo = econ
        self.bound = fractures.gen_domain(s)
        self.fracs = []
        self.activ = f
        self.faces = self.bound + self.fracs + self.activ 
        self.pipes = []
        self.nodes = []
        self.wells = w
        self.time = ts
        self.temp = ps_temp
        self.enth = h2
        self.ther = ther
        self.heat = heat
        self.bulk = bulk
        self.pump = pump
        self.enet = net
        self.w_T = np.asarray([np.ones(len(ps_temp))*tinj,ps_temp])
        self.w_h = np.asarray([np.ones(len(h2))*water.h_from_PT(P=p_whp,T=tinj),h2])
    
    def save(self):
        names = ['Time_Elapsed_yr']
        values = [self.time/yr]
        names += ['Production_Temperature_C']
        values += [self.temp-273.15]
        names += ['Production_Enthalpy_kJpkg']
        values += [self.enth]
        names += ['Energy_Thermal_kW']
        values += [self.ther]
        names += ['Energy_Gross_kW']
        values += [self.bulk]
        names += ['Energy_Pump_kW']
        values += [self.pump]
        names += ['Energy_Heat_kW']
        values += [self.heat]
        names += ['Economics_CapexOpex_USD']
        values += [self.econo.Cost_Capex_USD+self.econo.Cost_Opex_USD]
        names += ['Economics_EnergyCredit_USD']
        values += [self.econo.Sales_Energy_USD+self.econo.Cost_Energy_USD+self.econo.Tax_DirectPay_USD]
        names += ['Economics_CashFlowTaxed_USD']
        values += [self.econo.Profit_CashFlowTaxed_USD]
        names += ['Economics_Equity_USD']
        values += [self.econo.Profit_Equity_USD]
        iogt.save(s,values,names)
    
    def pad(self,s,visuals=True,save=True):
        #**********************************************************************************************************
        ### multi-bench operations (well count is per bench, pattern includes bi-directional drilling)
        #**********************************************************************************************************
        #timesteps
        nt = 1+int(s.Domain_TimeSteps_steps)
        ts = np.linspace(0,s.Power_LifeSpan_s,nt)
        dt = ts[1]-ts[0]
        
        #counters and trackers
        count = 0
        maxit = 10
        passed = 0.0
        start = 0
        activ = []
        faces = []
        w = []
        f = []
        div = np.zeros(nt)
        temp = np.zeros(nt)
        enth = np.zeros(nt)
        ther = np.zeros(nt)
        heat = np.zeros(nt)
        bulk = np.zeros(nt)
        pump = np.zeros(nt)
        enet = np.zeros(nt)
        w_T = [] #TODO figure this out later
        w_h = [] #TODO figure this out later
        
        #loop until limit reached
        s.Well_Azimuth_rad += -np.pi
        while 1:
            #breakers
            count += 1
            if (count > maxit) or (count > s.Well_Pattern_count):
                break
            if ((passed*60.0) > s.Power_LifeSpan_s):
                break
            
            #modify depth and direction
            if s.Strategy_Design_type.upper() in ['O&G','DCM']:
                #change direction
                s.Well_Azimuth_rad += np.pi
                if s.Well_Azimuth_rad > np.pi*2.0:
                    s.Well_Azimuth_rad += -2.0*np.pi
                #deepen every other count
                if ((count % 2) == 1) and (count > 2):
                    s.Strategy_Target_type = 'DEEP'
                    s.Well_TargetDepth_m += s.Well_Spacing_m
                    s.re_init()
            else:
                if (s.Well_Pattern_count > 1):
                    print('Well patterns not yet implemented for %s' %(s.Strategy_Design_type.upper()))
            print('count %i: azimuth %.2f' %(count,s.Well_Azimuth_rad))
            
            #execute bench
            self.bench(s,visuals,save)
            
            #compounded trackers
            div[start:] += np.ones(nt-start)
            activ += self.activ
            faces = self.bound + self.fracs + activ
            w += self.wells
            temp[start:] += self.temp[:nt-start]
            enth[start:] += self.enth[:nt-start]
            ther[start:] += self.ther[:nt-start]
            heat[start:] += self.heat[:nt-start]
            bulk[start:] += self.bulk[:nt-start]
            pump[start:] += self.pump[:nt-start]
            enet[start:] += self.enet[:nt-start]
            w_T += list(self.w_T)
            w_h += list(self.w_h)
            
            # #additional information
            # Well_Clusters_count += s.Well_Clusters_count
            # Seismicity_MaxQuake_Mw = np.max([Seismicity_MaxQuake_Mw,s.Seismicity_MaxQuake_Mw])
            # Seismicity_NumQuake_ea += s.Seismicity_NumQuake_ea
            # Fractures_ShearStim_ea += s.Fractures_ShearStim_ea
            # Fractures_HydroStim_ea += s.Fractures_HydroStim_ea
            # Fractures_Hydroprop_ea += s.Fractures_Hydroprop_ea
            # Cost_Seismic_USD += s.Cost_Seismic_USD
            # Circulation_InjectionPressure_Pa = np.max([Circulation_InjectionPressure_Pa,s.Circulation_InjectionPressure_Pa])
            # Circulation_ProductionPressure_Pa = np.max([Circulation_ProductionPressure_Pa,s.Circulation_ProductionPressure_Pa])
            # Circulation_FrictionWell_Pa = np.max([Circulation_FrictionWell_Pa,s.Circulation_FrictionWell_Pa])
            # Circulation_FrictionFracture_Pa = np.max([Circulation_FrictionFracture_Pa,s.Circulation_FrictionFracture_Pa])
            # Circulation_FractureConductivity_m2m = np.min([Circulation_FractureConductivity_m2m,s.Circulation_FractureConductivity_m2m])
            # Circulation_InjectionRate_kgps = Circulation_InjectionRate_kgps + s.Circulation_InjectionRate_kgps
            # Circulation_ProductionRate_kgps = Circulation_ProductionRate_kgps +s.Circulation_ProductionRate_kgps
            # Circulation_InjectionRate_m3ps = Circulation_InjectionRate_m3ps + s.Circulation_InjectionRate_m3ps
            # Circulation_ProductionRate_m3ps = Circulation_ProductionRate_m3ps + s.Circulation_ProductionRate_m3ps
            # Circulation_Recovery_ratio = Circulation_ProductionRate_kgps/Circulation_InjectionRate_kgps
            # Power_RecoveryFactor_ratio = np.max([Power_RecoveryFactor_ratio,s.Power_RecoveryFactor_ratio])
            
            # #**********************************************************************************************************
            # ### economics using inflation, debt, equity, and fees to get LCOE and IRR
            # #**********************************************************************************************************
            # #costs
            # econ = economics.evaluate(s,w,ts,ther,bulk,pump,heat,visuals=False)
            # s.Land_PowerDensity_Wpm2 = 1e3*s.Power_RunningNet_kW/(np.max([s.Well_DeviatedLength_m*np.cos(s.Well_Dip_rad),2.0*s.Stimulation_Radius_m])
            #                             *np.max([s.Well_Spacing_m*np.cos(s.Well_RotationPhase_rad)*(s.Well_ProducerCount_wells+s.Well_InjectorCount_wells),2.0*s.Stimulation_Radius_m]))
        
            #move forward time for next iteration
            passed += s.Cost_Drilling_h+s.Cost_Stimulation_h #TODO add factor for rig utilization
            start = np.min([nt,int(60.0*60.0*passed/dt)])
        
        #averaging
        temp = temp/div
        enth = enth/div
        ther = ther/div
        heat = heat/div
        bulk = bulk/div
        pump = pump/div
        enet = enet/div
        
        #power information
        cutoff = nt
        k = np.where(enet<0)[0]
        if len(k)>0:
            cutoff = k[0] + 1
        s.Power_PeakGross_kW = np.max(bulk[0:cutoff])
        s.Power_PeakPump_kW = np.min(pump[0:cutoff])
        s.Power_PeakNet_kW = np.max(enet[0:cutoff])
        s.Power_MinHeat_kW = np.min(heat[0:cutoff])
        running = np.zeros(len(bulk)) + bulk
        running = np.mean(running[0:cutoff])
        s.Power_Running_kW = running
        output = np.zeros(len(enet)) + enet
        output = np.mean(output[0:cutoff])
        s.Power_NetOutput_kW = output
            
        #recompute economics #TODO - be cognizant of capacity expansion CAPEX versus duplicated costs for redundant hardware
        
        #visuals
        if visuals: wells.vtk(w,'%i_wells' %(s.PIN_PIN_PIN),0.003*s.Domain_Size_m)
        if visuals: fractures.vtk(faces,'%i_hydr' %(s.PIN_PIN_PIN),0.003*s.Domain_Size_m)

        
#**********************************************************************************************************
### execute
#**********************************************************************************************************
def RUN(s,report=False,visuals=False,save=True):
    geom = core(s,run=True,visuals=visuals,save=save)
    if save: pickler.packup(geom,geom.pin)
    if report: return geom

#**********************************************************************************************************
### testing
#**********************************************************************************************************
if True:
    s = iogt.setup()
    s.Well_ProducerCount_wells = 4
    s.Well_Pattern_count = 6
    s.Strategy_Design_type = 'DCM'
    s.re_init()
    geom = RUN(s,report=True,visuals=True,save=True)
    # s.Well_ProducerCount_wells = 2
    # s.re_init()
    # geom = RUN(s,report=True,visuals=True)
    # s.Well_ProducerCount_wells = 3
    # s.re_init()
    # geom = RUN(s,report=True,visuals=True)

