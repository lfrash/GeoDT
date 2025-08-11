# -*- coding: utf-8 -*-
"""
Created on Sep 29 00:00:00 2024

@author: Luke P Frash
"""

#imports
import numpy as np
from plugins.solvers import geodt as gt
from plugins.solvers import gringarten as gg
from plugins.solvers.libs import gt_properties as properties
water = properties.water()

#master loop
if True:
    #dependencies
    import tkinter as tk
    from tkinter import filedialog
    import os
    
    #select replay file
    root = tk.Tk()
    root.title('GeoDT')
    root.geometry('250x50')
    print("select a replay file (e.g., 100000000.geodt)")
    filename = filedialog.askopenfilename(initialdir=os.getcwd())
    if filename == '':
        print("no file selected")
    root.destroy()    
    
    #fetch results from replay file
    geom = gt.core()
    geom.pin = int(filename[-15:-6])
    path = filename[:-23]
    print(geom.pin)
    gt.unpack(geom,geom.pin)
    # gt.unpack(geom,filename=filename)
    
    ########################################################################################################################################
    ### export vtk files for this run
    ########################################################################################################################################
    if True:
        directory = path #os.getcwd()
        try:
            os.makedirs('vtk')
        except:
            pass
        geom.build_vtk('%s/vtk/%i' %(os.getcwd(),geom.pin),[1,1,1,1,1,1]) #get everything [wells, fracs, fnets, nodes, pipes, bound]
    
    ########################################################################################################################################
    ### solve gringarten using the model's input parameters
    # - assumes uniform flow
    # - includes effects of fracture spacing
    # - no natural fractures
    # - save result as *.csv
    ########################################################################################################################################
    if True:
        ti = geom.rock.TinjK
        rt = geom.rock.BH_T
        kt = geom.rock.ResKt
        sv = geom.rock.ResSv
        nw = geom.rock.w_count + 1
        ws = geom.rock.w_spacing
        nf = geom.rock.setup.Stimulation_Clusters_clusters
        if nf == 1:
            fs = 1e9
        else:
            fs = 2.0*geom.rock.w_proportion*geom.rock.w_length/(nf*2.0-1.0)
        bf = geom.rock.setup.Circulation_InjectionRate_m3ps
        tt = geom.rock.f_roughness[1]
        pr = geom.rock.BH_P
        water.set_Pref(P=pr)
        fv = (water.h_from_T(rt)-water.h_from_T(ti))/((rt-ti)*geom.v5)
        tp = geom.ts
        
        #get geodt temps
        gt_T = np.zeros(len(geom.p_hm)) 
        gt_T = gt_T + geom.p_hm
        for i in range(0,len(gt_T)):
            gt_T[i] = water.T_from_h(gt_T[i])
        
        #run gringarten solution
        grin = gg.gringarten()
        grin.basic(rocktemp=rt,rockKt=kt,rockSv=sv,
                  numwells=nw,wellspacing=ws,numfracs=nf,fracspacing=fs,
                  tinj=ti,bulkflow=bf,tortuosity=tt,poreSv=fv,
                  timepassed=tp)
        
        #save results
        data = np.asarray([grin.time,gt_T,grin.Tout])
        data = np.core.records.fromarrays(data,names='time_s, geodt_K, gringarten_K',formats='f8,f8,f8')
        header = data.dtype.names
        gt.build_csv(header,data,filename='gringarten_%i.csv' %(geom.pin),append=False)
            
        # #more advanced functions are a work in progress
        # #test internal data
        # for f in range(0,len(geom.faces)):
        #     print(gt.typ(geom.faces[f].typ))
        # for p in range(0,geom.pipes.num):
        #     print(gt.typ(geom.pipes.typ[p]))
        #     rho = geom.v5
        #     if gt.typ(geom.pipes.typ[p]) in ['propped']: #'choke','propped'
        #         L = geom.pipes.L[p]
        #         W = geom.pipes.W[p]
        #         q = geom.pipes.q[p]
        #         print('index %i, length %.3e, width %.3e, q %.3e' %(p,L,W,q))
        
        # #test gringarten
        # from plugins.solvers import gringarten as gg
        # Darcy = gg.Darcy
        # MPa = gg.MPa
        # gg_solver = gg.gringarten()
        # gg_solver.EGS(numfracs=50,fracspacing=50,fracaperture=0.001,
        #               wellspacing=110,welldiameter=8*0.0254,roughness=80.0,
        #               depth=2300,gradient=80.0,homogeneity=0.1,numwells=2,
        #               proppantconductivity=50*Darcy,backpressure=1*MPa)
    
    ########################################################################################################################################
    ### solve injection flowrate vs pressure using the model's input parameters
    # - computes geomechanics (i.e., stress dependent permeability)
    # - injection pressure at wellhead
    # - injection rate at wellhead
    # - save result as *.csv
    ########################################################################################################################################
    if False:
        print("compute injection rate from pressure")
        #get injection pressure and rate information
        pmax = (70.0*3)*gt.MPa
        pmin = pmax/1000.0
        points = 40
        tips = 10.0**np.linspace(np.log10(pmin),np.log10(pmax),points) + geom.rock.BH_P + geom.rock.p_whp
        visuals = False
        Pis, Qis, Vis, Pcl, Pch, Psc, Pba, Ps3 = geom.rate_vs_pressure(tips=tips,points=points,visuals=visuals)
        Pi = np.max(Pis,axis=1)
        Qi = np.sum(Qis,axis=1)
        
        # #get intermediate point information
        # p_clusters = []
        # p_screens = []
        # for i in range(0,geom.pipes.num):
        #     if gt.typ(geom.pipes.typ[i]) == 'perfcluster':
        #         p_clusters += [geom.nodes.p[geom.pipes.n1]]
        #     elif gt.typ(geom.pipes.typ[i]) == 'screen':
        #         p_screens += [geom.nodes.p[geom.pipes.n1]]
        
        
        
        #save results
        data = np.asarray([range(0,len(Qi)),Qi,Pi,Vis,Pcl,Pch,Psc, Pba, Ps3])
        data = np.core.records.fromarrays(data,names='Point_Count_Num,Injection_Rate_m3ps, Injection_Pressure_Pa, Stored_Volume_m3, ' +
                                          'Entry_Pressure_Pa, Fracture_Pressure_Pa, Production_Pressure_Pa, ' +
                                          'Wellhead_Pressure_Pa, Stress_Opening_Pa',formats='i8,f8,f8,f8,f8,f8,f8,f8,f8')
        header = data.dtype.names
        gt.build_csv(header,data,filename='injectionprofile_%i.csv' %(geom.pin),append=False)
    
    ########################################################################################################################################
    ### generate graphics to illustrate modeling process
    # - well drilling
    # - stimulation
    # - production
    ########################################################################################################################################
    if True:
        if False:
            #wiggle fractures
            for i in range(6,len(geom.faces)):
                geom.faces[i].dia = geom.faces[i].dia * (1.0 + (np.random.uniform(0,1)-0.5)*0.2)
                #geom.faces[i].c0[1:] = geom.faces[i].c0[1:] + (geom.faces[i].dia * (np.random.uniform(0,1,3)[1:]-0.5)*0.1)
            
            print('video creation')
            #general vtk
            directory = path #os.getcwd()
            try:
                os.makedirs('vtk/video')
            except:
                pass
            geom.build_vtk('%s/vtk/video/%i' %(os.getcwd(),geom.pin),[1,1,1,1,1,1]) #get temps [wells, fracs, fnets, nodes, pipes, bound]
        
        #drilling
        
        #stimulation back tracking
        if False:
            #original values
            Rf = []
            Hf = []
            count = -1
            scale = 1.2
            #stimulate toe to heel
            for i in range(6,len(geom.faces)):
                stim = geom.faces[i].stim
                if stim > 0:
                    #first build
                    count += 1
                    geom.build_vtk('%s/vtk/video/%i_f%i' %(os.getcwd(),geom.pin,count),[0,0,1,0,0,0]) #get temps [wells, fracs, fnets, nodes, pipes, bound]
                    #modify radii but skew growth to be slower at larger radius
                    init = geom.faces[i].dia/(scale**stim)
                    for s in range(0,stim):
                        count += 1
                        diff = init*(scale**(stim-s+1) - scale**(stim-s))
                        geom.faces[i].bh = geom.faces[i].bh * (geom.faces[i].dia - diff)/(geom.faces[i].dia)
                        geom.faces[i].dia = geom.faces[i].dia - diff
                        geom.build_vtk('%s/vtk/video/%i_f%i' %(os.getcwd(),geom.pin,count),[0,0,1,0,0,0]) #get temps [wells, fracs, fnets, nodes, pipes, bound]
                    #final radius
                    if geom.faces[i].typ == 2: #nullify hydrofracs
                        count += 1
                        geom.faces[i].bh = geom.faces[i].bh * (geom.faces[i].dia - diff)/(geom.faces[i].dia)
                        geom.faces[i].dia = geom.rock.rc
                        geom.build_vtk('%s/vtk/video/%i_f%i' %(os.getcwd(),geom.pin,count),[0,0,1,0,0,0]) #get temps [wells, fracs, fnets, nodes, pipes, bound]
        
        #temp gradients
        if False:
            spacing = 50.0
            geom.nodes.T = geom.nodes.Tr
            geom.build_pts(spacing,'%s/vtk/video/%i_xback' %(os.getcwd(),geom.pin),'x')
            geom.build_pts(spacing,'%s/vtk/video/%i_yback' %(os.getcwd(),geom.pin),'y')
            geom.build_pts(spacing,'%s/vtk/video/%i_zback' %(os.getcwd(),geom.pin),'b')
            geom.build_pts(spacing,'%s/vtk/video/%i_xfron' %(os.getcwd(),geom.pin),'t')
            geom.build_pts(spacing,'%s/vtk/video/%i_zsurf' %(os.getcwd(),geom.pin),'s')
        
        #initial temps
        if False:
            spacing = 3.3
            geom.nodes.T = geom.nodes.Tr
            geom.pipes.R0 = geom.Rt[0,:]
            geom.build_vtk('%s/vtk/video/%i_t%i' %(os.getcwd(),geom.pin,0),[0,0,0,1,0,0]) #get temps [wells, fracs, fnets, nodes, pipes, bound]
            geom.build_pts(spacing,'%s/vtk/video/%i_t%i' %(os.getcwd(),geom.pin,0),'z')
            #production
            for i in range(0,len(geom.ts)-1):
                print('production %i' %(i))
                
                #backtrack nodal temperature information
                geom.nodes.T = geom.Tt[i,:]
                geom.pipes.R0 = geom.Rt[i,:]
                #produce vtk and heat map
                geom.build_vtk('%s/vtk/video/%i_t%i' %(os.getcwd(),geom.pin,i+1),[0,0,0,1,0,0]) #get temps [wells, fracs, fnets, nodes, pipes, bound]
                geom.build_pts(spacing,'%s/vtk/video/%i_t%i' %(os.getcwd(),geom.pin,i+1),'z')
                # if (i == len(geom.ts)-2):
                #     geom.build_pts(spacing,'%s/vtk/video/%i_t%i' %(os.getcwd(),geom.pin,i+1))
            
        #scale up
        if True:
            arr = np.zeros(3)
            s = geom.rock.w_spacing
            c0s = [arr + np.asarray([0.0, 2*s, 0.0]),
                   arr + np.asarray([0.0, -2*s, 0.0]),
                   arr + np.asarray([0.0, 1.5*s, -s]),
                   arr + np.asarray([0.0, -1.5*s, -s]),
                   ]
            c0i = []
            for s in range(0,len(c0s)):
                for i in range(0,len(geom.wells)):
                    if s == 0:
                        c0i += [geom.wells[i].c0]
                    geom.wells[i].c0 = c0i[i] + c0s[s]
                    geom.build_vtk('%s/vtk/video/%i_w%i' %(os.getcwd(),geom.pin,s+1),[1,0,0,0,0,0]) #get temps [wells, fracs, fnets, nodes, pipes, bound]
                
        
        
        
        
        
        
        