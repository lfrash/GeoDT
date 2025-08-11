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
import matplotlib.pyplot as plt
import pylab

### gringarten solver validation test using Gringarten, 1975, Fig. 5's inputs
# - assumes uniform flow
# - includes effects of fracture spacing
# - no natural fractures
# - plot results
if True:
    #setup
    fig = plt.figure(figsize=(10.0,5.0),dpi=100)
    ax = fig.add_subplot(111)
    ty = np.linspace(0,100,101)
    tp = np.asarray(ty*365.25*24*60*60)
    grin = gg.gringarten()
    #10 fractures, 160 m spacing
    fn = 10; fs = 160.0
    grin.basic(rocktemp=300.0,
               rockKt=2.594,
               rockSv=1.046*2650,
               numwells=2,
               wellspacing=1000.0,
               numfracs=fn,
               fracspacing=fs,
               tinj=65.0,
               bulkflow=0.145,
               tortuosity=1.0,
               poreSv=4.184*1000,
               timepassed=tp)
    ax.plot(ty,grin.Tout,label='%i fractures, %.0f m spacing, %0.2e Xe' %(fn,fs,grin.Xe)) 
    #10 fractures, 80 m spacing
    fn = 10; fs = 80.0
    grin.basic(rocktemp=300.0,
                rockKt=2.594,
                rockSv=1.046*2650,
                numwells=2,
                wellspacing=1000.0,
                numfracs=fn,
                fracspacing=fs,
                tinj=65.0,
                bulkflow=0.145,
                tortuosity=1.0,
                poreSv=4.184*1000,
                timepassed=tp)
    ax.plot(ty,grin.Tout,label='%i fractures, %.0f m spacing, %0.2e Xe' %(fn,fs,grin.Xe)) 
    #10 fracture, 40 m spacing
    fn = 10; fs = 40.0
    grin.basic(rocktemp=300.0,
               rockKt=2.594,
               rockSv=1.046*2650,
               numwells=2,
               wellspacing=1000.0,
               numfracs=fn,
               fracspacing=fs,
               tinj=65.0,
               bulkflow=0.145,
               tortuosity=1.0,
               poreSv=4.184*1000,
               timepassed=tp)
    ax.plot(ty,grin.Tout,label='%i fractures, %.0f m spacing, %0.2e Xe' %(fn,fs,grin.Xe)) 
    #1 fractures, 500 m spacing
    fn = 1; fs = 500.0
    grin.basic(rocktemp=300.0,
               rockKt=2.594,
               rockSv=1.046*2650,
               numwells=2,
               wellspacing=1000.0,
               numfracs=fn,
               fracspacing=fs,
               tinj=65.0,
               bulkflow=0.145,
               tortuosity=1.0,
               poreSv=4.184*1000,
               timepassed=tp)
    ax.plot(ty,grin.Tout,label='%i fractures, %.0f m spacing, %0.2e Xe' %(fn,fs,grin.Xe)) 
    #labels
    ax.set_xlabel('Time (yr)',fontsize=10)
    ax.set_ylabel('Temp (K)',fontsize=10)
    ax.set_title('%0.2e C1, %.02e C2' %(grin.C1,grin.C2))
    ax.legend(loc='lower center',ncol=2,fontsize=10)
    plt.tight_layout()
    pylab.show()

### solve gringarten using the model's input parameters
# - assumes uniform flow
# - includes effects of fracture spacing
# - no natural fractures
# - save result as *.csv
if False:
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
    
    #export vtk
    if True:
        directory = path #os.getcwd()
        try:
            os.makedirs('vtk')
        except:
            pass
        geom.build_vtk('%s/vtk/%i' %(os.getcwd(),geom.pin),[1,1,1,1,1,1]) #get everything [wells, fracs, fnets, nodes, pipes, bound]
    
    #extract key model parameters
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
    #fv = (water.h_from_T(rt)-water.h_from_T(ti))/((rt-ti)*water.v_from_PT(P=pr,T=0.5*(rt-ti)))
    fv = (water.h_from_T(rt)-water.h_from_T(ti))/((rt-ti)*geom.v5)
    tp = geom.ts
    
    # print(ti)
    # print(rt)
    # print(kt)
    # print(sv)
    # print(nw)
    # print(ws)
    # print(nf)
    # print(fs)
    # print(bf)
    # print(tt)
    # print(pr)
    # print(fv)
    # print(tp)
    
    #get geodt temps
    gt_T = np.zeros(len(geom.p_hm)) 
    gt_T = gt_T + geom.p_hm
    for i in range(0,len(gt_T)):
        gt_T[i] = water.T_from_h(gt_T[i])
    
    #run gringarten solution (default)
    grin = gg.gringarten()
    grin.basic(rocktemp=rt,rockKt=kt,rockSv=sv,
              numwells=nw,wellspacing=ws,numfracs=nf,fracspacing=fs,
              tinj=ti,bulkflow=bf,tortuosity=tt,poreSv=fv,
              timepassed=tp)
    
    #run gringarten solution (ideal)
    idea = gg.gringarten()
    idea.basic(rocktemp=rt,rockKt=kt,rockSv=sv,
              numwells=nw,wellspacing=ws,numfracs=nf,fracspacing=fs,
              tinj=ti,bulkflow=bf,tortuosity=1.0,poreSv=fv,
              timepassed=tp)
    
    #save results
    data = np.asarray([grin.time,gt_T,grin.Tout,idea.Tout])
    data = np.core.records.fromarrays(data,names='time_s, geodt_K, gringarten_K, ideal_K',formats='f8,f8,f8,f8')
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
