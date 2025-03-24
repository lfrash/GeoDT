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

### solve gringarten using the model's input parameters
# - assumes uniform flow
# - includes effects of fracture spacing
# - no natural fractures
# - save result as *.csv
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
    
    print(ti)
    print(rt)
    print(kt)
    print(sv)
    print(nw)
    print(ws)
    print(nf)
    print(fs)
    print(bf)
    print(tt)
    print(pr)
    print(fv)
    print(tp)
    
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
        