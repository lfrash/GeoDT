# -*- coding: utf-8 -*-
"""
A script to run geodt solvers in a batch with the inputs structured in a grid
"""

#import packages
import sys
import numpy as np
path = ''
if __package__ is None or __package__ == '':
    # from solvers import geodt as gt
    from solvers import basic as gg
    from solvers.libs import csv
    from solvers.libs import iogt
    path = sys.path[0]
else:
    # from .solvers import geodt as gt
    from .solvers import basic as gg
    from .solvers.libs import csv
    from .solvers.libs import iogt
    path = sys.path[0]+'\\plugins'

#**********************************************************************************************************
### execute batch run of GeoDT with stochastic inputs
#**********************************************************************************************************
if True:
    #randomize inputs using settings file
    try:
        names, data = csv.load_set(path+'\\runme.set')
    except:
        names, data = csv.load_set(path+'\\runme.gts')
        
    #default setup
    setup = iogt.setup()
    for name in names:
        if name in ['Strategy_SourceDirectory_path','Strategy_TargetDirectory_path','Strategy_DiscreteWells_path','Strategy_DiscreteFractures_path',
                    'Strategy_Cluster_type','Strategy_Solver_type','Strategy_Design_type','Strategy_Target_type','Strategy_Stim_type']:
            val = data[name][1]
            setattr(setup,name,val)
        else:
            nom = data[name][1]
            if nom == '': nom = '0.0'
            nom = float(nom)
            setattr(setup,name,nom)
    orig = iogt.setup()
    iogt.transfer_attributes(setup,orig)
    setup.re_init()

    #columns to ignore
    ignore = []
    if 'Strategy_Target_type' in names:
        if data['Strategy_Target_type'][1].upper() in ['DEEP','DEPTH']:
            ignore += ['Well_TargetTemp_K']
        else: #['TEMP','TEMPERATURE']
            ignore += ['Well_TargetDepth_m']
    if 'Strategy_Stim_type' in names:
        if data['Strategy_Stim_type'][1].upper() in ['RAD','RADI','RADIUS']:
            ignore += ['Stimulation_TargetVolume_m3pfrac']
        else: #['VOL','VOLU','VOLUME']
            ignore += ['Stimulation_TargetRadius_m']
    if 'Strategy_Cluster_type' in names:
        if data['Strategy_Cluster_type'][1].upper() in ['NUM','COUNT','NUMBER']:
            ignore += ['Stimulation_TargetSpacing_m']
        else: #['SPA','SPACING','VAR']
            ignore += ['Stimulation_TargetClusters_count']
    
    print('ignore list:')
    print(ignore)
    
    #count the number of randomized inputs
    iters = np.max([1,int(float(data['Strategy_Iterations_units'][1]))])
    dims = 0
    draw = []
    who = []
    for name in names:
        if name.split('_') != 'Strategy':
            if name != 'PIN_PIN_PIN':
                if (data[name][3] != 'none') and (not(name in ignore)):
                    dims += 1
                    who += [name]
                    seed = np.zeros(iters)
                    #set first value as nominal input
                    low = float(data[name][0])
                    nom = float(data[name][1])
                    hig = float(data[name][2])
                    seed[0] = nom
                    #set next values as using structured array
                    if iters > 1:
                        if data[name][3] in ['loguniform']:
                            seed[1:] = 10.0**(np.linspace(0,1,iters-1)*(np.log10(hig)-np.log10(low))+np.log10(low))
                        else: #data[name][3] in ['uniform','normal']:
                            seed[1:] = np.linspace(0,1,iters-1)*(hig-low)+low
                    draw += [seed]
    draw = np.asarray(draw)
    idex = np.zeros(dims,dtype=int)
    total = iters**dims
    print('Running %i variables with %i iterations, giving %i total simulations' %(dims, iters, total))
    
    #prevent absurdities
    total = np.min([1000000,total])
    
    #run all simulations
    idex[0] += -1
    for i in range(0,total):
        #reset to default
        iogt.transfer_attributes(orig,setup)
        
        #calculate index with rollover
        idex[0] += 1
        for j in range(0,dims):
            #check for rollover
            if idex[j] > (iters-1):
                idex[j] = 0
                idex[j+1] += 1
            # print(who[j])
            # print(draw[j])
            # print(idex[j])
            
            #draw sample
            setattr(setup,who[j],draw[j][idex[j]])
        #print(i, idex)
            
        #reinitialize
        setup.re_init()
        
        #run
        if setup.Strategy_Solver_type == 'gt':
            print('*******************************************************************')
            print('                    running realization %i' %(i))
            print('*******************************************************************\n\n')
            # gt.RUN(setup)
        #basic
        else:
            gg.RUN(setup)
else:
    print('Sorry, runme.set was not created as it should have been.')

