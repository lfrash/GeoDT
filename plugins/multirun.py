# -*- coding: utf-8 -*-
"""
A script to run geodt solvers in a batch with stochastic inputs
"""

#import packages
import sys
path = ''
if __package__ is None or __package__ == '':
    # from solvers import geodt as gt
    from solvers import basic as gg
    from solvers.libs import csv
    from solvers.libs import iogt
    from solvers.libs import stats
    path = sys.path[0]
else:
    # from .solvers import geodt as gt
    from .solvers import basic as gg
    from .solvers.libs import csv
    from .solvers.libs import iogt
    from .solvers.libs import stats
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
    setup = iogt.setup()
    def sample():
        for name in names:
            if name in ['Strategy_SourceDirectory_path','Strategy_TargetDirectory_path','Strategy_DiscreteWells_path','Strategy_DiscreteFractures_path',
                        'Strategy_Cluster_type','Strategy_Solver_type','Strategy_Design_type','Strategy_Target_type','Strategy_Stim_type']:
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
                    val = stats.uniform(1,low,hig)[0]
                elif dis == 'normal':
                    val = stats.norm_trunc(1,mu=0.5*(hig+low),
                                        dev=0.25*(hig-low),
                                        lo=low,hi=hig)[0]
                elif dis == 'loguniform':
                    val = stats.loguniform(1,low,hig)[0]
            setattr(setup,name,val)
        setup.re_init()
    
    counter = 0
    iters = int(float(data['Strategy_Iterations_units'][1]))
    if 'Strategy_Iterations_units' in names:
        iters = int(float(data['Strategy_Iterations_units'][1]))
        for i in range(0,iters):
            #randomize inputs
            sample()
            #standard
            if setup.Strategy_Solver_type == 'gt':
                print('*******************************************************************')
                print('                    running realization %i' %(i))
                print('*******************************************************************\n\n')
                # gt.RUN(setup)
            #basic
            else:
                gg.RUN(setup)
    else:
        #randomize inputs
        sample()
        #standard
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

