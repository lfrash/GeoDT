#import packages
import sys
path = ''
if __package__ is None or __package__ == '':
    from solvers import geodt as gt
    path = sys.path[0]
else:
    from .solvers import geodt as gt
    path = sys.path[0]+'\\plugins'
import numpy as np
import os

if True:
    names, data = gt.load_csv(path+'\\runme.gts')
    setup = gt.setup()
    
    def sample():
        for name in names:
            if name in ['Strategy_SourceDirectory_path','Strategy_TargetDirectory_path','Strategy_DiscreteFractures_path',
                        'Strategy_Design_type','Strategy_Target_type','Strategy_Stim_type']:
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
                    val = gt.norm_trunc(1,mu=0.5*(hig+low),
                                        dev=0.25*(hig-low),
                                        lo=low,hi=hig)[0]
                elif dis == 'loguniform':
                    # val = 10.0**(np.random.uniform(np.log10(float(low)),np.log10(float(hig))))
                    val = 10.0**(np.random.uniform()*(np.log10(hig)-np.log10(low))+np.log10(low))
            setattr(setup,name,val)
            # if name in ['Strategy_TargetDirectory_path']:
            #     print('inital working directory %s' %(os.getcwd()))
            #     os.chdir(data[name][1])
            #     print('modified working directory %s' %(os.getcwd()))
    
    counter = 0
    iters = int(float(data['Strategy_Iterations_units'][1]))
    if 'Strategy_Iterations_units' in names:
        iters = int(float(data['Strategy_Iterations_units'][1]))
        for i in range(0,iters):
            print('*******************************************************************')
            print('                    running realization %i' %(i))
            print('*******************************************************************\n\n')
            sample()
            gt.RUN(setup)
    else:
        print('*******************************************************************')
        print('                    running realization %i' %(i))
        print('*******************************************************************\n\n')
        sample()
        gt.RUN(setup)
        
else:
    print('Sorry, runme.gts was not created as it should have been.')

