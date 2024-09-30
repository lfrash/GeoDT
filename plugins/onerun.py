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

try:
    names, data = gt.load_csv(path+'\\runme.gts')
    setup = gt.setup()
    
    def fetch_musty_old_toy():
        for name in names:
            val = data[name][1]
            try:
                setattr(setup,name,float(val))
            except:
                setattr(setup,name,val)
    fetch_musty_old_toy()
    gt.RUN(setup)
except:
    print('Sorry, runme.gts was not created as it should have been.')

