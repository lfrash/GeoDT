# -*- coding: utf-8 -*-
""" 
enumeration of pipe types for cross-referencing
"""
#initializations
import numpy as np

#definitions and cross referencing for pipe types
def typ(key=6):
    ret = []
    choices = np.asarray([
            ['procluster','-6'],
            ['shunt','-5'],
            ['screen','-4'],
            ['perfcluster','-3'],
            ['producer', '-2'],
            ['injector', '-1'],
            ['pipe', '0'],
            ['fracture', '1'],
            ['propped', '2'],
            ['darcy', '3'],
            ['choke', '4'],
            ['perf', '5'],
            ['boundary', '6'] #last choice should always be boundary for [-1] referencing
            ])
    key = str(key)
    ret = np.where(choices == key)
    if ret[1] == 0:
        ret = int(choices[ret[0],1][0])
    elif ret[1] == 1:
        ret = str(choices[ret[0],0][0])
    else:
        print( '**invalid pipe type defined**')
        ret = []
    return ret
