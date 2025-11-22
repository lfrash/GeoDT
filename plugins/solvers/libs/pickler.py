""" 
pickler
"""
import os
import pickle
if __package__ is None or __package__ == '':
    import iogt
    import wells
    import fractures
    import stress
else:
    from . import iogt
    from . import wells
    from . import fractures
    from . import stress

class databrick:
    def __init__(self):
        pass

def fill(source,target):
    for item in vars(source):
        setattr(target,item,getattr(source,item))

def empty(source,target):
    for item in vars(source):
        if hasattr(target,item):
            setattr(target,item,getattr(source,item))

def packup(obj,pin=100000000):
    try:
        os.makedirs('replay')
    except:
        pass
    brick = databrick()
    fill(obj,brick)
    with open('replay\%i.geodt' %(pin),'wb') as file:
        pickle.dump(brick,file)

#correct for path issues between which module does the pickling and which module unpickles
#  note: check errors and read pickle as text to debug
class repickle(pickle.Unpickler):
    def find_class(self, module, name):
        remodule = module
        #print(module)
        if (module == "solvers.geodt") or (module == "GeoDT"):
            remodule = "plugins.solvers.geodt"
        elif (module == "libs.pickler") or (module == "solvers.libs.pickler"):
            remodule = "plugins.solvers.libs.pickler"
        elif (module == "libs.iogt") or (module == "solvers.libs.iogt"):
            remodule = "plugins.solvers.libs.iogt"
        elif (module == "libs.fractures") or (module == "solvers.libs.fractures"):
            remodule = "plugins.solvers.libs.fractures"
        elif (module == "libs.wells") or (module == "solvers.libs.wells"):
            remodule = "plugins.solvers.libs.wells"
        elif (module == "libs.stress") or (module == "solvers.libs.stress"):
            remodule = "plugins.solvers.libs.stress"
        return super(repickle, self).find_class(remodule, name)

def unpickle(file):
    return repickle(file).load()

def unpack(obj,pin=100000000):
    if True:
        with open('replay\%i.geodt' %(pin),'rb') as file:
            minime = unpickle(file)
            # minime = pickle.load(file)
            for item in vars(minime):
                setattr(obj,item,getattr(minime,item)) #will add if attribute is missing, if this works poorly, consider "hasattr(target,item)"
    else:
        print('pin %s replay file not found' %(pin))