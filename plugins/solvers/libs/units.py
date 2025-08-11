""" 
Useful unit conversions for GeoDT
"""
#initializations
import numpy as np

#to metric
lpm=(0.1*0.1*0.1)/(60.0)
ft=12*25.4e-3#m
m=(1.0/ft)#ft
deg=1.0*np.pi/180.0
#rad=1.0/deg
gal=3.785*0.1*0.1*0.1#m^3=3.785 l
gpm=gal/60.0
liter=0.1*0.1*0.1
lps=liter/1.0# = liters per second
cP=10.0**-3#Pa-s
g=9.81#m/s2
MPa=10.0**6.0#Pa
GPa=10.0**9.0#Pa
darcy=9.869233*10**-13#m2
mD=darcy*10.0**3.0#m2
yr=365.2425*24.0*60.0*60.0#s
mLmin = 1.66667e-8 #m3/s
mDft = 3.324e15 #m2m

#from metric