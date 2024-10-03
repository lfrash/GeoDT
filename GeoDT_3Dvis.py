# -*- coding: utf-8 -*-
"""
Created on Sep 29 00:00:00 2024

@author: Luke P Frash
"""

#import solvers
from plugins.solvers import geodt as gt

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

#fetch results from replay file
geom = gt.core()
geom.pin = int(filename[-15:-6])
path = filename[:-23]
print(geom.pin)
gt.unpack(geom,geom.pin)
#export vtk
directory = path #os.getcwd()
try:
    os.makedirs('vtk')
except:
    pass
geom.build_vtk('%s/vtk/%i' %(os.getcwd(),geom.pin),[1,1,1,1,1,1]) #get everything [wells, fracs, fnets, nodes, pipes, bound]
root.destroy()