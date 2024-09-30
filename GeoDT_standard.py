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
import sys

#select settings file (*.gts)
root = tk.Tk()
root.title('GeoDT')
root.geometry('250x50')
filename = filedialog.askopenfilename(initialdir=os.getcwd())
if filename == '':
    filename = sys.path[0]+'\\settings\\DEFAULT.gts'
head, data = gt.load_csv(filename)
gt.special_csv(head,data,sys.path[0],os.getcwd(),sys.path[0]+'\\plugins\\runme.gts',append=False)
root.destroy()

#execute GeoDT
import plugins.multirun