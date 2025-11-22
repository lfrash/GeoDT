# -*- coding: utf-8 -*-
"""
Convert gts files to set files (row vs column format)
"""

#imports
import numpy as np
# from plugins.solvers.libs import csv
from plugins.solvers.libs import iogt
import tkinter as tk
from tkinter import filedialog
import os

#select settings file
root = tk.Tk()
root.title('GeoDT')
root.geometry('250x50')
print("select a settings file (*.gts)")
filename = filedialog.askopenfilename(initialdir=os.getcwd())
if filename == '':
    print("no file selected")
root.destroy()

# #load settings
# names, data = csv.load_csv(filename)
# out = []
# for i in range(0,len(names)):
#     line = [names[i]]
#     # for j in range(0,len(data)):
#         # if np.isvoid(data[i][j]):
#         #     data[i][j] = ''
#         # line += [data[i][j]]
#     out += [[names[i], data[names[i]][0], data[names[i]][1], data[names[i]][2], data[names[i]][3]]]
# csv.save_set(out,filename[:-4]+'.set')

# csv.convert(filename)

iogt.convert_5_to_6(filename)


