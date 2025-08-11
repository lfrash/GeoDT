# -*- coding: utf-8 -*-
print('prospectivity mapping tool')

# ****************************************************************************
# Calculate geothermal prospectivity from color maps and GeoDT simulations
# Author: Luke P. Frash
#
# Notation:
#  !!! or TODO for code needing to be updated (e.g, limititations or work in progress)
#  *** for breaks betweeen key sections
# ****************************************************************************

#imports
import numpy as np
from plugins.solvers import geodt as gt
from plugins.solvers import gringarten as gg
from plugins.solvers.libs import gt_properties as properties
water = properties.water()
import matplotlib.pyplot as plt
import pylab
import tkinter as tk
from tkinter import filedialog
import os

#load image file
if True:
    #file browser
    if False:
        #select settings file
        root = tk.Tk()
        root.title('GeoDT')
        root.geometry('250x50')
        print("select an image (*.bmp)")
        filename = filedialog.askopenfilename(initialdir=os.getcwd())
        if filename == '':
            print("no file selected")
        root.destroy()
        print(filename)
    
    #heat gradient
    if True:
        filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - geothermal gradient.png'
        # filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - Andersen Stress Coefficient.png'
        # filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - satellite view.png'
        # filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - seismic hazard.png'
        grad = pylab.imread(filename)
        pylab.figure()
        view = pylab.imshow(grad)
        colors = [grad[678,136],grad[696,136],grad[713,136],grad[729,136],grad[747,136],grad[763,136]]
        values = [15,25,35,45,55,65]
        grad = grad[125:920,450:1775,:]
        s = np.shape(grad)
        yx_g = np.zeros((s[0],s[1],1))
        for x in range(0,s[1]):
            for y in range(0,s[0]):
                if grad[y,x,2] > 0.55:
                    yx_g[y,x,0] = 2.95e-1*grad[y,x,0]**2+3.60e-3*grad[y,x,0]+1.43e-1
                else:
                    yx_g[y,x,0] = -1.40e-1*grad[y,x,1]**2-1.22e-1*grad[y,x,1]+6.63e-1
        pylab.figure()
        view = pylab.imshow(yx_g,cmap='Greys_r')
        pylab.savefig('gradient.png')
        np.save('gradient',yx_g)
        
    #stress regime
    if True:
        #stress domain
        # filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - geothermal gradient.png'
        filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - Andersen Stress Coefficient.png'
        # filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - satellite view.png'
        # filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - seismic hazard.png'
        stress = pylab.imread(filename)
        pylab.figure()
        view = pylab.imshow(stress)
        colors = [stress[170,677],stress[226,634],stress[213,623],stress[267,539],stress[239,512],stress[377,137]]
        values = [0.40,           0.55,           0.70,           0.90,           1.75,           3.00]
        stress = stress[125:920,450:1775,:]
        s = np.shape(stress)
        yx_s = np.zeros((s[0],s[1],1))
        for x in range(0,s[1]):
            for y in range(0,s[0]):
                if (stress[y,x,0]+stress[y,x,1]) < 0.70:
                    yx_s[y,x,0] = 3.64*stress[y,x,1]**2-2.08*stress[y,x,1]+0.687
                else:
                    yx_s[y,x,0] = 4.32*stress[y,x,0]**2-2.21*stress[y,x,0]+0.920
        pylab.figure()
        view = pylab.imshow(yx_s,cmap='Greys_r')
        pylab.savefig('stress.png')
        np.save('stress',yx_s)
    
    #basemap
    if False:
        #stress domain
        # filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - geothermal gradient.png'
        # filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - Andersen Stress Coefficient.png'
        filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - satellite view.png'
        # filename = 'C:/Python/scratch/Feb2025/Basin and Range/GeoMap - seismic hazard.png'
        basemap = pylab.imread(filename)
        pylab.figure()
        view = pylab.imshow(basemap)
        basemap = basemap[125:920,450:1775,:]
        pylab.figure()
        view = pylab.imshow(basemap)
        pylab.savefig('basemap.png')
        np.save('basemap',basemap)
    
    #get prospectivity
    grid = np.load('grid.npy')
    gradient = np.load('gradient.npy')*0.1
    stress = np.load('stress.npy')
    basemap = np.load('basemap.npy')
    prospect = np.zeros(np.shape(gradient))
    market = np.zeros(np.shape(gradient))
    capital = np.zeros(np.shape(gradient))
    leadtime = np.zeros(np.shape(gradient))
    powers = np.zeros(np.shape(gradient))
    bounds = [[grid[0,0,0,0],grid[-1,-1,0,0]],[grid[0,0,0,1],grid[-1,-1,0,1]]]
    for x in range(0,np.shape(gradient)[1]):
        for y in range(0,np.shape(gradient)[0]):
            gi = np.min([len(grid)-1,int(len(grid)*(gradient[y,x,0]-bounds[0][0])/(bounds[0][1]-bounds[0][0])+0.49)])
            si = np.min([len(grid)-1,int(len(grid)*(stress[y,x,0]-bounds[1][0])/(bounds[1][1]-bounds[1][0])+0.00)])
            lcoe  = grid[gi,si,1,:]
            capex = grid[gi,si,2,:]
            depth = grid[gi,si,3,:]
            power = grid[gi,si,4,:]
            time  = grid[gi,si,5,:]
            prospect[y,x,:] = np.min([lcoe[2],295]) * capex[2] #prospecitivity is best when capex and lcoe are both low
            market[y,x,:] = np.min([lcoe[1],295]) #also plot the cost of electricity
            capital[y,x,:] = capex[1] #also plot the capital cost
            leadtime[y,x,:] = time[1]/24
            powers[y,x,:] = power[1]*1e-3
            # print(y)
        print(x)
        
    # visu = gt.visualization_SP('gringartenMap0.csv')
    # stuff = (visu.data['Drilling_BottomTemperature_K']-283.15)/visu.data['Well_TargetDepth_m']
    # visu.data['Well_Dip_rad'] = stuff
    # visu.orig['Well_Dip_rad'] = stuff
    # for x in range(0,np.shape(gradient)[1]):
    #     for y in range(0,np.shape(gradient)[0]):
    #         visu.multifilter(target='Well_Dip_rad',vs=[gradient[y,x]*0.9,gradient[y,x]*1.1]) #assume 10% uncertainty
    #         visu.multifilter(target='StressMin_Coefficient_ratio',vs=[stress[y,x]*0.9,stress[y,x]*1.1]) #assume 10% uncertainty
    #         lcoe = visu.quantile(visu.data['LCOE_Net_USDpMWh'],[10,50,90])[1]
    #         capex = visu.quantile(visu.data['Economics_Capital_USD'],[10,50,90])[1]
    #         depth = visu.quantile(visu.data['Well_TargetDepth_m'],[10,50,90])[1]
    #         power = visu.quantile(visu.data['Power_Bulk_kW'],[10,50,90])[1]
    #         time = visu.quantile((visu.data['Drilling_Time_h']+visu.data['Stimulation_Time_h']),[10,50,90])[1]
    #         prospect[y,x,:] = lcoe[2] * capex[2] #prospecitivity is best when capex and lcoe are both low
    #         market[y,x,:] = lcoe[1] #also plot the cost of electricity
    #         capital[y,x,:] = capex[1] #also plot the capital cost
    #         visu.reset()
    #         print(y)
    #     print(x)
    fig = pylab.figure(figsize=(13.0,7.0),dpi=200)
    view = pylab.imshow(basemap)
    pylab.savefig('basemap.png')
    fig = pylab.figure(figsize=(13.0,7.0),dpi=200)
    view = pylab.imshow(gradient*1e3,cmap='Spectral_r')
    fig.colorbar(mappable=view,label='Thermal Gradient (C/km)',orientation='vertical')
    fig = pylab.figure(figsize=(13.0,7.0),dpi=200)
    view = pylab.imshow(stress)
    fig.colorbar(mappable=view,label='Stress (0-1 normal, 1-2 strike-slip, 2-3 reverse)',orientation='vertical')
    fig = pylab.figure(figsize=(13.0,7.0),dpi=200)
    view = pylab.imshow(prospect,cmap='rainbow_r',vmin=np.min(prospect),vmax=10*np.min(prospect))
    fig.colorbar(mappable=view,label='Prospectivity ($USDx$USD/MWh)',orientation='vertical')
    pylab.savefig('prospect.png')
    fig = pylab.figure(figsize=(13.0,7.0),dpi=200)
    view = pylab.imshow(market,cmap='rainbow',vmin=0,vmax=300)
    fig.colorbar(mappable=view,label='Cost of Electricity ($USD/MWh)',orientation='vertical')
    pylab.savefig('market.png')
    fig = pylab.figure(figsize=(13.0,7.0),dpi=200)
    view = pylab.imshow(capital,cmap='rainbow')
    fig.colorbar(mappable=view,label='Capital Cost ($USD)',orientation='vertical')
    pylab.savefig('capital.png')
    fig = pylab.figure(figsize=(13.0,7.0),dpi=200)
    view = pylab.imshow(leadtime,cmap='rainbow')
    fig.colorbar(mappable=view,label='Lead Time (d)',orientation='vertical')
    pylab.savefig('leadtime.png')
    fig = pylab.figure(figsize=(13.0,7.0),dpi=200)
    view = pylab.imshow(powers,cmap='rainbow_r')
    fig.colorbar(mappable=view,label='Power (MW)',orientation='vertical')
    pylab.savefig('powers.png')
    np.save('prospect',prospect)
    np.save('market',market)
    np.save('capital',capital)
    np.save('powers',powers)
    np.save('leadtime',leadtime)
            
            
            
            
            
            
            
            
            
    
    
    # img = img[125:920,450:1775,:] #img[-910:-145,-1800:-465,:] #
    # view = pylab.imshow(img)


pylab.show()