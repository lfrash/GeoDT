# -*- coding: utf-8 -*-
"""
Tools for data vizualization and statistical analysis
1. Load data
2. Filter data (direct and relative)
3. 3D scatterplots
4. Timeseries
5. Quantiles
6. Percentiles
Luke P. Frash
"""

#requirements
import numpy as np
import pylab
import matplotlib.pyplot as plt

#2 axis scatter plot with color per 3rd axis
def scatter_3D(x,y,c,name=['x','y','c','save'],cmap='rainbow_r',
               xlog=False,ylog=False,vrange=[],view=True,ax=[]):
    if not(ax):
        fig = pylab.figure(figsize=(10.0,5.0),dpi=100)
        ax1 = fig.add_subplot(111)
    else:
        ax1 = ax
    if vrange:
        for d in range(0,len(c)):
            if c[d] == np.inf:
                c[d] = vrange[1] 
        sc = ax1.scatter(x=x,y=y,c=c,s=8,vmin=vrange[0],vmax=vrange[1],cmap=plt.cm.get_cmap(cmap))
    else:
        sc = ax1.scatter(x=x,y=y,c=c,s=8,cmap=plt.cm.get_cmap(cmap))
    plt.colorbar(ax=ax1,mappable=sc,label=name[2],orientation='vertical')
    ax1.set_xlabel(name[0],fontsize=10)
    ax1.set_ylabel(name[1],fontsize=10)
    if xlog: ax1.set_xscale('log')
    if ylog: ax1.set_yscale('log')
    plt.tight_layout()
    if name[3]:
        plt.savefig(name[3]+'.png', format='png')
    if not(view): 
        pylab.close()

#multi-timeseries
def timeseries(ts,ys,name=['t','y','save'],
               xlog=False,ylog=False,vrange=[],view=False,ax=[]):
    if not(ax):
        fig = pylab.figure(figsize=(10.0,5.0),dpi=100)
        ax1 = fig.add_subplot(111)
    else:
        ax1 = ax
    for s in range(0,len(ys)):
        ax1.plot(ts,ys[s,:])
    ax1.set_xlabel(name[0],fontsize=10)
    ax1.set_ylabel(name[1],fontsize=10)
    plt.tight_layout()
    if name[2]:
        plt.savefig(name[2]+'.png', format='png')
    if not(view): 
        pylab.close()

#percentile
def quantile(x,pcts=[],ax=[]):
    q = np.copy(x)
    q = np.sort(q)
    if len(pcts) < 0:
        pcts = [0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99,100]
    nums = len(q)    
    qs = []
    for p in pcts:
        qs += [q[np.min([nums-1,int(1.0*nums*p/100)])]]
    if ax:
        for i in range(0,len(qs)):
            ax.bar(pcts[1:-1],qs[1:-1])
        ax.set_xlabel('Quantiles',fontsize=10)
        ax.set_ylabel('Value',fontsize=10)
        plt.tight_layout()
    return np.asarray([pcts,qs])

# # Visualization tools
# class visualization_SP:
#     def __init__(self,filename='setup_payoff.csv'):
#         #import data and store unaltered copy
#         try:
#             self.names, self.data = load_csv(filename)
#         except:
#             geom = core()
#             setu = setup()
#             geom.rock.fetch(setu)
#             geom.setup_payoff('setup_payoff_placeholder.csv',geom.pin,aux=[],append=False)
#             self.names, self.data = load_csv('setup_payoff_placeholder.csv')
#         # if np.size(self.data) > 1:
#         self.orig = np.copy(self.data)
#         self.normalize()
#     def normalize(self):
#         #normalize data
#         self.data_n = np.copy(self.data)
#         for n in self.names:
#             try:
#                 np.nan_to_num(self.data_n[n],copy=False,nan=0.0)
#                 #TODO: Add special cases
#                 # if n == 'Well_RotationSkew_rad':
#                 #     self.data_n[n] = np.abs(self.data_n[n])
#                 span = np.max(self.data_n[n])-np.min(self.data_n[n])
#                 if span > 0:
#                     self.data_n[n] = (self.data_n[n]-np.min(self.data_n[n]))/span
#                 else:
#                     self.data_n[n] = 0.0
#             except:
#                 self.data_n[n] = np.zeros(len(self.data[n]))


#     #data filter
#     def multifilter(self,target='pin',compare=None,vs=[0.0,1e12]):
#         keep = np.ones(len(self.data),dtype=bool)
#         try:
#             if compare:
#                 keep = (self.data[target] < vs[1]*self.data[compare]) * (self.data[target] > vs[0]*self.data[compare])
#             else:
#                 keep = (self.data[target] < vs[1]) * (self.data[target] > vs[0])
#         except:
#             print('%s not in setup_payoff' %(target))
#         self.data = self.data[keep]
#         self.data_n = self.data_n[keep]
#     #reset
#     def reset(self):
#         self.data = np.copy(self.orig)
#         self.normalize()

#     #percentiles
#     def percentiles(self,filename='percentiles'):
#         pcts = [0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99,100]
#         first = True
#         for p in pcts:
#             out = [['Quantile','%i' %(p)]]
#             out += [['Percentile','%i' %(100-p)]]
#             for n in self.names:
#                 try:
#                     pv = '%.3e' %(self.quantile(self.data[n],[p])[1][0])
#                 except:
#                     pv = ''
#                 out += [[n, pv]]
#             if first:
#                 save_csv(out,filename,False)
#                 first = False
#             else:
#                 save_csv(out,filename,True)
        