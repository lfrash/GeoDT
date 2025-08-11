# -*- coding: utf-8 -*-
"""
A solver to estimate maximum injection pressure and rate for EGS & CGS
"""
print('circuit_2.0')

#required libraries
import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib.collections import LineCollection

#Plot lines with colored by 3rd value
def multiline(xs, ys, c, ax=None, **kwargs):
    #find axes
    ax = plt.gca() if ax is None else ax
    #create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)
    #set coloring of line segments
    lc.set_array(np.asarray(c))
    #add lines to axes and rescale 
    ax.add_collection(lc)
    ax.autoscale()
    return lc

class cagelimit:
    def __init__(self,w_S=np.asarray([500.0]), w_R=0.1016, RtoS=1.33, LtoS=10.0, 
                 f_rho=980.0, f_mu=0.9e-3, w_f=80.0, p_whp=101e3, f_Pc=2.3e6, 
                 f_s=1, f_n=1, w_n=2, r_E=60.0e9, r_v=0.32, r_rho=2800.0, r_Kn=0.5,
                 f_k=1.0e-10, f_e=3.0e-8, f_b=0.001):
        #constants
        g = 9.81 #m/s2
        #input variables
        self.w_S = np.asarray(w_S) #m, well spacing
        self.w_R = w_R #m, well radius (hydraulic)
        self.RtoS = RtoS #m/m, fracture radius vs spacing
        self.LtoS = LtoS #m/m, well length & depth vs spacing
        self.f_rho = f_rho #kg/m3, fluid density
        self.f_mu = f_mu #Pa-s, water viscosity
        self.w_f = w_f #metric, Hazen-Williams pipe roughness - varies from 150 (smooth) to 80 (rough)
        self.p_whp = p_whp #Pa, production pressure at wellhead
        self.f_Pc = f_Pc #Pa, critical pressure for fracture propagation
        self.f_s = f_s #number of fracture stages
        self.f_n = f_n #number of fracture strands
        self.w_n = w_n #number of boundary wells
        self.r_E = r_E #Pa, Young's modulus
        self.r_v = r_v #Poisson's ratio
        self.r_rho = r_rho #kg/m3, rock density
        self.r_Kn = r_Kn #Pa/Pa, normal stress to vertical stress ratio
        self.f_k = f_k #m2, proppant permeability
        self.f_e = f_e #1/Pa, proppant compressibility
        self.f_b = f_b #m, uncompressed propped fracture aperture
        #intermediate variables
        self.f_R = w_S*RtoS                                        #m, fracture radius
        self.f_w = 8.0*f_Pc*(1.0-r_v**2.0)*self.f_R / (np.pi*r_E)  #m, fracture total aperture
        self.f_h = (f_s*f_n*(self.f_w/f_n)**3.0)**(1.0/3.0)        #m, fracture hydraulic apertue
        self.c_r = 6.0                                             #m/m, ratio of choke radius to separation distance
        self.c_R = w_S/self.c_r                                    #m, choke radius
        self.w_l = LtoS*w_S                                        #m, well length
        self.r_d = LtoS*w_S                                        #m, reservoir depth
        self.bhp = self.w_l*self.f_rho*9.81 #Pa, bottomhole pressure
        self.s3 = self.w_l*r_rho*9.81*r_Kn #Pa, minimum principal stress
        self.p_pro = self.p_whp - self.s3 + self.bhp #Pa, outlet pressure adjusted for hydrostatic effects
        self.f_c = f_k*f_b #m2-m, fracture conductivity
        #output variables
        self.qinj = [] #m3/s, critical injection rate
        self.Pn = [] #Pa, critical net pressure
        self.Ps = [] #Pa, position and pressure of each keypoint
        self.Pq = [] #Pa, pressure at reduced flow rates
        self.Qq = [] #m3/s, reduced flow rates
        #pressure losses
        C1 = (10.7/0.9e-3)*(self.w_l*f_mu*f_rho*g*(1.0/w_n)**1.852)/(w_f**1.852*(2*w_R)**4.87) 
        C2 = 6.0*f_mu*(0.5*self.c_R**2*np.log(self.c_R/w_R) - 0.25*self.c_R**2.0)/(np.pi*self.f_h**3.0)
        C3 = 6.0*f_mu*(1.0+1.0/w_n)*(((1.0-2.0/self.c_r)*w_S)**2.0)/(2.0*self.f_h**3.0) 
        C4 = 6.0*f_mu*np.log(self.c_R/w_R)/(np.pi*self.f_h**3.0) 
        C5 = 6.0*f_mu*(1.0+1.0/w_n)*(1.0-2.0/self.c_r)/self.f_h**3.0 
        C6 = 2.0*np.pi*self.c_R**2.0 + (1.0-2.0/self.c_r)*w_S**2.0
        #first estimate of max cageable flow by quadratic equation
        A = C1*(2.0*np.pi*self.c_R**2.0 + (1.0-2.0/self.c_r)*w_S**2.0)
        B = 2.0*C2 + C3 + C4*(np.pi*self.c_R**2.0 + (1.0-2.0/self.c_r)*w_S**2.0) + C5*(np.pi*self.c_R**2.0)
        C = C6*(self.p_pro - f_Pc)
        D = B**2.0 - 4.0*A*C
        self.qinj = np.max(np.asarray([(-B+D**0.5)/(2.0*A),(-B-D**0.5)/(2.0*A)]),axis=0)
        #critical net pressure binary search
        qs = np.zeros((3,len(self.qinj)),dtype=float)
        qs[0,:] = self.qinj*0.1
        qs[1,:] = self.qinj*1.0
        qs[2,:] = self.qinj*1.9
        for j in range(0,20):
            ck = ((C2*qs[1,:] + (C1*qs[1,:]**1.852 + self.p_pro + C4*qs[1,:] + C5*qs[1,:])*np.pi*self.c_R**2.0 + 
                      C3*qs[1,:] + (C1*qs[1,:]**1.852 + self.p_pro + C4*qs[1,:])*(1.0-2.0/self.c_r)*self.w_S**2.0 +
                      C2*qs[1,:] + (C1*qs[1,:]**1.852 + self.p_pro)*np.pi*self.c_R**2.0)/C6)
            ck[ck <= f_Pc] = -1
            ck[ck > f_Pc] = 0
            qs[0,:] = qs[0,:]*(1+ck[:]) + qs[1,:]*(-ck[:])
            qs[2,:] = qs[1,:]*(1+ck[:]) + qs[2,:]*(-ck[:])
            qs[1,:] = 0.5*(qs[0,:]+qs[2,:])
        self.qinj = qs[1,:]
        self.Pn = ((C2*self.qinj + (C1*self.qinj**1.852 + self.p_pro + C4*self.qinj + C5*self.qinj)*np.pi*self.c_R**2.0 + 
                C3*self.qinj + (C1*self.qinj**1.852 + self.p_pro + C4*self.qinj)*(1.0-2.0/self.c_r)*self.w_S**2.0 +
                C2*self.qinj + (C1*self.qinj**1.852 + self.p_pro)*np.pi*self.c_R**2.0)/C6)
        #frictional losses
        def losses(qinj):
            #pressure drops
            dPi = 6.0*qinj*f_mu*np.log(self.c_R/self.w_R)/(np.pi*self.f_h**3.0)
            dPf = 6.0*qinj*f_mu*(1+1/self.w_n)*(1.0-2.0/self.c_r)/self.f_h**3.0
            dPp = dPi/self.w_n
            dPwi = (10.7/0.9e-3)*(self.w_l*f_mu*f_rho*g*qinj**1.852)/(self.w_f**1.852*(2*self.w_R)**4.87)
            dPwp = (10.7/0.9e-3)*(self.w_l*f_mu*f_rho*g*(qinj/self.w_n)**1.852)/(self.w_f**1.852*(2*self.w_R)**4.87) #
            #spacing must be greater than well diameter
            if np.size(self.c_R) > 1:
                invalid = np.zeros(len(self.c_R[self.c_R<2.0*self.w_R]),dtype=float)
                invalid[:] = np.nan
                self.Pn[0:len(invalid)] = invalid
                qinj[0:len(invalid)] = invalid
                dPi[0:len(invalid)] = invalid
                dPf[0:len(invalid)] = invalid
                dPp[0:len(invalid)] = invalid
                dPwi[0:len(invalid)] = invalid
                dPwp[0:len(invalid)] = invalid
            #head pressures
            P6 = self.p_pro*np.ones(len(qinj)) #production surface
            P5 = P6 + dPwp #fracture exit entry
            P4 = self.Pn - 0.5*dPf #fracture net
            P3 = P4 + dPf #near well
            P2 = P3 + dPi #fracture entry
            P1 = P2 + dPwi #injection wellhead
            head = np.asarray([P1,P2,P3,P4,P5,P6])
            #absolute pressures
            real = head + self.s3
            real[0] = real[0] - self.bhp
            real[5] = real[5] - self.bhp
            return head, real
        #critical pressure profile
        a1 = np.ones(len(self.qinj))
        self.Ps = np.zeros((3,6,len(self.qinj)))
        self.Ps[0,:,:] = np.asarray([-self.w_l*a1,
                                self.w_R*a1,
                                (1.0/self.c_r)*self.w_S*a1,
                                (1.0 - 1.0/self.c_r)*self.w_S*a1,
                                self.w_S*a1-self.w_R,
                                self.w_S*a1+self.w_l])
        head, real = losses(self.qinj)
        self.Ps[1,:,:] = head
        self.Ps[2,:,:] = real
        #compressible propped fracture flow
        def propped(qinj):
            #setup binary search
            pn_guess = np.zeros(3)
            pn_guess[0] = self.p_pro #maximum compressive stress (compression negative)
            pn_guess[2] = self.Pn #maximum positive stress (tensile opening)
            pn_guess[1] = 0.5*(pn_guess[2]+pn_guess[0]) #binary search first guess
            pn_calc = 0.0
            for k in range(0,22):
                #calculate hydraulic aperture
                self.f_p = (self.f_s * 12.0 * self.f_k * self.f_b*np.exp(self.f_e*pn_guess[1]))**(1.0/3.0)
                #pressure drops
                dPi = 6.0*qinj*f_mu*np.log(self.c_R/self.w_R)/(np.pi*self.f_p**3.0)
                dPf = 6.0*qinj*f_mu*(1+1/self.w_n)*(1.0-2.0/self.c_r)/self.f_p**3.0
                dPp = dPi/self.w_n
                dPwi = (10.7/0.9e-3)*(self.w_l*f_mu*f_rho*g*qinj**1.852)/(self.w_f**1.852*(2*self.w_R)**4.87)
                dPwp = (10.7/0.9e-3)*(self.w_l*f_mu*f_rho*g*(qinj/self.w_n)**1.852)/(self.w_f**1.852*(2*self.w_R)**4.87) #
                #spacing must be greater than well diameter
                if np.size(self.c_R) > 1:
                    invalid = np.zeros(len(self.c_R[self.c_R<2.0*self.w_R]),dtype=float)
                    invalid[:] = np.nan
                    self.Pn[0:len(invalid)] = invalid
                    qinj[0:len(invalid)] = invalid
                    dPi[0:len(invalid)] = invalid
                    dPf[0:len(invalid)] = invalid
                    dPp[0:len(invalid)] = invalid
                    dPwi[0:len(invalid)] = invalid
                    dPwp[0:len(invalid)] = invalid
                #head pressures
                P6 = self.p_pro*np.ones(len(qinj)) #production surface
                P5 = P6 + dPwp #fracture exit entry
                P4 = P5 + dPi/w_n #fracture net
                P3 = P4 + dPf #near well
                pn_calc = 0.5*(P4+P3) #net pressure
                P2 = P3 + dPi #fracture entry
                P1 = P2 + dPwi #injection wellhead
                head = np.asarray([P1,P2,P3,P4,P5,P6])
                #next guess
                if pn_calc > pn_guess[1]:
                    pn_guess[0] = pn_guess[1]
                else:
                    pn_guess[2] = pn_guess[1]
                pn_guess[1] = 0.5*(pn_guess[2]+pn_guess[0])
            #absolute pressures
            real = head + self.s3
            real[0] = real[0] - self.bhp
            real[5] = real[5] - self.bhp
            return head, real
        # #propped fracture flow
        # div = 100
        # Pi = 10.0**np.linspace(np.log10(0.0001),np.log10(1),div)
        # self.Hq = np.zeros((div,6,len(self.qinj)))
        # self.Pq = np.zeros((div,6,len(self.qinj)))
        # self.Qq = np.zeros((div,1,len(self.qinj)))
        # for i in range(0,len(Pi)):
        #     head, real = propped(self.qinj*Pi[i])
        #     self.Hq[i,:,:] = head + self.s3 - self.bhp
        #     self.Pq[i,:,:] = real
        #     self.Qq[i,0,:] = self.qinj*Pi[i]       
        # #hydropropped flow (fracture aperture assumed constant at critical opening)
        # div = 100
        # Pi = 10.0**np.linspace(np.log10(0.0001),np.log10(1),div)
        # self.Hq = np.zeros((div,6,len(self.qinj)))
        # self.Pq = np.zeros((div,6,len(self.qinj)))
        # self.Qq = np.zeros((div,1,len(self.qinj)))
        # for i in range(0,len(Pi)):
        #     head, real = losses(self.qinj*Pi[i])
        #     self.Hq[i,:,:] = head + self.s3 - self.bhp
        #     self.Pq[i,:,:] = real
        #     self.Qq[i,0,:] = self.qinj*Pi[i]
        
        #subcritical fracture flow
        div = 100
        Pi = 10.0**np.linspace(np.log10(0.001),np.log10(1),div)
        self.Hq = np.zeros((div,6,len(self.qinj)))
        self.Pq = np.zeros((div,6,len(self.qinj)))
        self.Qq = np.zeros((div,1,len(self.qinj)))
        self.Fk = np.zeros((div,1,len(self.qinj)))
        self.Fh = np.zeros((div,1,len(self.qinj)))
        self.Qopen = np.zeros(len(self.qinj))
        for i in range(0,len(Pi)):
            head_p, real_p = propped(self.qinj*Pi[i])
            head_h, real_h = losses(self.qinj*Pi[i])
            for j in range(0,len(self.qinj)):
                if head_p[0,j] < head_h[0,j]:
                    self.Hq[i,:,:] = head_p + self.s3 - self.bhp
                    self.Pq[i,:,:] = real_p
                    self.Qq[i,0,:] = self.qinj*Pi[i]
                    self.Fk[i,0,:] = self.f_p**2.0/12.0
                    self.Fh[i,0,:] = self.f_p
                    self.Qopen[:] = self.qinj*Pi[i]
                else:
                    self.Hq[i,:,:] = head_h + self.s3 - self.bhp
                    self.Pq[i,:,:] = real_h
                    self.Qq[i,0,:] = self.qinj*Pi[i]
                    self.Fk[i,0,:] = self.f_h**2.0/12.0
                    self.Fh[i,0,:] = self.f_h
                
        
        #psuedo step rate test
        div = 6
        Pi = 10.0**np.linspace(np.log10(0.003),np.log10(0.3),div)
        self.HT = np.zeros((div,6,len(self.qinj)))
        self.PT = np.zeros((div,6,len(self.qinj)))
        self.QT = np.zeros((div,1,len(self.qinj)))
        self.Fk = np.zeros((div,1,len(self.qinj)))
        self.Fh = np.zeros((div,1,len(self.qinj)))
        for i in range(0,len(Pi)):
            head_p, real_p = propped(self.qinj*Pi[i])
            head_h, real_h = losses(self.qinj*Pi[i])
            for j in range(0,len(self.qinj)):
                if head_p[0,j] < head_h[0,j]:
                    self.HT[i,:,:] = head_p + self.s3 - self.bhp
                    self.PT[i,:,:] = real_p
                    self.QT[i,0,:] = self.qinj*Pi[i]
                    self.Fk[i,0,:] = self.f_p**2.0/12.0
                    self.Fh[i,0,:] = self.f_p
                else:
                    self.HT[i,:,:] = head_h + self.s3 - self.bhp
                    self.PT[i,:,:] = real_h
                    self.QT[i,0,:] = self.qinj*Pi[i]
                    self.Fk[i,0,:] = self.f_h**2.0/12.0
                    self.Fh[i,0,:] = self.f_h
        
    def p_vs_q(self):
        pass
            
#**********************************************************************
### Frash, 2025, Stanford Geothermal Workshop
#**********************************************************************
if True:
    #setup and run model    
    p_whp = 5.0e6
    w_R = np.asarray([0.5*8*0.0254])
    w_S = 250.0; f_Pc=0.1e6; f_mu=0.20e-3
    RtoS=1.10; LtoS=20.0; f_rho=958.0; w_f=80.0; f_s=50; f_n=2; w_n=2; r_E=60.0e9; r_v=0.30; r_rho=2700.0; r_Kn=0.5
    f_k=1.0e-11; f_e=3.0e-8; f_b=0.001
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    print(m.f_c)

    #plot injection pressure versus flow rate
    dpi = 100
    figsize = (10,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,5,0]*1e-6,color='red',label='production wellhead')
    ax1.plot(m.Qq[:,0,0],m.Hq[:,4,0]*1e-6,color='pink',label='production friction')
    ax1.plot(m.Qq[:,0,0],(0.5*m.Hq[:,2,0]+0.5*m.Hq[:,3,0])*1e-6,color='blue',label='fracture net')
    ax1.plot(m.Qq[:,0,0],m.Hq[:,1,0]*1e-6,color='green',label='injection friction')
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='black',label='injection wellhead')
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.00014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='steelblue',label='fracture opening')
    # xs = [0.005,1.0]
    # ys1 = [5.0,5.0]
    # ys2 = [105.0,105.0]
    # ax1.fill_between(xs,ys1,ys2,color='black',alpha=0.3,label='feasibility range')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='k',label='caging limit')
    ax1.set_ylabel('Pressure Head (MPa)',fontsize=10)
    ax1.set_xlabel('Flow Rate (m3/s)',fontsize=10)
    ax1.set_xscale('log')
    plt.title('Subcritical injection at %.0f m depth' %(m.r_d))
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()
    
    #plot step rate test
    dpi = 100
    figsize = (10,12)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(511)
    t = np.linspace(0,len(m.QT)-0.001,2*len(m.QT),dtype=int)
    t = np.asarray(list(t) + list(len(m.QT)-t-1))
    ts = np.linspace(0,12,len(t))
    p0 = []
    p1 = []
    p4 = []
    p5 = []
    qs = []
    hs = [] 
    for i in range(0,len(ts)):
        p0 += [m.HT[t[i],0,0]*1e-6]
        p1 += [m.HT[t[i],1,0]*1e-6]
        p4 += [m.HT[t[i],4,0]*1e-6]
        p5 += [m.HT[t[i],5,0]*1e-6]
        qs += [m.QT[t[i],0,0]]
        hs += [m.Fh[t[i],0,0]]
    p0 = np.asarray(p0)
    p1 = np.asarray(p1)
    p4 = np.asarray(p4)
    p5 = np.asarray(p5)
    qs = np.asarray(qs)
    plt.title('Step rate injection at %.0f m depth' %(m.r_d))
    ax1.plot(ts,p5,color='red',label='production wellhead')
    ax1.plot(ts,p4,color='pink',label='production bottomhole')
    ax1.plot(ts,p1,color='green',label='injection bottomhole')
    ax1.plot(ts,p0,color='black',label='injection wellhead')
    ax1.plot([ts[0],ts[-1]],[(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6],':',color='steelblue',label='fracture opening')
    ax1.set_ylabel('Pressure Head\n(MPa)',fontsize=10)
    ax1.set_xlabel('Time (h)',fontsize=10)
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)
    ax2 = fig.add_subplot(512)
    ax2.plot(ts,qs,color='blue',label='injection wellhead')
    ax2.plot([ts[0],ts[-1]],[m.Qopen[0],m.Qopen[0]],':',color='steelblue',label='fracture opening')
    ax2.plot([ts[0],ts[-1]],[m.qinj[0],m.qinj[0]],'--',color='blueviolet',label='caging limit')
    ax2.set_ylabel('Injection Rate\n(m3/s)',fontsize=10)
    ax2.set_xlabel('Time (h)',fontsize=10)
    ax2.set_yscale('log')
    ax2.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)
    ax3 = fig.add_subplot(513)
    ax3.plot(ts,qs/(p0-m.p_whp*1e-6),color='orange',label='inferred injectivity')
    ax3.set_ylabel('Injectivity\n(m3/s/MPa)',fontsize=10)
    ax3.set_xlabel('Time (h)',fontsize=10)
    ax3.set_yscale('log')
    ax3.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)
    ax4 = fig.add_subplot(514)
    # ax4.plot(ts,m.f_mu*(qs/(p0*1e6-m.p_whp))/(2.0*np.pi*0.1*m.r_d),color='darkorange',label='inferred permeability')
    ax4.plot(ts,m.f_mu*(qs/(p0*1e6-m.p_whp))/(2.0*np.pi*m.f_b*m.f_s),color='darkorange',label='inferred from injectivity')
    ax4.plot(ts,qs*m.f_mu/(m.f_b*m.f_s*(p1-p4)*1e6),color='darkgreen',label='corrected from pressure and flow')
    ax4.plot(ts,m.f_k*qs/qs,color='darkblue',label='true model input')
    ax4.set_ylabel('Permeability\n(m2)',fontsize=10)
    ax4.set_xlabel('Time (h)',fontsize=10)
    ax4.set_yscale('log')
    ax4.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)
    ax5 = fig.add_subplot(515)
    ax5.plot(ts,hs,color='chocolate',label='hydraulic aperture')
    ax5.set_ylabel('Hydraulic Aperture\n(m)',fontsize=10)
    ax5.set_xlabel('Time (h)',fontsize=10)
    ax5.set_yscale('log')
    ax5.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()    

    #variants - low perm small well deeper - high perm large well shallower
    dpi = 100
    figsize = (10,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    w_R = np.asarray([0.5*6*0.0254]); f_k=1.0e-12; LtoS=28.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='purple',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='orchid',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='purple',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='black',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='steelblue',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='black',label='caging limit')
    w_R = np.asarray([0.5*10*0.0254]); f_k=1.0e-10; LtoS=12.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='darkgreen',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='mediumseagreen',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='darkgreen',label='caging limit')
    ax1.set_ylabel('Pressure Head (MPa)',fontsize=10)
    ax1.set_xlabel('Flow Rate (m3/s)',fontsize=10)
    ax1.set_xscale('log')
    plt.title('Effect of key parameters')
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()

    #variants - well size
    dpi = 100
    figsize = (10,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    w_R = np.asarray([0.5*6*0.0254]); f_k=1.0e-11; LtoS=20.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='purple',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='orchid',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='purple',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='black',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='steelblue',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='black',label='caging limit')
    w_R = np.asarray([0.5*10*0.0254]); f_k=1.0e-11; LtoS=20.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='darkgreen',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='mediumseagreen',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='darkgreen',label='caging limit')
    ax1.set_ylabel('Pressure Head (MPa)',fontsize=10)
    ax1.set_xlabel('Flow Rate (m3/s)',fontsize=10)
    ax1.set_xscale('log')
    plt.title('Effect of key parameters')
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()

    #variants - permeability
    dpi = 100
    figsize = (10,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-12; LtoS=20.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='purple',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='orchid',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='purple',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='black',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='steelblue',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='black',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-10; LtoS=20.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='darkgreen',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='mediumseagreen',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='darkgreen',label='caging limit')
    ax1.set_ylabel('Pressure Head (MPa)',fontsize=10)
    ax1.set_xlabel('Flow Rate (m3/s)',fontsize=10)
    ax1.set_xscale('log')
    plt.title('Effect of key parameters')
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()
    
    #variants - depth
    dpi = 100
    figsize = (10,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=28.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='purple',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='orchid',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='purple',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='black',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='steelblue',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='black',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=12.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='darkgreen',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='mediumseagreen',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='darkgreen',label='caging limit')
    ax1.set_ylabel('Pressure Head (MPa)',fontsize=10)
    ax1.set_xlabel('Flow Rate (m3/s)',fontsize=10)
    ax1.set_xscale('log')
    plt.title('Effect of key parameters')
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()

    #variants - stages
    dpi = 100
    figsize = (10,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0; f_s=25
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='purple',label='%.1f" well, %.1e D, %.0f m, %i fractures' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d,f_s))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='orchid',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='purple',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0; f_s=50
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='black',label='%.1f" well, %.1e D, %.0f m, %i fractures' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d,f_s))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='steelblue',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='black',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0; f_s=100
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='darkgreen',label='%.1f" well, %.1e D, %.0f m, %i fractures' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d,f_s))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='mediumseagreen',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='darkgreen',label='caging limit')
    ax1.set_ylabel('Pressure Head (MPa)',fontsize=10)
    ax1.set_xlabel('Flow Rate (m3/s)',fontsize=10)
    ax1.set_xscale('log')
    plt.title('Effect of key parameters')
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()

    #variants - stranding
    dpi = 100
    figsize = (10,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0; f_s=50; f_n=1
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='purple',label='%.1f" well, %.1e D, %.0f m, %i fractures, %i strands' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d,f_s,f_n))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='orchid',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='purple',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0; f_s=50; f_n=2
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='black',label='%.1f" well, %.1e D, %.0f m, %i fractures, %i strands' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d,f_s,f_n))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='steelblue',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='black',label='caging limit')
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0; f_s=50; f_n=20
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Qq[:,0,0],m.Hq[:,0,0]*1e-6,color='darkgreen',label='%.1f" well, %.1e D, %.0f m, %i fractures, %i strands' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d,f_s,f_n))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(ys,xs,':',color='mediumseagreen',label='fracture opening')
    ax1.plot(m.qinj[0],m.Ps[2,0,0]*1e-6,marker='*',linestyle='',color='darkgreen',label='caging limit')
    ax1.set_ylabel('Pressure Head (MPa)',fontsize=10)
    ax1.set_xlabel('Flow Rate (m3/s)',fontsize=10)
    ax1.set_xscale('log')
    plt.title('Effect of key parameters')
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()

    #variants - flow vs pressure
    dpi = 100
    figsize = (10,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    w_R = np.asarray([0.5*8*0.0254]); f_k=1.0e-11; LtoS=20.0; f_s=50; f_n=2
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_s,f_n,w_n,r_E,r_v,r_rho,r_Kn,f_k,f_e,f_b)
    ax1.plot(m.Hq[:,0,0]*1e-6,m.Qq[:,0,0],color='black',label='%.1f" well, %.1e D, %.0f m' %(2*w_R/0.0254,f_k*1.0132e12,m.r_d))
    xs = [(m.s3-m.bhp)*1e-6,(m.s3-m.bhp)*1e-6]
    ys = [0.0014,np.nanmax(m.Qq[:,0,0])]
    ax1.plot(xs,ys,':',color='steelblue',label='fracture opening')
    ax1.plot(m.Ps[2,0,0]*1e-6,m.qinj[0],marker='*',linestyle='',color='black',label='caging limit')
    ax1.set_xlabel('Pressure Head (MPa)',fontsize=10)
    ax1.set_ylabel('Flow Rate (m3/s)',fontsize=10)
    ax1.set_xscale('log')
    plt.title('Effect of key parameters')
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()
    
    


#**********************************************************************
###Frash et al., 2025, Renewable Energy Journal function
#**********************************************************************
if False:
    #general settings
    dpi = 150
    figsize = (5,5)
    w_S = 10.0**np.linspace(np.log10(0.010),np.log10(1000.0),1000,dtype=float) #m, well spacing
    RtoS=1.05; LtoS=10.0; f_rho=980.0; f_mu=0.9e-3; w_f=80.0; f_mu=0.9e-3; r_rho=2800.0; r_Kn=0.5
    #material
    #r_E=2.55e9; r_v=0.402; p_whp=101e3; f_Pc=2.3e6; w_R=0.00076 #acrylic
    r_E=50.0e9; r_v=0.30; p_whp=1e6; f_Pc=1.0e6; w_R=0.5*10*0.0254 #granite 
    #well and fracture details
    f_n=1; w_n=3 #standard
    #f_n=3; w_n=4 #challenging
    
    #**********************************************************************
    #base model
    #**********************************************************************
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    
    #solution check
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    ax1.plot(m.w_S,m.Pn,color='k',label='Nominal Estimate')
    xs = [np.min(m.w_S),np.max(m.w_S)]
    ys = [m.f_Pc,m.f_Pc]
    ax1.plot(xs,ys,'--',color='grey',label='f_Pc')
    ax1.set_xlabel('Well Spacing (m)',fontsize=10)
    ax1.set_ylabel('Fracture Net Pressure (kPa)',fontsize=10)
    ax1.set_xscale('log')
    ax1.set_ylim([0.0,5.0e6])
    ax1.legend(loc='lower right', prop={'size':10}, ncol=2, numpoints=1)   
    plt.tight_layout()

    #pressure multi-profiles
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    div = 50
    rs = np.linspace(0,1,div,dtype=float)
    xs = np.zeros((div,len(m.Ps[0])),dtype=float)
    ys = np.zeros((div,len(m.Ps[0])),dtype=float)
    cs = np.zeros(div,dtype=float)
    for r in range(0,div):
        i = int(rs[r]*(len(m.qinj)-1))
        xs[r] = 100.0*m.Ps[0,:,i]/np.max(m.Ps[0,:,i])
        ys[r] = m.Ps[1,:,i]
        cs[r] = np.log10(m.w_S[i])
    lc = multiline(xs,ys*10**-6,cs,ax1,cmap='cool')
    axcb = fig.colorbar(lc,label='Log well spacing (m)')
    ax1.set_xlabel('Percent distance along flowline (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.tight_layout()
    
    #**********************************************************************
    #labscale
    #**********************************************************************
    #example single profile (labscale - granite & water)
    w_R=0.00076; f_Pc=2.3e6; p_whp=0.500e6; f_mu=0.9e-3
    RtoS=1.33; LtoS=35.0; f_rho=980.0; w_f=80.0; f_n=1; w_n=4; r_E=50.0e9; r_v=0.30; r_rho=2800.0; r_Kn=0.5
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    i = np.where(m.w_S>0.0381)[0][0] #int(0.30*nel)
    xs = m.Ps[0,:,i]
    ys = m.Ps[1,:,i]
    ax1.plot(xs,ys*10**-6,color='grey',label='pressure head')
    ax1.plot(xs,ys*10**-6,'.',color='grey')
    xs = m.Ps[0,:,i]
    ys = m.Ps[2,:,i]
    ax1.plot(xs,ys*10**-6,color='black',label='absolute pressure')
    ax1.plot(xs,ys*10**-6,'.',color='k')
    ax1.set_xlabel('Distance along flowline (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.title('Labscale granite with water %.1e mL/min' %(m.qinj[i]*6e7))
    ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()

    #example single profile (labscale - granite & oil)
    w_R=0.00076; f_Pc=2.3e6; p_whp=0.500e6; f_mu=404.0e-3
    RtoS=1.33; LtoS=35.0; f_rho=980.0; w_f=80.0; f_n=1; w_n=4; r_E=50.0e9; r_v=0.30; r_rho=2800.0; r_Kn=0.5
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    i = np.where(m.w_S>0.0381)[0][0] #int(0.30*nel)
    xs = m.Ps[0,:,i]
    ys = m.Ps[1,:,i]
    ax1.plot(xs,ys*10**-6,color='grey',label='pressure head')
    ax1.plot(xs,ys*10**-6,'.',color='grey')
    xs = m.Ps[0,:,i]
    ys = m.Ps[2,:,i]
    ax1.plot(xs,ys*10**-6,color='black',label='absolute pressure')
    ax1.plot(xs,ys*10**-6,'.',color='k')
    ax1.set_xlabel('Distance along flowline (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.title('Labscale granite with oil %.1e mL/min' %(m.qinj[i]*6e7))
    ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()
    
    #example single profile (labscale - acrylic & water)
    w_R=0.00076; f_Pc=2.3e6; p_whp=0.500e6; f_mu=0.9e-3
    RtoS=1.05; LtoS=35.0; f_rho=980.0; w_f=80.0; f_n=1; w_n=2; r_E=3.0e9; r_v=0.37; r_rho=2800.0; r_Kn=0.5
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    i = np.where(m.w_S>0.0381)[0][0] #int(0.30*nel)
    xs = m.Ps[0,:,i]
    ys = m.Ps[1,:,i]
    ax1.plot(xs,ys*10**-6,color='grey',label='pressure head')
    ax1.plot(xs,ys*10**-6,'.',color='grey')
    xs = m.Ps[0,:,i]
    ys = m.Ps[2,:,i]
    ax1.plot(xs,ys*10**-6,color='black',label='absolute pressure')
    ax1.plot(xs,ys*10**-6,'.',color='k')
    ax1.set_xlabel('Distance along flowline (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.title('Labscale acrylic with water %.1e mL/min' %(m.qinj[i]*6e7))
    ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()
    
    #example single profile (labscale - acrylic & viscous oil)
    w_R=0.00076; f_Pc=2.3e6; p_whp=0.500e6; f_mu=404.0e-3
    RtoS=1.00; LtoS=35.0; f_rho=980.0; w_f=80.0; f_n=1; w_n=2; r_E=3.0e9; r_v=0.37; r_rho=2800.0; r_Kn=0.5
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    i = np.where(m.w_S>0.0381)[0][0] #int(0.30*nel)
    xs = m.Ps[0,:,i]
    ys = m.Ps[1,:,i]
    ax1.plot(xs,ys*10**-6,color='grey',label='pressure head')
    ax1.plot(xs,ys*10**-6,'.',color='grey')
    xs = m.Ps[0,:,i]
    ys = m.Ps[2,:,i]
    ax1.plot(xs,ys*10**-6,color='black',label='absolute pressure')
    ax1.plot(xs,ys*10**-6,'.',color='k')
    ax1.set_xlabel('Distance along flowline (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.title('Labscale acrylic with oil %.1e mL/min' %(m.qinj[i]*6e7))
    ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()
    
    #+ example single profile (labscale - acrylic & viscous oil)
    w_R=0.00076; f_Pc=2.3e6; p_whp=0.500e6; f_mu=404.0e-3
    RtoS=1.3; LtoS=35.0; f_rho=980.0; w_f=80.0; f_n=1; w_n=2; r_E=3.0e9; r_v=0.37; r_rho=2800.0; r_Kn=0.5
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    i = np.where(m.w_S>0.0381)[0][0] #int(0.30*nel)
    xs = m.Ps[0,:,i]
    ys = m.Ps[1,:,i]
    ax1.plot(xs,ys*10**-6,color='grey',label='pressure head')
    ax1.plot(xs,ys*10**-6,'.',color='grey')
    xs = m.Ps[0,:,i]
    ys = m.Ps[2,:,i]
    ax1.plot(xs,ys*10**-6,color='black',label='absolute pressure')
    ax1.plot(xs,ys*10**-6,'.',color='k')
    ax1.set_xlabel('Distance along flowline (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.title('Labscale acrylic with oil %.1e mL/min x1.3s' %(m.qinj[i]*6e7))
    ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()
    
    #++ example single profile (labscale - acrylic & viscous oil)
    w_R=0.00076; f_Pc=2.3e6; p_whp=0.500e6; f_mu=404.0e-3
    RtoS=2.0; LtoS=35.0; f_rho=980.0; w_f=80.0; f_n=1; w_n=2; r_E=3.0e9; r_v=0.37; r_rho=2800.0; r_Kn=0.5
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    i = np.where(m.w_S>0.0381)[0][0] #int(0.30*nel)
    xs = m.Ps[0,:,i]
    ys = m.Ps[1,:,i]
    ax1.plot(xs,ys*10**-6,color='grey',label='pressure head')
    ax1.plot(xs,ys*10**-6,'.',color='grey')
    xs = m.Ps[0,:,i]
    ys = m.Ps[2,:,i]
    ax1.plot(xs,ys*10**-6,color='black',label='absolute pressure')
    ax1.plot(xs,ys*10**-6,'.',color='k')
    ax1.set_xlabel('Distance along flowline (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.title('Labscale acrylic with oil %.1e mL/min x2.0s' %(m.qinj[i]*6e7))
    ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()
    
    #=== example single profile (labscale - acrylic & viscous oil)
    w_R=0.00076; f_Pc=2.3e6; p_whp=0.500e6; f_mu=404.0e-3
    RtoS=4.81/3.81; LtoS=35.0; f_rho=980.0; w_f=80.0; f_n=1; w_n=2; r_E=3.0e9; r_v=0.37; r_rho=2800.0; r_Kn=0.5
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    i = np.where(m.w_S>0.0381)[0][0] #int(0.30*nel)
    xs = m.Ps[0,:,i]
    ys = m.Ps[1,:,i]
    ax1.plot(xs,ys*10**-6,color='grey',label='pressure head')
    ax1.plot(xs,ys*10**-6,'.',color='grey')
    xs = m.Ps[0,:,i]
    ys = m.Ps[2,:,i]
    ax1.plot(xs,ys*10**-6,color='black',label='absolute pressure')
    ax1.plot(xs,ys*10**-6,'.',color='k')
    ax1.set_xlabel('Distance along flowline (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.title('Labscale acrylic with oil %.1e mL/min 4.78' %(m.qinj[i]*6e7))
    ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()
    
    #**********************************************************************
    #labscale vs fracture radius
    #**********************************************************************
    #acrylic
    figsize = (4,4)
    w_R=0.00076; f_Pc=2.3e6; p_whp=0.500e6; f_mu=404.0e-3
    RtoS=1.05; LtoS=35.0; f_rho=980.0; w_f=80.0; f_n=1; w_n=2; r_E=3.0e9; r_v=0.37; r_rho=2800.0; r_Kn=0.5
    #critical flow rate plot as function of fracture radius
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    # plt.axvline(0.0381,color='grey')
    xs = [0.0381,0.0381]
    ys = [0.5,16.0]
    ax1.plot(xs,ys,'-x',color='grey',label='acrylic test')
    RtoS = 1.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj*6e7,':',color='r',label='Rf = 1.0*S')
    RtoS = 1.3
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj*6e7,color='k',label='Rf = 1.3*S')
    RtoS = 2.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj*6e7,'--',color='b',label='Rf = 2.0*S')
    ax1.set_xlabel('Well Spacing (m)',fontsize=10)
    ax1.set_ylabel('Injection Rate (mL/min)',fontsize=10)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    plt.title('Acrylic with oil and varied fracture radius')
    ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()
    figsize = (5,5)
    
    #granite
    w_R=0.00076; f_Pc=2.3e6; p_whp=0.500e6; f_mu=0.9e-3
    RtoS=1.05; LtoS=35.0; f_rho=980.0; w_f=80.0; f_n=1; w_n=4; r_E=50.0e9; r_v=0.30; r_rho=2800.0; r_Kn=0.5
    #critical flow rate plot as function of fracture radius
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    # plt.axvline(0.0381,color='grey')
    xs = [0.0381,0.0381]
    ys = [0.5,16.0]
    ax1.plot(xs,ys,'-x',color='grey',label='acrylic test')
    RtoS = 1.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj*6e7,':',color='r',label='Rf = 1.0*S')
    RtoS = 1.3
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj*6e7,color='k',label='Rf = 1.3*S')
    RtoS = 2.0
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj*6e7,'--',color='b',label='Rf = 2.0*S')
    ax1.set_xlabel('Well Spacing (m)',fontsize=10)
    ax1.set_ylabel('Injection Rate (mL/min)',fontsize=10)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    plt.title('Granite with water and varied fracture radius')
    ax1.legend(loc='lower right', prop={'size':10}, ncol=2, numpoints=1)   
    plt.tight_layout()
    
    #**********************************************************************
    #granite with hot water
    #**********************************************************************
    #with variable well diameter
    figsize = (7.5,6)
    f_Pc=1.0e6; p_whp=0.5e6; f_mu=0.1e-3
    RtoS=1.00; LtoS=10.0; f_rho=790.0; w_f=80.0; f_n=3; w_n=3; r_E=50.0e9; r_v=0.30; r_rho=2800.0; r_Kn=0.5
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    w_R = 0.00076
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj,color='m',label='Rw = 0.0008 (1/8" tubing)')
    w_R = 0.0028
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj,color='k',label='Rw = 0.0028 (7/32" lab openhole)')
    w_R = 0.0480
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj,color='b',label='Rw = 0.0480 (HQ EGS Collab)')
    w_R = 0.5*10*0.0254
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.w_S,m.qinj,color='g',label='Rw = 0.1270 (10" commercial EGS)')
    xs = [0.0381,0.0381]
    ys = [8.333e-9,1.333e-7]
    ax1.plot(xs,ys,'-x',color='grey',label='labscale')
    xs = np.asarray([200.0,800.0,800.0,200.0,200.0])
    ys = np.asarray([0.005,0.005,0.08,0.08,0.005])
    ax1.fill_between(xs[0:2],ys[0:2],ys[2:4],color='lightseagreen',alpha=0.5, label='commercial scale EGS')
    # ax1.plot(xs,ys,':',color='grey',label='field scale EGS')
    ax1.set_xlabel('Well Spacing (m)',fontsize=10)
    ax1.set_ylabel('Injection Rate (m3/s)',fontsize=10)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    plt.title('Granite and 260 C water with varied well diameter')
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()
    
    #example single profile (fieldscale)
    figsize = (2.5,2.5)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    i = np.where(m.w_S>500.0)[0][0] #int(0.30*nel)
    xs = m.Ps[0,:,i]
    ys = m.Ps[1,:,i]
    ax1.plot(xs,ys*10**-6,color='grey',label='pressure head')
    ax1.plot(xs,ys*10**-6,'.',color='grey')
    xs = m.Ps[0,:,i]
    ys = m.Ps[2,:,i]
    ax1.plot(xs,ys*10**-6,color='black',label='absolute pressure')
    ax1.plot(xs,ys*10**-6,'.',color='k')
    ax1.set_xlabel('Flow Length (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.title('Granite with 260 C water %.0f m spacing %.1e m3/s' %(m.w_S[i],m.qinj[i]*6e7))
    # ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()
    
    #example single profile (labscale)
    figsize = (2.5,2.5)
    w_R = 0.00076
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    i = np.where(m.w_S>0.0381)[0][0] #int(0.30*nel)
    xs = m.Ps[0,:,i]
    ys = m.Ps[1,:,i]
    ax1.plot(xs,ys*10**-6,color='grey',label='pressure head')
    ax1.plot(xs,ys*10**-6,'.',color='grey')
    xs = m.Ps[0,:,i]
    ys = m.Ps[2,:,i]
    ax1.plot(xs,ys*10**-6,color='black',label='absolute pressure')
    ax1.plot(xs,ys*10**-6,'.',color='k')
    ax1.set_xlabel('Flow Length (m)',fontsize=10)
    ax1.set_ylabel('Pressure (MPa)',fontsize=10)
    plt.title('Granite with 260 C water %.0f m spacing %.1e m3/s' %(m.w_S[i],m.qinj[i]*6e7))
    # ax1.legend(loc='upper right', prop={'size':10}, ncol=1, numpoints=1)
    plt.tight_layout()
    
    #with variable backpressure
    figsize = (7.5,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    p_whp = 0.1e6; w_R = 0.5*10*0.0254
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.r_d,m.qinj,color='m',label='0.1 MPa (atmospheric vent)')
    p_whp = 0.5e6
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.r_d,m.qinj,color='k',label='0.5 MPa (laboratory control)')
    p_whp = 5.0e6
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.r_d,m.qinj,color='b',label='5.0 MPa (compressed water)')
    # xs = [0.0381,0.0381]
    # ys = [8.333e-9,1.333e-7]
    # ax1.plot(xs,ys,'-x',color='grey',label='labscale')
    xs = np.asarray([2000.0,8000.0,8000.0,2000.0,2000.0])
    ys = np.asarray([0.005,0.005,0.08,0.08,0.005])
    ax1.fill_between(xs[0:2],ys[0:2],ys[2:4],color='lightseagreen',alpha=0.5,label='commercial scale EGS')
    # ax1.plot(xs,ys,':',color='grey',label='field scale EGS')
    ax1.set_xlabel('Reservoir Depth (m)',fontsize=10)
    ax1.set_ylabel('Injection Rate (m3/s)',fontsize=10)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    plt.title('Granite with 260 C water and varied backpressure')
    ax1.legend(loc='lower right', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()
    figsize = (5,5)
    
    #with variable critical pressure
    figsize = (7.5,6)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    p_whp = 5.0e6; f_Pc = 0.1e6
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.r_d,m.qinj,color='m',label='0.1 MPa (low)')
    f_Pc = 1.0e6
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.r_d,m.qinj,color='k',label='1.0 MPa (medium)')
    f_Pc = 2.3e6
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.r_d,m.qinj,color='g',label='2.3 MPa (lab)')
    f_Pc = 10.0e6
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)
    ax1.plot(m.r_d,m.qinj,color='b',label='10.0 MPa (high)')
    # xs = [0.0381,0.0381]
    # ys = [8.333e-9,1.333e-7]
    # ax1.plot(xs,ys,'-x',color='grey',label='labscale')
    xs = np.asarray([2000.0,8000.0,8000.0,2000.0,2000.0])
    ys = np.asarray([0.005,0.005,0.08,0.08,0.005])
    ax1.fill_between(xs[0:2],ys[0:2],ys[2:4],color='lightseagreen',alpha=0.5,label='commercial scale EGS')
    # ax1.plot(xs,ys,':',color='grey',label='field scale EGS')
    ax1.set_xlabel('Reservoir Depth (m)',fontsize=10)
    ax1.set_ylabel('Injection Rate (m3/s)',fontsize=10)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    plt.title('Granite with 260 C water and varied critical pressure')
    ax1.legend(loc='lower right', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()
    figsize = (5,5)
    
    #**********************************************************************
    #fieldscale subcritical flow rates
    #**********************************************************************
    p_whp = np.asarray([0.1e6, 0.5e6, 5.0e6])
    w_S = 500.0; w_R = 0.5*10*0.0254; f_Pc=1.0e6; f_mu=0.1e-3
    RtoS=1.05; LtoS=10.0; f_rho=790.0; w_f=80.0; f_n=3; w_n=3; r_E=50.0e9; r_v=0.30; r_rho=2800.0; r_Kn=0.5
    m = cagelimit(w_S,w_R,RtoS,LtoS,f_rho,f_mu,w_f,p_whp,f_Pc,f_n,w_n,r_E,r_v,r_rho,r_Kn)

    figsize = (4,4)
    fig = pylab.figure(figsize=figsize,dpi=dpi)
    ax1 = fig.add_subplot(111)
    ax1.plot(m.Qq[:,0,2],m.Pq[:,0,2]*1e-6,color='k',label='injection surface 10" well')
    ax1.plot([0.001,m.Qq[0,0,2]],[5,m.Pq[0,0,2]*1e-6],color='k')
    ax1.plot(m.Qq[:,0,2],m.Pq[:,1,2]*1e-6,color='r',label='fracture pressure')
    ax1.plot([0.001]+list(m.Qq[:,0,2]),[5]+list(m.Pq[:,1,2]*1e-6),color='r')
    ax1.plot([0.001]+list(m.Qq[:,0,2]),[5]+list(m.Pq[:,5,2]*1e-6),color='b',label='producer surface')
    xs = [m.s3*1e-6,m.s3*1e-6]
    ys = [0.001,np.nanmax(m.Qq[:,0,2])]
    ax1.plot(ys,xs,color='pink',label='fracture opening')
    xs = [0.005,1.0]
    ys1 = [5.0,5.0]
    ys2 = [105.0,105.0]
    ax1.fill_between(xs,ys1,ys2,color='black',alpha=0.3,label='feasibility range')
    ax1.plot(m.qinj[2],m.Ps[2,0,2]*1e-6,marker='*',linestyle='',color='k',label='caging limit')

    ax1.set_ylabel('Absolute Pressure (MPa)',fontsize=10)
    ax1.set_xlabel('Flow Rate (m3/s)',fontsize=10)
    # ax1.set_yscale('log')
    ax1.set_xscale('log')
    plt.title('Caged sub-limit-flow at %.0f m depth' %(m.r_d))
    ax1.legend(loc='upper left', prop={'size':10}, ncol=1, numpoints=1)   
    plt.tight_layout()    
    figsize = (5,5)

pylab.show()
