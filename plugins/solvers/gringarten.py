# -*- coding: utf-8 -*-
"""
A lookup table based solver for (Gringarten et al, 1975: Journal of Geophyscial Research)
"""
print('gringarten_1975')

import numpy as np
import pylab
import matplotlib.pyplot as plt

yr=365.2425*24.0*60.0*60.0 #s
Darcy = 1.0e-12 #m2
MPa=10.0**6.0#Pa
cP=10.0**-3#Pa-s

class gringarten:
    #solution table
    def __init__(self):
        self.Twd = np.asarray([0.001]+list(np.linspace(0.01,0.99,41-2))+[0.999])
        self.Xed = np.asarray([0.5,1,2,4,8,16])
        self.tD = np.asarray([
               [2.51720714e-01, 3.21921259e-01, 4.08331559e-01, 4.65821224e-01,
                5.08740971e-01, 5.51660719e-01, 5.85676667e-01, 6.18587289e-01,
                6.51497911e-01, 6.78794807e-01, 7.02143669e-01, 7.25492532e-01,
                7.48841395e-01, 7.73145231e-01, 7.97829264e-01, 8.22513296e-01,
                8.47197329e-01, 8.71881361e-01, 8.96565394e-01, 9.21178786e-01,
                9.45766924e-01, 9.70355062e-01, 9.94943200e-01, 1.01953134e+00,
                1.04411948e+00, 1.06870761e+00, 1.09882831e+00, 1.13325529e+00,
                1.16768227e+00, 1.20210925e+00, 1.23653624e+00, 1.27096322e+00,
                1.31418686e+00, 1.36027226e+00, 1.40635766e+00, 1.45244306e+00,
                1.50315367e+00, 1.55610266e+00, 1.70803421e+00, 2.10915991e+00,
                2.54626769e+00],
               [2.51720714e-01, 3.21921259e-01, 4.56261377e-01, 5.71801668e-01,
                6.68483509e-01, 7.56497677e-01, 8.42183343e-01, 9.03757431e-01,
                9.69105215e-01, 1.03445300e+00, 1.09980078e+00, 1.15987242e+00,
                1.21634130e+00, 1.27281019e+00, 1.32927907e+00, 1.39139075e+00,
                1.45926797e+00, 1.52714519e+00, 1.59502241e+00, 1.66289964e+00,
                1.73077686e+00, 1.80823562e+00, 1.89006313e+00, 1.97189063e+00,
                2.05371814e+00, 2.13554564e+00, 2.21737315e+00, 2.31459748e+00,
                2.41287478e+00, 2.51115209e+00, 2.60942940e+00, 2.70770670e+00,
                2.83522293e+00, 2.98359129e+00, 3.13195965e+00, 3.28032801e+00,
                3.49183451e+00, 3.95683297e+00, 4.42183142e+00, 6.28759212e+00,
                9.27279569e+00],
               [2.51720714e-01, 3.21921259e-01, 4.56261377e-01, 5.71801668e-01,
                6.68483509e-01, 7.82583276e-01, 8.99817649e-01, 1.02359912e+00,
                1.14738059e+00, 1.28070443e+00, 1.42248485e+00, 1.56426527e+00,
                1.70267065e+00, 1.84096412e+00, 1.97925759e+00, 2.11875348e+00,
                2.28779308e+00, 2.45683268e+00, 2.62587228e+00, 2.79491189e+00,
                2.99584138e+00, 3.20968661e+00, 3.42353185e+00, 3.63737708e+00,
                3.85122232e+00, 4.11229951e+00, 4.39946799e+00, 4.68663646e+00,
                4.97380494e+00, 5.26097342e+00, 5.54814190e+00, 5.96652357e+00,
                6.42203581e+00, 6.87754804e+00, 7.33306028e+00, 7.78857252e+00,
                8.78232394e+00, 9.88785663e+00, 1.21894446e+01, 1.82024539e+01,
                2.33885711e+01],
               [ 0.25172071,  0.32192126,  0.45626138,  0.57180167,  0.66848351,
                0.78258328,  0.89981765,  1.02359912,  1.14738059,  1.28070443,
                1.42248485,  1.56426527,  1.76808917,  1.97397044,  2.18366964,
                2.44919547,  2.71472129,  3.02972873,  3.44477476,  3.8598208 ,
                4.27943537,  4.74797176,  5.21650815,  5.68504454,  6.15358093,
                6.82894963,  7.61857405,  8.40819848,  9.1978229 , 10.2470986 ,
                11.39968659, 12.55227458, 13.70486258, 15.26561583, 17.47019012,
                19.67476442, 22.09323585, 25.4888283 , 31.29907978, 49.82152011,
                72.40667386],
               [2.51720714e-01, 3.21921259e-01, 4.56261377e-01, 5.71801668e-01,
                6.68483509e-01, 7.82583276e-01, 8.99817649e-01, 1.02359912e+00,
                1.14738059e+00, 1.28070443e+00, 1.42248485e+00, 1.56426527e+00,
                1.76808917e+00, 1.97397044e+00, 2.18366964e+00, 2.44919547e+00,
                2.71472129e+00, 3.02972873e+00, 3.44477476e+00, 3.85982080e+00,
                4.29702902e+00, 4.97155933e+00, 5.64608963e+00, 6.32061994e+00,
                7.35253593e+00, 8.45258550e+00, 9.55263507e+00, 1.13725872e+01,
                1.33483132e+01, 1.54757882e+01, 1.84910775e+01, 2.15063668e+01,
                2.45216561e+01, 2.98863640e+01, 3.55108240e+01, 4.11352839e+01,
                5.17659489e+01, 6.44229953e+01, 8.56800985e+01, 1.76698341e+02,
                3.32669834e+02],
               [2.51720714e-01, 3.21921259e-01, 4.56261377e-01, 5.71801668e-01,
                6.68483509e-01, 7.82583276e-01, 8.99817649e-01, 1.02359912e+00,
                1.14738059e+00, 1.28070443e+00, 1.42248485e+00, 1.56426527e+00,
                1.76808917e+00, 1.97397044e+00, 2.18366964e+00, 2.44919547e+00,
                2.71472129e+00, 3.02972873e+00, 3.44477476e+00, 3.85982080e+00,
                4.29702902e+00, 4.97155933e+00, 5.64608963e+00, 6.32061994e+00,
                7.35253593e+00, 8.45258550e+00, 9.55263507e+00, 1.13725872e+01,
                1.33483132e+01, 1.56764281e+01, 2.00662096e+01, 2.44559911e+01,
                3.34785899e+01, 4.32031938e+01, 6.37846668e+01, 1.00802797e+02,
                1.78380331e+02, 3.67364987e+02, 2.49606527e+03, 2.80000000e+04,
                6.58000000e+04]])
    #Temperature vs time for plotting
    # def Tvt(self,fracspacing=3.519,Qinj=7e-5,Kt_r=3.2,Sv_r=2357.6,
    #         Sv_f=4033.7,L=12.5,W=12.5,Tinj=10.0,Trock=30.7):
    def Tvt(self,fracspacing=50.0,Qinj=0.1/12,Kt_r=2.5,Sv_r=2063.0,
            Sv_f=4030.0,L=200.0,W=200.0,Tinj=323.1,Trock=523.1):
        #scaling parameters
        self.C1 = 1e3*Sv_f**2.0/(Kt_r*Sv_r)
        self.C2 = 1e3*Sv_f/Kt_r
        self.Xe = 0.5*fracspacing*self.C2*Qinj/(L*W)
        #interpolated scaled-time vs spacing
        # for t in range(0,len(self.Twd)):
        if self.Xe > self.Xed[-1]:
            self.td = self.tD[-1,:]
        elif self.Xe < self.Xed[0]:
            self.td = self.tD[0,:]
        else:
            j = np.where(self.Xe < self.Xed)[0][0]
            i = j-1
            self.td = ((self.Xe-self.Xed[j])*(self.tD[i,:]-self.tD[j,:])/(self.Xed[i]-self.Xed[j])+self.tD[j,:])
        #real units
        self.ts = self.td/(self.C1*(Qinj/(L*W))**2)
        self.Ts = Trock-self.Twd*(Trock-Tinj)
    def plotTd(self):
        fig = plt.figure(figsize=(10.0,5.0),dpi=100)
        ax1 = fig.add_subplot(111)
        for i in range(0,len(self.Xed)):
            ax1.plot(self.tD[i],self.Twd,label='Xe %.1f' %(self.Xed[i]))
        ax1.plot(self.td,self.Twd,':')
        ax1.set_xlabel('Dimensionless Time',fontsize=10)
        ax1.set_xscale('log')
        ax1.set_ylabel('Dimensionless Temp',fontsize=10)
        ax1.set_ylim([1,0])
        ax1.legend(loc='lower left',ncol=2,fontsize=9)
        plt.tight_layout()
    def plotTs(self):
        fig = plt.figure(figsize=(10.0,5.0),dpi=100)
        ax1 = fig.add_subplot(111)
        ax1.plot(self.ts/yr,self.Ts,label='pre-interpolation')
        ax1.plot(self.time/yr,self.temp,label='post-interpolation')
        ax1.set_xlabel('Time (yr)',fontsize=10)
        ax1.set_xscale('log')
        ax1.set_ylabel('Temp (K)',fontsize=10)
        ax1.legend(loc='lower left',ncol=2,fontsize=9)
        plt.tight_layout()
    def Tvts(self,time):
        #interpolated temperature vs passed in time
        self.time = time
        self.temp = np.zeros(len(time))
        for t in range(0,len(time)):
            if time[t] > self.ts[-1]:
                self.temp[t] = self.Ts[-1]
            elif time[t] < self.ts[0]:
                self.temp[t] = self.Ts[0]
            else:
                j = np.where(time[t] < self.ts)[0][0]
                i = j-1
                self.temp[t] = ((time[t]-self.ts[j])*(self.Ts[i]-self.Ts[j])/
                           (self.ts[i]-self.ts[j])+self.Ts[j])        
    def plottime(self):
        fig = plt.figure(figsize=(10.0,5.0),dpi=100)
        ax1 = fig.add_subplot(111)
        ax1.plot(self.time/yr,self.temp)
        ax1.set_xlabel('Time (yr)',fontsize=10)
        ax1.set_xscale('log')
        ax1.set_ylabel('Temp (K)',fontsize=10)
        plt.tight_layout()
    def hvT(self,T):
        h = 0.0025*T**2 + 2.27*T - 770
        return h
    def Pvh(self,h,q,P1):
        effic = 7.8795*np.log(h) - 45.651 #(Zarrouk and Moon, 2014: Geothermics)
        power = (effic/100)*(h-self.hvT(50+273.15))*q*965 - P1*q*1e-3
        power[power < 0] = 0.0
        return power
    # def simple_EGS_demo(self,numfracs=7,fracspacing=50,fracaperture=0.001,
    #                     wellspacing=250,welldiameter=6*0.0254,roughness=80.0,
    #                     depth=3500,gradient=80,homogeneity=0.25,numwells=3,
    #                     proppantconductivity=100*Darcy,backpressure=5*MPa):
    # def simple_EGS_demo(self,numfracs=4,fracspacing=200,fracaperture=0.001,
    #                     wellspacing=400,welldiameter=10*0.0254,roughness=80.0,
    #                     depth=4000,gradient=62.5,homogeneity=1.0,numwells=3,
    #                     proppantconductivity=100*Darcy,backpressure=5*MPa):
    def EGS_demo(self,numfracs=50,fracspacing=50,fracaperture=0.001,
            wellspacing=110,welldiameter=8*0.0254,roughness=80.0,
            depth=2300,gradient=80.0,homogeneity=0.1,numwells=2,
            proppantconductivity=50*Darcy,backpressure=1*MPa):
        L = wellspacing
        W = wellspacing
        mu = 0.1*cP
        g = 9.81
        rho = 965.0
        #subsurface conditions
        openingpressure = 0.5*2700*g*depth
        staticpressure = rho*g*depth
        rocktemperature = gradient*depth/1e3 + 273.15
        #pressure based on parallel circuit model with no chokes
        def PvQ(Qinj):
            heterogeneity = np.linspace(homogeneity,1,numfracs)
            frac_K = (mu*L)/(rho*g*W*proppantconductivity*fracaperture*heterogeneity)
            well_K = 10.7*(depth*mu/(0.9*cP))/(roughness**1.852*welldiameter**4.87)
            fracs_K = 1.0/np.sum(1.0/frac_K)
            dPinj = rho*g*well_K*Qinj**1.852
            dPpro = rho*g*well_K*(Qinj/(numwells-1))**1.852
            dPfra = rho*g*fracs_K*Qinj
            P4 = backpressure
            P3 = P4+dPpro
            P2 = P3+dPfra
            P1 = P2+dPinj
            q_fracs = dPfra/(frac_K*rho*g)
            return P1, P2, P3, P4, q_fracs
        #solve heat and pressure for range of flow rates
        fig = plt.figure(figsize=(10.0,9.0),dpi=100)
        ax0 = fig.add_subplot(221)
        ax = fig.add_subplot(222)
        ax2 = fig.add_subplot(223)
        ax3 = fig.add_subplot(224)
        Qs = 10**(np.linspace(np.log10(0.01),np.log10(1.0),14))
        #Qs = np.asarray([0.01,0.04,0.08])
        timepassed = np.linspace(0,1.5,101)*yr
        for q in range(0,len(Qs)):
            #solve flow
            P1, P2, P3, P4, q_fracs = PvQ(Qs[q])
            ax2.plot([1,2,3,4],[P1/MPa,(P2+staticpressure)/MPa,(P3+staticpressure)/MPa,P4/MPa],label='test %i: %.3f m3/s' %(q, Qs[q]))
            text = 'Case %i (%.3f m3/s):' %(q,Qs[q])
            text += ' P1=%.1f' %(P1/MPa)
            text += ' P2=%.1f' %(P2/MPa)
            text += ' P3=%.1f' %(P3/MPa)
            text += ' P4=%.1f' %(P4/MPa)
            #text += ', sub-flows'
            T_fracs = []
            T_mixed = np.zeros(len(timepassed))
            for f in range(0,len(q_fracs)):
                #text += ' %.1e' %(q_fracs[f])
                #solve heat
                self.Tvt(fracspacing=fracspacing,Qinj=q_fracs[f]/(numwells-1),L=L,
                                    W=W,Tinj=50.0+273.15,Trock=rocktemperature)
                self.Tvts(timepassed)
                T_fracs += [self.temp]
                T_mixed += self.temp*(q_fracs[f]/Qs[q])
                if q == 11:
                    print(self.temp)
                    ax0.plot(self.time/yr,self.temp,label='test %i: frac %i' %(q, f))
                    
            ax.plot(self.time/yr,T_mixed,label='test %i: mixed' %(q))
            if q == 11:
                ax0.plot(self.time/yr,T_mixed,'--',label='test %i: mixed' %(q))
                td = self.time[-1]*self.C1*(q_fracs[-1]/(L*W))**2
                print('max dimensionless time [%i] %.2e and dimensionless spacing of %.2e' %(q,td,2*self.Xe))
            #solve power
            h = self.hvT(T_mixed)
            P = self.Pvh(h,Qs[q],P1)
            ax3.plot(self.time,P,label='test %i: %.3f m3/s' %(q, Qs[q]))
            text += ', average net-power %.3f kWe' %(np.sum(P)/len(P))
            print(text)
        ax.plot(self.time/yr,np.ones(len(self.time))*(273.15+80),'--',label='minimum for electricity')
        ax.set_xlabel('Time (yr)',fontsize=10)
        ax.set_ylabel('Temp (K)',fontsize=10)
        ax.legend(loc='upper right',ncol=2,fontsize=9)
        ax2.plot([1,4],[openingpressure/MPa,openingpressure/MPa],'--',label='max allowable at P2')
        ax2.set_xlabel('Position',fontsize=10)
        ax2.set_ylabel('Abs Pressure (MPa)',fontsize=10)
        ax2.legend(loc='upper right',ncol=2,fontsize=9)
        ax0.set_xlabel('Time (yr)',fontsize=10)
        ax0.set_ylabel('Temperature (K)',fontsize=10)
        ax0.legend(loc='upper right',ncol=2,fontsize=9)
        ax3.plot(self.time/yr,np.ones(len(self.time))*(1e3*numwells),'--',label='minimum for viability')
        ax3.set_xlabel('Time (yr)',fontsize=10)
        ax3.set_ylabel('Net Power (kWe)',fontsize=10)
        ax3.set_ylim([0,10e3])
        ax3.legend(loc='upper right',ncol=2,fontsize=9)
        plt.tight_layout()
        plt.show()
    
    #standard Gringarten solution without pressure or heterogeneity modeling
    def basic(self,rocktemp=473.15,rockKt=2.306,rockSv=2133,
              numwells=2,wellspacing=110,numfracs=50,fracspacing=50,
              tinj=274.0,bulkflow=0.04,tortuosity=1.0,poreSv=3320.0,
              timepassed=[0.0e0,1.6e7,3.2e7,4.7e7,6.3e7]):
        #setup
        L = wellspacing
        W = wellspacing
        #flow distribution from cubic law
        heterogeneity = np.linspace(0,1,int(numfracs))**3*(1-tortuosity)+tortuosity
        q_fracs = (bulkflow*heterogeneity**3.0)/np.sum(heterogeneity**3.0)
        #solve heat transfer for each fracture
        T_mixed = np.zeros(len(timepassed))
        for f in range(0,len(q_fracs)):
            self.Tvt(fracspacing=fracspacing,Qinj=q_fracs[f]/(numwells-1),
                     Kt_r=rockKt,Sv_r=rockSv,Sv_f=poreSv,
                     L=L,W=W,Tinj=tinj,Trock=rocktemp)
            self.Tvts(timepassed)
            T_mixed += self.temp*(q_fracs[f]/bulkflow)
        self.Tout = T_mixed
        # self.hout = self.hvT(T_mixed)
        # self.Pout = self.Pvh(self.hout,bulkflow,P1)
    
    #gringarten solver validation test using Gringarten, 1975, Fig. 5's inputs
    # - assumes uniform flow
    # - includes effects of fracture spacing
    # - no natural fractures
    # - plot results
    def verify(self,tortuosity=1.0,visuals=True):
        #setup
        fig = plt.figure(figsize=(10.0,5.0),dpi=100)
        ax = fig.add_subplot(111)
        ty = np.linspace(0,100,101)
        tp = np.asarray(ty*365.25*24*60*60)
        #10 fractures, 160 m spacing
        fn = 10; fs = 160.0
        self.basic(rocktemp=300.0,
                   rockKt=2.594,
                   rockSv=1.046*2650,
                   numwells=2,
                   wellspacing=1000.0,
                   numfracs=fn,
                   fracspacing=fs,
                   tinj=65.0,
                   bulkflow=0.145,
                   tortuosity=tortuosity,
                   poreSv=4.184*1000,
                   timepassed=tp)
        ax.plot(ty,self.Tout,label='%i fractures, %.0f m spacing, %0.2f Xe' %(fn,fs,self.Xe)) 
        #10 fractures, 80 m spacing
        fn = 10; fs = 80.0
        self.basic(rocktemp=300.0,
                    rockKt=2.594,
                    rockSv=1.046*2650,
                    numwells=2,
                    wellspacing=1000.0,
                    numfracs=fn,
                    fracspacing=fs,
                    tinj=65.0,
                    bulkflow=0.145,
                    tortuosity=tortuosity,
                    poreSv=4.184*1000,
                    timepassed=tp)
        ax.plot(ty,self.Tout,label='%i fractures, %.0f m spacing, %0.2f Xe' %(fn,fs,self.Xe)) 
        #10 fractures, 40 m spacing
        fn = 10; fs = 40.0
        self.basic(rocktemp=300.0,
                   rockKt=2.594,
                   rockSv=1.046*2650,
                   numwells=2,
                   wellspacing=1000.0,
                   numfracs=fn,
                   fracspacing=fs,
                   tinj=65.0,
                   bulkflow=0.145,
                   tortuosity=tortuosity,
                   poreSv=4.184*1000,
                   timepassed=tp)
        ax.plot(ty,self.Tout,label='%i fractures, %.0f m spacing, %0.2f Xe' %(fn,fs,self.Xe)) 
        #1 fractures, 500 m spacing
        fn = 1; fs = 500.0
        self.basic(rocktemp=300.0,
                   rockKt=2.594,
                   rockSv=1.046*2650,
                   numwells=2,
                   wellspacing=1000.0,
                   numfracs=fn,
                   fracspacing=fs,
                   tinj=65.0,
                   bulkflow=0.145,
                   tortuosity=tortuosity,
                   poreSv=4.184*1000,
                   timepassed=tp)
        ax.plot(ty,self.Tout,label='%i fractures, %.0f m spacing, %0.2f Xe' %(fn,fs,self.Xe))
        #labels
        ax.set_xlabel('Time (yr)',fontsize=10)
        ax.set_ylabel('Temp (K)',fontsize=10)
        ax.set_title('%0.2e C1, %.02e C2' %(self.C1,self.C2))
        ax.legend(loc='lower center',ncol=2,fontsize=10)
        plt.tight_layout()
        if visuals:
            pylab.show()
        
#validation plots
if False:
    x = gringarten()
    x.verify(tortuosity=1.0,visuals=True)

#EGS demo
if False:
    x = gringarten()
    x.Tvt()
    x.plotTs()
    times = 10**np.linspace(np.log10(0.1),np.log10(15),41)*yr
    x.Tvts(times)
    x.plottime()
    x.simple_EGS_demo()
    plt.show()
        