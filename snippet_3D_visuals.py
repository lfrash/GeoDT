# -*- coding: utf-8 -*-
"""
Code snippets to create 3D visuals from setup_payoff.csv
"""

#***************************************************************************************************************
### load file
#***************************************************************************************************************
from plugins.solvers.libs import csv
from plugins.solvers.libs import stress
from plugins.solvers.libs import stats
names, data = csv.load_csv('setup_payoff.csv')

#***************************************************************************************************************
#from: wells (i.e., line segments that represent wells)
#***************************************************************************************************************
#units
import numpy as np
deg=1.0*np.pi/180.0

#line object
class line:
    def __init__(self,x0=0.0,y0=0.0,z0=0.0,length=1.0,azn=0.0*deg,dip=0.0*deg,w_type='pipe',
                 ra=0.0254*3.0,rb=0.0254*3.5,rc=0.0254*3.5,rough=80.0,pID=0):
        #position geometry
        self.c0 = np.asarray([x0, y0, z0]) #origin
        self.leg = length #length
        self.azn = azn #axis azimuth north
        self.dip = dip #axis dip from horizontal
        self.type = w_type #type of well
        #flow geometry
        self.ra = ra #radius, m
        self.rb = rb #radius, m
        self.rc = rc #radius, m
        self.rgh = rough
        #stimulation traits
        self.hydrofrac = False #was well already hydrofraced?
        self.completed = False #was well stimulation process completed?
        self.stabilize = False #was well flow stabilied?
        #parent well index
        self.pID = pID
        #dependent geometry
        self.get_c1()
    def get_c1(self):
        self.vAxi = np.asarray([np.sin(self.azn)*np.cos(-self.dip), np.cos(self.azn)*np.cos(-self.dip), np.sin(-self.dip)])
        self.vRht = np.asarray([np.sin(self.azn+np.pi/2),np.cos(self.azn+np.pi/2),0.0])
        self.vNor = np.cross(self.vRht,self.vAxi)
        self.vNor = self.vNor/np.linalg.norm(self.vNor)
        self.c1 = self.c0 + self.vAxi*self.leg

#rotation matrix
def rotationMatrix(axis, theta):
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

#rotate points about axis
def rotatePoints(points, axis, theta):
    rot=rotationMatrix(axis,theta)
    #return rot*points
    return np.transpose( np.dot(rot, np.transpose( points ) ) )

#azn and dip from endpoints
def azn_dip(x0,x1):
    dx = x1[0] - x0[0]
    dy = x1[1] - x0[1]
    dz = x1[2] - x0[2]
    dr = (dx**2.0+dy**2.0)**0.5
    dis = (dx**2.0+dy**2.0+dz**2.0)**0.5
    azn = []
    dip = []
    if dx == 0 and dy >= 0:
        azn = 0.0
    elif dx == 0:
        azn = np.pi
    else:
        azn = np.sign(dx)*np.arccos(dy/dr) + (1 - np.sign(dx))*np.pi
    if dr == 0.0:
        dip = -np.sign(dz)*np.pi/2.0
    else:
        dip = -np.arctan(dz/dr)
    return azn, dip, dis

#well geometry (this is essentially a list of well segments (i.e., list of line objects))
def gen_wells(data,row):
    #initialization
    wells = []
    n = len(wells)
    #centerpoint, target depth (z0) = 0.0
    i0 = np.asarray([0.0, 0.0, 0.0])
    #ref axes (parallel injector and 90 horizontal to the right)
    azn = data['Well_Azimuth_rad'][row]
    dip = data['Well_Dip_rad'][row]
    vInj = np.asarray([np.sin(azn)*np.cos(-dip),
             np.cos(azn)*np.cos(-dip),np.sin(-dip)])
    vRht = np.asarray([np.sin(azn+np.pi/2),np.cos(azn+np.pi/2),0.0])
    # vNor = np.asarray([np.sin(azn)*np.sin(-dip),
    #          np.cos(azn)*np.sin(-dip),np.cos(-dip)])
    vNor = np.cross(vRht,vInj)
    vNor = vNor/np.linalg.norm(vNor)
    #offset center
    p0 = i0 + vRht*data['Well_Spacing_m'][row]
    #lengths of wells
    length = data['Well_DeviatedLength_m'][row]
    proportion = data['Well_ProducerProportion_ratio'][row]
    iLen = length*proportion
    pLen = length
    cLen = pLen - iLen
    #segments
    depth = data['Well_TargetDepth_m'][row]
    segle = depth/data['Well_SegmentCount_units'][row]
    num = data['Well_ProducerCount_wells'][row]
    phase = -1.0*data['Well_RotationPhase_rad'][row]
    seg = data['Well_Intervals_count'][row]
    leg = iLen/(seg+seg-1)
    rgh = data['Well_Roughness_metric'][row]
    #conductor segment
    h = np.asarray([0.0,0.0,-segle])
    #radii
    prb = data['Well_ProductionRadius_m'][row]
    pra = prb - data['Well_WallThickness_m'][row]
    prc = prb + data['Well_AnnularThickness_m'][row]
    ira = prc
    irb = ira + data['Well_WallThickness_m'][row]
    irc = irb + data['Well_AnnularThickness_m'][row]
    sra = irc
    srb = sra + data['Well_WallThickness_m'][row]
    src = srb + data['Well_AnnularThickness_m'][row]
    
    # ************************************************************************
    # O&G - wine-rack style wells with casing to bottom in all wells
    # ************************************************************************
    if data['Strategy_Design_type'][row].upper() in ['O&G','DCM']:
        #make wells parallel
        vPro = vInj
        #populate rack of production and injection wells
        i0s = []
        p0s = []
        for i in range(0,num):
            if i == 0: 
                #first production well includes a mirrored well
                p0s += [rotatePoints([p0],vInj,(np.pi - phase))[0]]
                i0s += [i0]
            else: 
                #place additional producers
                p0s += [p0s[-1] + 2.0*vRht*data['Well_Spacing_m'][row]*np.cos(phase)]
            #don't place injector in second step
            if i > 1:
                i0s += [i0s[-1] + 2.0*vRht*data['Well_Spacing_m'][row]*np.cos(phase)]
        #recenter
        c0 = np.zeros(3)
        if num < 2:
            c0 = c0 - 0.5*vRht*data['Well_Spacing_m'][row]*np.cos(phase)
        else:
            c0 = c0 + (num-2)*vRht*data['Well_Spacing_m'][row]*np.cos(phase)
        for i in range(0,len(i0s)):
            i0s[i] += -1.0*c0
        for i in range(0,len(p0s)):
            p0s[i] += -1.0*c0  
        #for each injection well
        i1ss = [] #list of startpoints
        i2ss = [] #list of endpoints
        ibss = [] #list of bottomhole vertical segments
        isss = [] #list of surface vertical segments
        for i in range(0,len(i0s)):
            #injection well reservoir segments
            i1s = [] #startpoints
            i2s = [] #endpoints
            i1s += [i0s[i] - 0.5*vInj*iLen]
            i2s += [i1s[0] + leg*vInj]
            for i in range(1,seg):
                i1s += [i2s[i-1] + leg*vInj]
                i2s += [i1s[i] + leg*vInj]
            #injection well segments to surface
            ibs = i1s[0] - 0.5*cLen*vInj
            iss = np.asarray([ibs[0],ibs[1],depth])
            #store in lists
            i1ss += [i1s]
            i2ss += [i2s]
            ibss += [ibs]
            isss += [iss]
        #production well segments
        p1s = [] #startpoints
        p2s = [] #endpoints
        pss = [] #surface segments
        for i in range(0,num):
            p1s += [p0s[i] - 0.5*vPro*pLen]
            p2s += [p0s[i] + 0.5*vPro*pLen]    
            pss += [np.asarray([p1s[i][0],p1s[i][1],depth+0.1])]
        #place wellheads
        # h = np.asarray([0.0,0.0,s.Well_Spacing_m])
        iid = np.asarray(range(0,len(i0s)))+n
        pid = np.asarray(range(0,len(p0s)))+n+len(i0s)
        for i in range(0,len(i0s)):
            iss = isss[i]
            azn, dip, dis = azn_dip(iss,iss+h)
            wells += [line(iss[0],iss[1],iss[2],dis,azn,dip,'injector',sra,srb,src,rgh,iid[i])]
        for i in range(0,len(p0s)):
            azn, dip, dis = azn_dip(pss[i],pss[i]+h)
            wells += [line(pss[i][0],pss[i][1],pss[i][2],dis,azn,dip,'producer',sra,srb,src,rgh,pid[i])]
        #place surface segments
        for i in range(0,len(i0s)):
            iss = isss[i]
            ibs = ibss[i]
            azn, dip, dis = azn_dip(iss+h,ibs)
            wells += [line(iss[0]+h[0],iss[1]+h[1],iss[2]+h[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,iid[i])]
        for i in range(0,len(p0s)):
            azn, dip, dis = azn_dip(pss[i]+h,p1s[i])
            wells += [line(pss[i][0]+h[0],pss[i][1]+h[1],pss[i][2]+h[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,pid[i])]
        #place injector base segment
        for i in range(0,len(i0s)):
            ibs = ibss[i]
            i1s = i1ss[i]
            i2s = i2ss[i]
            azn, dip, dis = azn_dip(ibs,i1s[0])
            wells += [line(ibs[0],ibs[1],ibs[2],dis,azn,dip,'pipe',pra,prb,prc,rgh,iid[i])]
            #place injection internal segments
            for i in range(0,seg-1):
                azn, dip, dis = azn_dip(i1s[0],i2s[0])
                wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],dis,azn,dip,'perfcluster',pra,prb,prc,rgh,iid[i])]
                azn, dip, dis2 = azn_dip(i2s[0],i1s[1])
                wells += [line(i2s[i][0],i2s[i][1],i2s[i][2],dis2,azn,dip,'pipe',pra,prb,prc,rgh,iid[i])]
            azn, dip, dis = azn_dip(i1s[-1],i2s[-1])
            wells += [line(i1s[-1][0],i1s[-1][1],i1s[-1][2],dis,azn,dip,'perfcluster',pra,prb,prc,rgh,iid[i])]
        #place production wells
        for i in range(0,len(p0s)):
            azn, dip, dis = azn_dip(p1s[i],p2s[i])
            wells += [line(p1s[i][0],p1s[i][1],p1s[i][2],dis,azn,dip,'screen',pra,prb,prc,rgh,pid[i])]
        #stimulate all wells for O&G
        if data['Strategy_Design_type'][row].upper() == 'O&G':
            for w in wells:
                if w.type == 'screen':
                    w.type = 'procluster'
    #return result
    return wells

#***************************************************************************************************************
###from: fractures 
#***************************************************************************************************************
#fracture object
class surf:
    def __init__(self,c0=[0,0,0],dia=1.0,stk=0.0*deg,dip=90.0*deg,
                 ty='fracture', data=[], row=0):
        #node number of center point
        self.ci = -1
        self.wi = -1
        
        #geometry
        self.c0 = np.asarray(c0,dtype=float)
        self.dia = dia
        self.str = stk
        self.dip = dip
        self.type = ty
        
        #properties
        self.En = data['Rock_YoungsModulus_Pa'][row]
        self.vn = data['Rock_PoissonRatio_ratio'][row]
        self.roughness = data['Fracture_Tortuosity_ratio'][row] #stats.loguniform(1,s.Fracture_Tortuosity_ratio,1.0)[0]
        self.alpha = data['Tension_FractureCompressibility_1pPa'][row] #* (1.0 + stats.norm_trunc(1,0.0,0.05,-0.12,0.12)[0]*s.Strategy_SecondOrderNoise_ratio)
        self.prop_alpha = data['Stimulation_ProppantCompressibility_1pPa'][row] #* (1.0 + stats.norm_trunc(1,0.0,0.05,-0.12,0.12)[0]*s.Strategy_SecondOrderNoise_ratio)
        self.kf = 10.0**(np.log10(data['Stimulation_ProppantPermeability_m2'][row]) ) #* (1.0 + stats.norm_trunc(1,0.0,0.05,-0.12,0.12)[0]*s.Strategy_SecondOrderNoise_ratio))
        self.gamma = 10.0**(np.log10(data['Shear_SlipLength_ratio'][row]) ) #* (1.0 + stats.norm_trunc(1,0.0,0.2,-0.4,0.5)[0]*s.Strategy_SecondOrderNoise_ratio))
        self.psi = data['Shear_DilationSlip_ratio'][row] #* (1.0 + stats.norm_trunc(1,0.0,0.2,-0.46,0.35)[0]*s.Strategy_SecondOrderNoise_ratio) #a
        self.eta = data['Shear_HydraulicDilation_ratio'][row] #* (1.0 + ((stats.contact_trunc(1,0.15,1.0,0.5,0.5,0.15,0.0,2.0))[0]-1.0)*s.Strategy_SecondOrderNoise_ratio) #n
        self.bd = data['JointSets_InitialDilation_m'][row] #* (1.0 + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.Strategy_SecondOrderNoise_ratio)
        
        #friction
        if ty in ['propped']:
            self.mu = data['Tension_FrictionCoefficient_ratio'][row] #* (1.0 + stats.norm_trunc(1,0.0,0.04,-0.08,0.10)[0]*s.Strategy_SecondOrderNoise_ratio)
            self.phi = np.arctan(self.mu)
            self.mcc = data['Tension_Cohesion_Pa'][row] #* (1.0 + stats.norm_trunc(1,0.0,0.2,-0.5,1.5)[0]*s.Strategy_SecondOrderNoise_ratio)
            # self.bd = 0.0 #TODO: is this needed?
        else:
            self.mu = np.arctan(data['Shear_FrictionCoefficient_ratio'][row]) #* (1.0 + stats.norm_trunc(1,0.0,0.12,-0.28,0.36)[0]*s.Strategy_SecondOrderNoise_ratio))
            self.phi = np.arctan(self.mu)
            self.mcc = data['Shear_Cohesion_Pa'][row] #* (1.0 + stats.norm_trunc(1,0.0,0.2,-0.5,1.5)[0]*s.Strategy_SecondOrderNoise_ratio)
        
        #stress state #TODO: eliminate subclass?
        self.tensor = stress.cauchy()
        self.tensor.sigG = np.asarray([[s.Tensor_00_Pa,s.Tensor_01_Pa,s.Tensor_02_Pa],
                                       [s.Tensor_10_Pa,s.Tensor_11_Pa,s.Tensor_12_Pa],
                                       [s.Tensor_20_Pa,s.Tensor_21_Pa,s.Tensor_22_Pa]])
        self.Pc, self.sn, self.tau = self.tensor.Pc_frac(self.str, self.dip, self.phi, self.mcc)
        
        #trackers and dependants
        self.stim = 0 #count
        self.Pmax = 0.0 #Pa
        self.Pcen = 0.0 #Pa
        self.Mws = [] #Mw
        self.arup = 1.0 #m2
        self.hydroprop = False
        self.prop_load = 0.0 #m3
        self.bh = self.bd*self.eta
        self.bd0 = self.bd
        self.bd0p = 0.0
        self.vol = (4.0/3.0)*np.pi*0.25*self.dia**2.0*0.5*self.bd
        self.arup = 0.25*np.pi*self.dia**2.0
        self.fc = 3.28084*1.01324e15*(self.bh**3.0)/12.0

    #adjust fracture cohesion to prevent runaway stimulation at specified conditions
    def check_integrity(self,pres=0.0):
        #originally assigned values
        Pc0 = self.Pc
        mcc0 = self.mcc
        phi0 = self.phi
        #shear stability limit
        mcc_c = self.tau - (self.sn - pres)*self.mu
        #tensile stability limit
        mcc_t = pres - self.sn
        #update fracture cohesion to ensure stability at the input conditions
        self.mcc = np.max([self.mcc,mcc_c,mcc_t])
        #recompute critical conditions
        self.Pc, self.sn, self.tau = self.tensor.Pc_frac(self.str, self.dip, self.phi, self.mcc)
        #output message
        if self.Pc > Pc0:
            print('alert: critical fracture at Tau = %.2e, sn = %.2e, and Pp %.2e having phi = %.2e rad, mcc = %.2e Pa, Pc = %.2e Pa adjusted to phi = %.2e rad, mcc = %.2e Pa, Pc = %.2e Pa'
                  %(self.tau,self.sn,pres,phi0,mcc0,Pc0,self.phi,self.mcc,self.Pc))
            
    #set fracture cohesion to critical value at specified conditions
    def make_critical(self,pres=0.0):
        #shear stability limit
        mcc_c = self.tau - (self.sn - pres)*self.mu
        #tensile stability limit
        mcc_t = pres - self.sn
        #update fracture cohesion to enforce critical stability at the input conditions
        self.mcc = np.max([mcc_c,mcc_t])
        #recompute critical conditions
        self.Pc, self.sn, self.tau = self.tensor.Pc_frac(self.str, self.dip, self.phi, self.mcc)
        #output message
        print('         hydrofrac critical cohesion calculated as %.2e Pa to obtain Pc = %.2e Pa' %(self.mcc,self.Pc))
        #error case
        if self.Pc <= self.sn:
            print('         *** error: solution gives closed supercritical shear instead of hydrofracture')

def gen_hydrofracs(data,row,w):
    hydfs = []
    for i in range(0,len(w)):
        if w[i].type in ['procluster','perfcluster']:
            #hydraulic fractures
            num = data['Well_Clusters_count'][row]
            #location information
            spa = w[i].leg/(num+1)
            c0 = w[i].c0
            azn = w[i].azn
            dip = w[i].dip
            vAxi = w[i].vAxi
            #build list of fractures
            for n in range(0,num):
                #fracture parameters
                dia = 2.0*data['Stimulation_Radius_m'][row] #* (1.0 + stats.norm_trunc(1,0.0,0.2,-0.5,0.5)[0]*s.Strategy_SecondOrderNoise_ratio
                azn = np.pi/2.0+data['Stress_ShminAzimuth_rad'][row] #+ stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad*s.Strategy_SecondOrderNoise_ratio
                dip = np.pi/2.0-data['Stress_ShminDip_rad'][row] #+ stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad*s.Strategy_SecondOrderNoise_ratio
                c0 = c0 + spa*vAxi
                #add to model domain
                hydfs += [surf(c0,dia,azn,dip,'propped',data,row)]
                hydfs[-1].wi = w[i].pID
    #result
    return hydfs


#***************************************************************************************************************
###from: basic
#***************************************************************************************************************

#record information into fractures
f = fractures.gen_hydrofracs(s,w)
heterogeneity = np.linspace(0,1,int(s.Well_Clusters_count))**3*(1-s.Fracture_Tortuosity_ratio)+s.Fracture_Tortuosity_ratio
q_fracs = (qinj*heterogeneity**3.0)/np.sum(heterogeneity**3.0)
roll = -1
fdex = 0
bh_var = (12.0*q_fracs*0.5*(i_mu+p_mu)/(dPf/(2.0+1.0/(w_n/w_i))))**(1.0/3.0)
bd_var = q_fracs*0.5*(i_mu+p_mu)/(f_k*dPf/(2.0+1.0/(w_n/w_i)))
for a in range(0,len(f)):
    if f[a].type == 'propped':
        if (fdex == 0) or (fdex == s.Well_Clusters_count):
            roll = a
            fdex = 0
        fdex += 1
        f[a].bh = bh_var[s.Well_Clusters_count-1-(a-roll)]
        f[a].bd = bd_var[s.Well_Clusters_count-1-(a-roll)] #f[a].bd0p*np.exp(f_e*pn_guess[1])
        f[a].bd0p = f[a].bd/np.exp(f_e*pn_guess[1])
        f[a].bd0 = f[a].bd0p
        f[a].Pcen = 0.5*(Ps[2,-1]+Ps[1,-1])
        f[a].Pmax = Ps[1,-1]
        f[a].prop_load = 0.64*f[a].bd0p*np.pi*f_R**2.0
        f[a].stim = s.Stimulation_Steps_count
        f[a].fc = 3.28084*1.01324e15*(f[a].bh**3.0)/12.0



#***************************************************************************************************************
###from: geodt
#***************************************************************************************************************
#oriented circular fracture
class fracture:
    #initializations
    def __init__(self,xyz=[0,0,0],azn=0,dip=0,dia=1,nodes=9,surf=[]):
        if surf:
            xyz = surf.c0
            azn = surf.str
            dip = surf.dip
            dia = surf.dia
        if dip > np.pi/2:
            dip = np.pi-dip
            azn = azn+np.pi
        if dip < -np.pi/2:
            dip = np.pi+dip
            azn = azn+np.pi
        #build pointset
        ang = np.linspace(0.25,2.25,nodes)*np.pi
        d_r = ((np.sin(ang)*np.cos(-dip))**2+(np.cos(ang))**2)**0.5
        d_a = np.arccos(np.cos(ang)/d_r)+2*np.pi-(np.sign(np.sin(ang))<0)*np.pi
        d_s = np.sign(np.sin(ang))<0
        d_a = (1-d_s)*np.arccos(np.cos(ang)/d_r) + (d_s)*(2.0*np.pi-np.arccos(np.cos(ang)/d_r))
        self.pts = 0.5*dia*np.transpose(np.asarray(
            [d_r*np.sin(d_a+azn),
             d_r*np.cos(d_a+azn),
             np.sin(ang)*np.sin(-dip)]))
        #offset for center
        self.pts[:,0] += xyz[0]
        self.pts[:,1] += xyz[1]
        self.pts[:,2] += xyz[2]
        #meshgrid
        self.xs = np.concatenate(([self.pts[:,0]],[np.zeros(nodes)+xyz[0]]),axis=0)
        self.ys = np.concatenate(([self.pts[:,1]],[np.zeros(nodes)+xyz[1]]),axis=0)
        self.zs = np.concatenate(([self.pts[:,2]],[np.zeros(nodes)+xyz[2]]),axis=0)
#draw on the class' canvas, default translucent grey
def draw_fracture(self,fracture,color=(0.85, 0.85, 0.85, 0.1)):
    self.ax.plot_surface(fracture.xs,fracture.ys,fracture.zs,
                         rstride=1,cstride=1,color=color,linewidth=0)
class pipe:
    def __init__(self,c0=[0,0,0],azn=0,dip=0,leg=1,c1=[],well=[]):
        if well:
            c0 =  well.c0
            azn = well.azn
            dip = well.dip
            leg = well.leg
        elif len(c1)==3:
            pass
        else:
            vAxi = np.asarray([np.sin(azn)*np.cos(-dip), np.cos(azn)*np.cos(-dip), np.sin(-dip)])
            c1 = c0 + vAxi*leg
        self.xs = [c0[0], c1[0]]
        self.ys = [c0[1], c1[1]]
        self.zs = [c0[2], c1[2]]
def draw_pipe(self,well,color=(0.15, 0.15, 0.15, 0.8)):
    self.ax.plot3D(well.xs,well.ys,well.zs,color=color)
#demo draw
def demo(self):
    #place fractures
    self.draw_fracture(self.fracture([.3,.1,.5], np.pi*1.75, -np.pi*0.15, 1.0))
    self.draw_fracture(self.fracture([-.1,-.3,-.7], np.pi*.1, np.pi*0.45, 2.0))
    #place line
    self.draw_pipe(self.pipe([-0.6,0.1,0.5], 0.6*np.pi, 0.1*np.pi, 1.0))
    self.draw_pipe(self.pipe([-0.6,0.3,0.5], 0.5*np.pi, 0.1*np.pi, 1.2))
#geodt plotting tool
def show_geodt(self,geom,cmap='rainbow',**kwargs): #examples: faces='Pc', pipes='q', nodes='T'
    #get colormap
    self.cmap = getattr(matplotlib.cm, cmap)
    #set initial viewing size
    self.size = geom.setup.Domain_Size_m
    self.ax.set_xlim([-self.size,self.size])
    self.ax.set_ylim([-self.size,self.size])
    self.ax.set_zlim([-0.8*self.size,0.8*self.size])
    #TODO add hiding features
    #hide non-flowing fractures (ci>=0)
    #hide non-flowing pipes
    #plot by keywords
    for key, item in kwargs.items():
        #fractures and planes
        if key == 'faces':
            if len(geom.faces) > 0:
                if hasattr(geom.faces[0],item):
                    #print('%s has item %s' %(key, item))
                    #get color range
                    cs = [] 
                    for face in geom.faces:
                        cs += [getattr(face, item)]
                    crange = [np.min(cs),np.max(cs)]
                    #plot fractures
                    for i in range(6,len(geom.faces)):
                        #get color
                        color = self.cmap((cs[i]-crange[0])/(crange[1]-crange[0]),0.20)
                        #draw fracture
                        self.draw_fracture(self.fracture(surf=geom.faces[i]),color=color)