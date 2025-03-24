# -*- coding: utf-8 -*-
print('GeoDT_4.2.0')

# ****************************************************************************
# Calculate economic potential of EGS & optimize borehole layout with caging
# Author: Luke P. Frash
#
# Notation:
#  !!! or TODO for code needing to be updated (e.g, limititations or work in progress)
#  *** for breaks betweeen key sections
# ****************************************************************************

# ****************************************************************************
#### libraries
# ****************************************************************************
#internal
# if __package__ is None or __package__ == '':
#     from .libs.gt_linalg import solve
#     from .libs import gt_vtk as sg
#     from .libs import gt_properties
#     water = gt_properties.water()
from .libs.gt_linalg import solve
from .libs import gt_vtk as sg
from .libs import gt_properties as properties
water = properties.water()

#external
import numpy as np
import pylab
import math
from scipy import stats
import matplotlib.pyplot as plt
import pickle
import os

# ****************************************************************************
#### unit conversions
# ****************************************************************************
# Unit conversions for convenience
lpm=(0.1*0.1*0.1)/(60.0)
ft=12*25.4e-3#m
m=(1.0/ft)#ft
deg=1.0*math.pi/180.0
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
pi=math.pi


# ****************************************************************************
#### classes, functions, and modules
# ****************************************************************************
# Place a radial hydraulic fracture of radius r at x0
def HF(r,x0, strikeRad, dipRad, h=0.5):
    # start with a disk
    disk=sg.diskObj(r,h)
    disk=sg.rotateObj(disk,[0.0,1.0,0.0],dipRad)
    disk=sg.rotateObj(disk,[0.0,0.0,1.0],-strikeRad)
    disk=sg.transObj(disk,x0)
    return disk

#definitions and cross referencing for pipe types
def typ(key):
    ret = []
    choices = np.asarray([
            ['shunt','-5'],
            ['screen','-4'],
            ['perfcluster','-3'],
            ['producer', '-2'],
            ['injector', '-1'],
            ['pipe', '0'],
            ['fracture', '1'],
            ['propped', '2'],
            ['darcy', '3'],
            ['choke', '4'],
            ['perf', '5'],
            ['boundary', '6']
            ])
    key = str(key)
    ret = np.where(choices == key)
    if ret[1] == 0:
        ret = int(choices[ret[0],1][0])
    elif ret[1] == 1:
        ret = str(choices[ret[0],0][0])
    else:
        print( '**invalid pipe type defined**')
        ret = []
    return ret

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
        azn = pi
    else:
        azn = np.sign(dx)*np.arccos(dy/dr) + (1 - np.sign(dx))*pi
    if dr == 0.0:
        dip = -np.sign(dz)*pi/2.0
    else:
        dip = -np.arctan(dz/dr)
    return azn, dip, dis

def exponential_trunc(nsam,bval=1.0,Mmax=5.0,Mwin=1.0,prob=0.1): # random samples from a truncated exponential distribution
    # zero-centered Mcap
    lamba = np.log(10)*bval
    Mcap = (-1.0/lamba)*np.log(prob/(np.exp(lamba*Mwin)-1.0+prob))
    # calculated Mmin
    Mmin = Mmax-Mcap
    # sample sets with resamples out-of-range values
    s0 = np.random.exponential(1.0/(np.log(10)*bval), nsam) + Mmin
    iters = 0
    while 1:
        iters += 1
        r_pl = s0 > Mmax
        if np.sum(r_pl) > 0:
            r_dr = np.random.exponential(1.0/(np.log(10)*bval), nsam) + Mmin
            s0 = s0*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            break
    return s0

def contact_trunc(nsam,weight=0.15,
                  bd_nom=0.002,stddev=0.5*0.002,
                  exp_B=0.5,exp_C=0.5/np.pi,
                  bd_min=0.0001,bd_max=3.0*0.002): # get random samples from a 'contact' distribution
    # binomial samples
    n1 = np.random.binomial(nsam,weight,(1))
    n2 = nsam - n1
    # exponential samples
    s1 = np.random.exponential(exp_B, n1)*exp_C*bd_nom + bd_min
    iters = 0
    while 1:
        iters += 1
        r_pl = s1 > bd_max
        if np.sum(r_pl) > 0:
            r_dr = np.random.exponential(exp_B, n1)*exp_C*bd_nom + bd_min
            s1 = s1*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            break
    # normal samples
    s2 = np.random.normal(0,1,n2)*stddev + bd_nom
    iters = 0
    while 1:
        iters += 1
        r_pl = (s2 > bd_max) + (s2 < bd_min)
        if np.sum(r_pl) > 0:
            r_dr = np.random.normal(0,1,n2)*stddev + bd_nom
            s2 = s2*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            break
    s0 = np.concatenate((s1,s2),axis=0)
    return s0

def lognorm_trunc(nsam,logmu=0.0,logdev=1.0,
                  loglo=-2.0,loghi=2.0): # get random samples from a log normal distribution
    # normal samples
    s0 = np.random.normal(0,1,nsam)*logdev + logmu
    iters = 0
    while 1:
        iters += 1
        r_pl = (s0 > loghi) + (s0 < loglo)
        if np.sum(r_pl) > 0:
            r_dr = np.random.normal(0,1,nsam)*logdev + logmu
            s0 = s0*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            break
    return 10.0**s0

def norm_trunc(nsam,mu=0.0,dev=1.0,
                  lo=-2.0,hi=2.0): # get random samples from a log normal distribution
    s0 = np.random.normal(0,1,nsam)*dev + mu
    iters = 0
    while 1:
        iters += 1
        r_pl = (s0 > hi) + (s0 < lo)
        if np.sum(r_pl) > 0:
            r_dr = np.random.normal(0,1,nsam)*dev + mu
            s0 = s0*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            r_dr = np.random.uniform(0,1,nsam)*(hi-lo) + lo
            s0 = s0*(1-r_pl) + r_dr*(r_pl)
            break
    return s0

def save_csv(out=[['null',0.0]],filename='save.csv',append=True):
    #output to file
    out = list(map(list,zip(*out)))
    head = out[0][0]
    for i in range(1,len(out[0])):
        head = head + ',' + out[0][i]
    try:
        data = '%i' %(out[1][0])
    except:
        data = '' + out[1][0]
    for i in range(1,len(out[1])):
        try:
            data = data + ',%.5e' %(out[1][i])
        except:
            data = data + ',' + out[1][i]
    if append:
        try:
            with open(filename,'r') as f:
                test = f.readline()
            f.close()
            if test != '':
                with open(filename,'a') as f:
                    f.write(data + '\n')
                f.close()
            else:
                with open(filename,'a') as f:
                    f.write(head + '\n')
                    f.write(data + '\n')
                f.close()
        except:
            with open(filename,'a') as f:
                f.write(head + '\n')
                f.write(data + '\n')
            f.close()
    else:
        with open(filename,'w') as f:
            f.write(head + '\n')
            f.write(data + '\n')
        f.close()

def build_csv(header=[],data=[],filename='save.csv',append=True):
    print('data length: %i' %(len(data)))
    for i in range(0,len(data)):
        out = []
        for j in range(0,len(header)):
            out += [[header[j],data[header[j]][i]]]
            # if j == 0:
            #     print(data[header[j]][i])
        if i == 0:
            save_csv(out,filename,append)
        else:
            save_csv(out,filename,True)

def special_csv(header=[],data=[],source='',target='',filename='save.csv',append=True):
    for i in range(0,len(data)):
        out = []
        for j in range(0,len(header)):
            if (i == 1) and (header[j] == 'Strategy_SourceDirectory_path'):
                out += [[header[j],source]]
            elif (i == 1) and (header[j] == 'Strategy_TargetDirectory_path'):
                out += [[header[j],target]]
            else:
                out += [[header[j],data[header[j]][i]]]
        if i == 0:
            save_csv(out,filename,append)
        else:
            save_csv(out,filename,True)

def copy_csv(sourcefile='load.csv',targetfile='save.csv',append=True):
    #read sourcefile
    data = []
    try:
        with open(sourcefile,'r') as f1:
            data = f1.readlines()
            f1.close()
        if not(append):
            with open(targetfile,'w') as f2:
                f2.writelines(data)
                f2.close()
        else:
            with open(targetfile,'r') as f2:
                test = f2.readline()
                f2.close()
            with open(targetfile,'a') as f2:
                if test != '':
                    f2.writelines(data[1:])
                else:
                    f2.writelines(data[:])
                f2.close()
    except:
        print('failed to copy file %s to %s' %(sourcefile, targetfile)) 

def load_csv(filename='save.csv'):
    names = []
    data = []
    if False: #TODO: replace recfromcsv with custom importer (genfromtxt fails to recognize text)
        data = np.recfromcsv(filename,delimiter=',',filling_values=np.nan,deletechars='()',
                              case_sensitive=True,names=True,encoding=None)
        names = data.dtype.names
    elif True: #note that genfromtxt doesn't work for a variety of reasons, so a custom function is needed
        #load data as list of rows
        with open(filename,'r') as f:
            file = f.readlines()
            f.close()
        #convert into arrays
        names = file[0].split(',')
        names[-1] = names[-1][:-1]
        for i in range(1,len(file)):
            data += [file[i][:-1].split(',')]
        #convert into structured array
        numbers = np.ones(len(names))
        dtype = []
        for j in range(0,len(names)):
            try:
                for i in range(0,len(data)):
                    float(data[i][j])
                if (data[i][j].find('.') > -1) or (data[i][j] == "inf"):
                    numbers[j] = 2
            except:
                numbers[j] = 0
            if numbers[j] == 2:
                dtype += [tuple([names[j],'<f8'])]
            elif numbers[j] == 1:
                dtype += [tuple([names[j],'<i8'])]
            else:
                dtype += [tuple([names[j],'<U600'])]
        for i in range(0,len(data)):
            data[i] = tuple(data[i])
        dtype = np.dtype(dtype)
        data = np.array(data,dtype=dtype)
    else:
        print('%s not found in %s' %(filename, os.getcwd()))
    return names, data

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

class repickle(pickle.Unpickler):
    def find_class(self, module, name):
        remodule = module
        if (module == "solvers.geodt") or (module == "GeoDT"):
            remodule = "plugins.solvers.geodt"
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

def enthalpy_temperature_poly(T5=298.0,Tr=851.0,Pr=67.38e6,hu=136):
    #****** enthalpy function linearization *******
    x = np.linspace(T5,Tr,100,dtype=float)
    y = np.zeros(100,dtype=float)
    for i in range(0,len(x)):
        y[i] = water.h_from_PT(P=Pr,T=x[i]) #kJ/kg
    hTP = np.polyfit(x,y,3)
    ThP = np.polyfit(y,x,3)
    fig = pylab.figure(figsize=(8, 8), dpi=96, facecolor='w', edgecolor='k',tight_layout=True) # Medium resolution
    font = {'family' : 'serif',   'size'   : 12}
    pylab.rc('font', serif='Arial')
    pylab.rc('font', **font)
    #temperature
    ax = fig.add_subplot(111)
    ax.plot(x,y,'.',color='k') #points
    ax.plot(x,hTP[0]*x**3.0 + hTP[1]*x**2.0 + hTP[2]*x**1.0 + hTP[3],color='r')
    ax.plot(ThP[0]*y**3.0 + ThP[1]*y**2.0 + ThP[2]*y**1.0 + ThP[3],y,color='b')
    ax.plot(ThP[0]*hu**3.0 + ThP[1]*hu**2.0 + ThP[2]*hu**1.0 + ThP[3],hu,'x',color='m')
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Enthalpy (kJ/kg)')
    pylab.show()

#analytical flow rates and transient heat production (Gringarten et al., 1975: Journal of Geophysical Research)
class gringarten:
    def __init__(self):
        self.solution = []

class cauchy: #functions modified from JPM
    def __init__(self):
        self.sigP = np.zeros((3,3),dtype=float)
        self.sigG = np.zeros((3,3),dtype=float)
        self.Sh = 0.0 #Pa
        self.SH = 0.0 #Pa
        self.SV = 0.0 #Pa
        self.Sh_azn = 0.0*deg #rad
        self.Sh_dip = 0.0*deg #rad
    # http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    def rotationMatrix(self, axis, theta):
        #return the rotation matrix associated with counterclockwise rotation about the given axis by theta radians.
        axis = np.asarray(axis)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
    #http://www.continuummechanics.org/stressxforms.html
    def rotateTensor(self, tensor, axis, theta):
        rot=self.rotationMatrix(axis,theta)
        return np.dot(rot,np.dot(tensor,np.transpose(rot)))
    #Projection of normal plane from dip directions
    def normal_from_dip(self, dip_direction,dip_angle):
            # dip_direction=0 is north -> nrmxH=0, nrmyH=1
            # dip_direction=90 is east -> nrmxH=1, nrmyH=0
            nrmxH=np.sin(dip_direction)
            nrmyH=np.cos(dip_direction)
            # The lateral components of the normals are corrected for the dip angle
            nrmx=nrmxH*np.sin(dip_angle)
            nrmy=nrmyH*np.sin(dip_angle)
            # The vertical
            nrmz=np.cos(dip_angle)
            return np.asarray([nrmx,nrmy,nrmz])
    #normal stress and critical slip or opening pressure
    def Pc(self, nrmP, phi, mcc):
        # Shear traction on the fault segment
        t=np.zeros([3])
        for i in range(3):
            for j in range(3):
                t[i] += self.sigG[j][i]*nrmP[j]
        # Normal component of the traction
        Sn=0.0
        for i in range(3):
            Sn += t[i]*nrmP[i]
        # Shear component of the traction
        tauV=np.zeros([3])
        for i in range(3):
            tauV[i]=t[i]-Sn*nrmP[i]
        tau=np.sqrt(tauV[0]*tauV[0]+tauV[1]*tauV[1]+tauV[2]*tauV[2])
        # Critical pressure for slip from mohr-coulomb
        Pc1 = Sn - (tau-mcc)/np.tan(phi)
        # Critical pressure for tensile opening
        Pc2 = Sn + mcc
        # Critical pressure for fracture activation
        Pc = np.min([Pc1,Pc2])
        return Pc, Sn, tau
    #critical slip given fracture strike and dip
    def Pc_frac(self, strike, dip, phi, mcc):
        #get fracture normal vector
        nrmG=self.normal_from_dip(strike+np.pi/2, dip)
        return self.Pc(nrmG, phi, mcc)
    #Set cauchy stress tensor from rotated principal stresses
    def set_sigG_from_Principal(self,Sh,SH,SV,ShAzn,ShDip):
        # Stresses in the principal stress directions
        # We have been given the azimuth of the minimum stress, ShAzimuthDeg, to compare with 90 (x-dir)
        # That is Sh==Sxx, SH=Syy, SV=Szz in the principal stress coord system
        sigP=np.asarray([
            [Sh,  0.0, 0.0],
            [0.0, SH,  0.0],
            [0.0, 0.0, SV ]
            ])
        # Rotate about z-axis
        deltaShAz=ShAzn-np.pi/2.0 #(ShAznDeg-90.0)*np.pi/180.0
        # Rotate about y-axis
        ShDip=-ShDip #-ShDipDeg*np.pi/180.0
        sigG=self.rotateTensor(sigP,[0.0,1.0,0.0],-ShDip)
        sigG=self.rotateTensor(sigG,[0.0,0.0,1.0],-deltaShAz)
        self.sigP = sigP
        self.sigG = sigG
        self.Sh = Sh #Pa
        self.SH = SH #Pa
        self.SV = SV #Pa
        self.Sh_azn = ShAzn#Deg*deg #rad
        self.Sh_dip = ShDip#Deg*deg #rad
        return sigG
    #Plot citical pressure
    def plot_Pc(self, phi, mcc, filename='Pc_stereoplot.png'):
        # Working variables
        nRad=100
        dip_angle_deg = np.asarray(range(nRad+1))*(90.0/nRad)
        nTheta=200
        dip_dir_radians = np.asarray(range(nTheta+1))*(2.0*np.pi/nTheta)
        png_dpi=128
        # Calculate critical dP
        criticalDelPpG=np.zeros([nRad,nTheta])
        for i in range(nRad):
            for j in range(nTheta):
                # Convert the x and y into a normal to the fracture using dip angle and dip direction
                nrmG=self.normal_from_dip(dip_dir_radians[j], dip_angle_deg[i]*deg)
                nrmx,nrmy,nrmz=nrmG[0],nrmG[1],nrmG[2]
                criticalDelPpG[i,j], h1, h2 = self.Pc(np.asarray([nrmx,nrmy,nrmz]), phi, mcc)
        # Plot critical slip pressure (lower hemisphere projection)
        fig = pylab.figure(figsize=(6,4.75),dpi=png_dpi,tight_layout=True,facecolor='w',edgecolor='k')
        ax = fig.add_subplot(111, projection='polar')
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        pylab.pcolormesh(dip_dir_radians,dip_angle_deg,np.ma.masked_where(np.isnan(criticalDelPpG),criticalDelPpG), 
                        #vmin=110.0*MPa, vmax=150.0*MPa,cmap='rainbow_r') #vmin=0.0*MPa, vmax=self.SH,cmap='rainbow_r')
                        vmin=self.Sh*0.5, vmax=np.max([self.SV,self.SH]),cmap='rainbow_r')
        ax.grid(True)
        ax.set_rgrids([0,30,60,90],labels=[])
        ax.set_thetagrids([0,90,180,270])
        pylab.colorbar()
        ax.set_title("Critical pressure for Mohr-Coulomb failure", va='bottom')
        pylab.savefig(filename, format='png',dpi=png_dpi)
        pylab.close()
    #convert vector to strike and dip
    def vnorm_to_strdip(self,vnorm=np.asarray([0.0,0.0,0.0])):
        dr = (vnorm[0]**2.0+vnorm[1]**2.0)**0.5
        azn = []
        dip = []
        if vnorm[0] == 0 and vnorm[1] >= 0:
            azn = 0.0
        elif vnorm[0] == 0:
            azn = np.pi
        else:
            azn = np.sign(vnorm[0])*np.arccos(vnorm[1]/dr) + (1 - np.sign(vnorm[0]))*np.pi
        if dr == 0.0:
            dip = -np.sign(vnorm[2])*np.pi/2.0
        else:
            dip = -np.arctan(vnorm[2]/dr)
        if dip < 0:
            dip = -dip
            azn = azn + np.pi
        azn = azn + np.pi/2
        dip = np.pi/2 - dip
        if azn >= 2.0*np.pi:
            azn = azn - 2.0*np.pi
        elif azn < 0:
            azn = azn + 2.0*np.pi
        return [azn, dip]
    #critical orientations
    def get_conjugates(self,plots=False):
        # Working variables
        phi = 30.0*deg
        mcc = 1.0*MPa
        nRad=18
        dip_angle_deg = np.asarray(range(nRad+1))*(90.0/nRad)
        nTheta=72
        dip_dir_radians = np.asarray(range(nTheta+1))*(2.0*np.pi/nTheta)
        # Calculate critical dP3
        dP3 = [1e10,0,0] #global minimum
        criticalDelPpG=np.zeros([nRad,nTheta])
        for i in range(nRad):
            for j in range(nTheta):
                # Convert the x and y into a normal to the fracture using dip angle and dip direction
                nrmG=self.normal_from_dip(dip_dir_radians[j], dip_angle_deg[i]*deg)
                nrmx,nrmy,nrmz=nrmG[0],nrmG[1],nrmG[2]
                criticalDelPpG[i,j], h1, h2 = self.Pc(np.asarray([nrmx,nrmy,nrmz]), phi, mcc)
                if dP3[0] > criticalDelPpG[i,j]:
                    dP3 = [criticalDelPpG[i,j], i, j]
        azn = dip_dir_radians[dP3[2]]
        dip = np.pi/2 - dip_angle_deg[dP3[1]]*deg
        v_dP3 = np.asarray([np.sin(azn)*np.cos(-dip), np.cos(azn)*np.cos(-dip), np.sin(-dip)])            
        # Rotate about shmin to get dP2
        v_sh = np.asarray([np.sin(self.Sh_azn)*np.cos(-self.Sh_dip), np.cos(self.Sh_azn)*np.cos(-self.Sh_dip), np.sin(-self.Sh_dip)])
        v_dP2 = np.dot(self.rotationMatrix(v_sh, np.pi),v_dP3)
        # Calculate dP1 using cross product of dP
        v_dP1 = np.cross(v_dP3,v_dP2)
        v_dP1 = v_dP1/np.linalg.norm(v_dP1)
        # Convert to strike and dip
        str_dip = np.asarray([self.vnorm_to_strdip(v_dP3),
                              self.vnorm_to_strdip(v_dP2),
                              self.vnorm_to_strdip(v_dP1)])
        # Plot critical slip pressure (lower hemisphere projection)
        if plots:
            fig = pylab.figure(figsize=(6,4.75),dpi=128,tight_layout=True,facecolor='w',edgecolor='k')
            ax = fig.add_subplot(111, projection='polar')
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            pylab.pcolormesh(dip_dir_radians,dip_angle_deg,np.ma.masked_where(np.isnan(criticalDelPpG),criticalDelPpG), 
                            #vmin=110.0*MPa, vmax=150.0*MPa,cmap='rainbow_r') #vmin=0.0*MPa, vmax=self.SH,cmap='rainbow_r')
                            vmin=self.Sh*0.5, vmax=np.max([self.SV,self.SH]),cmap='rainbow_r')
            ax.plot(dip_dir_radians[dP3[2]],dip_angle_deg[dP3[1]],'om')
            ax.plot(str_dip[0,0]-np.pi/2,str_dip[0,1]/deg,'or')
            ax.plot(str_dip[:,0]-np.pi/2,str_dip[:,1]/deg,'.k')
            ax.grid(True)
            ax.set_rgrids([0,30,60,90],labels=[])
            ax.set_thetagrids([0,90,180,270])
            pylab.colorbar()
            ax.set_title("Critical pressure for Mohr-Coulomb failure", va='bottom')
            pylab.savefig('check.png', format='png',dpi=128)
            pylab.close()
        return str_dip
        
#model definition
class setup:
    def __init__(self):
        self.Strategy_SourceDirectory_path = os.getcwd()
        self.Strategy_TargetDirectory_path = os.getcwd()
        self.Strategy_DiscreteFractures_path = os.getcwd()+'\\dfn.csv'
        self.Strategy_Iterations_units = 10
        self.Strategy_Design_type = 'ALL' #'EGS' 'AGS' 'CGS' 'ALL' 'O&G'
        self.Strategy_Target_type = 'Deep' #'Temp'
        self.Strategy_Stim_type = 'Volu' #'Radi'
        self.Well_TargetTemp_K = 273.15+325 #K
        self.Well_TargetDepth_m = 6000.0 #m
        self.Well_ProducerCount_wells = 4 #wells
        #self.Well_InjectionIntervals_intervals = 1 #breaks in well length
        self.Well_Spacing_m = 400.0 #m
        self.Well_DeviatedLength_m = 800.0 #m
        self.Well_ProducerProportion_ratio = 0.6 #m/m
        self.Well_Azimuth_rad = 90.0*deg #rad
        self.Well_Dip_rad = 30.0*deg #rad
        self.Well_RotationPhase_rad = 0.0*deg #rad
        self.Well_RotationToe_rad = 0.0*deg #rad
        self.Well_RotationSkew_rad = 0.0*deg #rad
        self.Well_HydraulicRadius_m = 0.0254*5.0 #m
        self.Well_CasingRadius_m = 0.0254*5.5 # m
        self.Well_BoreRadius_m = 0.0254*6.5 # m
        self.Well_Roughness_metric = 80.0
        self.Circulation_InjectionRate_m3ps = 0.10 #m3/s
        self.Stimulation_InjectionRate_m3ps = 0.4 #m3/s
        self.Stimulation_TargetVolume_m3 = 50000.0 #m3
        self.Stimulation_TargetRadius_m = 525.0 #m3
        self.Domain_Size_m = 1000.0 #m
        self.Rock_ThermalGradient_Kpm = 0.060 #K/m
        self.Rock_Density_kgpm3 = 2700.0 # kg/m3
        self.Rock_ThermalConductivity_WpmK = 2.5 # W/m-K
        self.Rock_HeatCapacity_kJpm3K = 2063.0 # kJ/m3-K
        self.Air_Temperature_K = 273.15 # K
        self.Air_Pressure_Pa = 0.101*MPa #Pa
        self.Rock_YoungsModulus_Pa = 50.0*GPa #Pa
        self.Rock_PoissonRatio_ratio = 0.3 #
        self.StressMin_Coefficient_ratio = 0.5 #Pa/Pa
        self.StressMax_Coefficient_ratio = 0.75 #Pa/Pa
        self.StressMin_Azimuth_rad = 90.0*deg #rad
        self.StressMin_Dip_rad = 0.0*deg #rad
        self.StressMin_Variance_rad = 0.5*deg #rad
        self.Fracture_SlipLengthMin_ratio = 10.0**-3.0 #m/m
        self.Fracture_SlipLengthNom_ratio = 10.0**-2.0 #m/m
        self.Fracture_SlipLengthMax_ratio = 10.0**-1.2 #m/m
        self.Fracture_DilationSlipMin_ratio = 0.000 #m/m
        self.Fracture_DilationSlipNom_ratio = 0.200 #m/m
        self.Fracture_DilationSlipMax_ratio = 0.800 #m/m
        self.Fracture_HydraulicDilationMin_ratio = 0.0 #m/m
        self.Fracture_HydraulicDilationNom_ratio = 0.6 #m/m
        self.Fracture_HydraulicDilationMax_ratio = 2.0 #m/m
        self.Fracture_CompressibilityMin_1pPa = 2.0e-9 #1/Pa
        self.Fracture_CompressibilityNom_1pPa = 2.9e-8 #1/Pa
        self.Fracture_CompressibilityMax_1pPa = 10.0e-8 #1/Pa
        self.Proppant_CompressibilityMin_1pPa = 2.0e-9 #1/Pa
        self.Proppant_CompressibilityNom_1pPa = 2.9e-8 #1/Pa
        self.Proppant_CompressibilityMax_1pPa = 10.0e-8 #1/Pa
        self.Fracture_InitialHydraulicApertureMin_m = 0.00005 #m
        self.Fracture_InitialHydraulicApertureNom_m = 0.00010 #m 
        self.Fracture_InitialHydraulicApertureMax_m = 0.00020 #m 
        self.Domain_BoundaryAperture_m = 0.003 #m
        self.Fracture_TortuosityMin_ratio = 0.80 #m/m
        self.Fracture_TortuosityNom_ratio = 0.90 #m/m
        self.Fracture_TortuosityMax_ratio = 1.00 #m/m
        self.Cement_ThermalConductivity_WpmK = 2.0 # W/m-K
        self.Cement_HeatCapacity_kJpm3K = 2000.0 # kJ/m3-K
        self.Power_GeneralEfficiency_ratio = 0.85 # kWe/kWt
        self.Power_LifeSpan_s = 20.5*yr #years
        self.Domain_TimeSteps_steps = 41 #steps
        self.Circulation_Backpressure_Pa = 1.0*MPa #Pa
        self.Circulation_InjectionTemperature_K = 95.0+273.15 #K
        self.Power_ConvectionCoefficient_kWpm2K = 3.0 #kW/m2-K
        self.Fluid_NominalDensity_kgpm3 = 965.0 #kg/m3
        self.Fluid_Viscosity_Pas = 0.2*cP #Pa-s
        self.Rock_Permeability_m2 = 0.1*mD #m2
        self.Proppant_InitialPermeabilityMin_m2 = 10.0*darcy #m2
        self.Proppant_InitialPermeabilityNom_m2 = 100.0*darcy #m2
        self.Proppant_InitialPermeabilityMax_m2 = 300.0*darcy #m2
        self.JointSets_Set1_fractures = 16 #fractures
        self.JointSets_Set2_fractures = 8 #fractures
        self.JointSets_Set3_fractures = 4 #fractures
        # self.JointSets_Set4_fractures = 0 #fractures
        # self.JointSets_Set5_fractures = 0 #fractures
        self.JointSets_Variance_rad = 7.0*deg #rad
        self.JointSets_DiameterMin_m = 10.0 #m
        self.JointSets_DiameterMax_m = 1000.0 #m
        self.Stimulation_Clusters_clusters = 3 #clusters
        self.Stimulation_PerfDiameter_m = 0.013 #m
        self.Stimulation_PerfPerCluster_perfs = 6 #holes
        self.Stimulation_ProppantConcentration_m3pm3 = 0.045 #m3/m3
        self.Domain_PressureIncrement_Pa = 0.1*MPa #Pa
        self.Domain_StimulationLimit_steps = 1
        self.Circulation_AllowableOverpressure_Pa = 0.99 #pinj/s3 #Pa
        self.Shear_FrictionAngleMin_rad = 20.0*deg #rad
        self.Shear_FrictionAngleNom_rad = 35.0*deg #rad
        self.Shear_FrictionAngleMax_rad = 45.0*deg #rad
        self.Shear_CohesionMin_Pa = 1.0*MPa #Pa
        self.Shear_CohesionNom_Pa = 7.0*MPa #Pa
        self.Shear_CohesionMax_Pa = 15.0*MPa #Pa
        self.Tension_Cohesion_Pa = 0.1*MPa #Pa
        self.Tension_FrictionAngle_rad = 30.0*deg #rad
        self.Tension_Toughness_Pasqrtm = 1.5*MPa #Pa-m**0.5
        self.Economics_ElectricitySales_USDpkWh = 0.1372 #$/kWh - customer electricity retail price                  
        self.Economics_DrillingCost_USDpm = 2763.06 #$/m - Lowry et al, 2017 large diameter well baseline
        self.Economics_AccessCost_USD = 590e3 #$ Lowry et al, 2017 large diameter well baseline
        self.Economics_EquipmentCost_USDpkW = 2025.65 #$/kWe simplified from GETEM model 
        self.Economics_ExplorationCost_USDpm = 2683.41 #$/m simplified from GETEM model
        self.Economics_Maintenance_USDpkWh = 0.03648 #$/kWh simplified from GETEM model
        self.Economics_Seismic_USDpMw = 2e-4 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re)
        self.Economics_Seismic_exp = 5.0 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re)

class payoff:
    def __init__(self):
        self.Power_Net_kW = 0.0
        self.Economics_NetPresentValue_USD = 0.0
        self.LCOE_Net_USDpMWh = 0.0
        self.Injection_Pressure_MPa = 0.0
        self.Injection_Rate_m3ps = 0.0
        
        self.Power_Bulk_kW = 0.0
        self.Power_Pump_kW = 0.0
        self.Power_Thermal_kW = 0.0
        self.Injection_Rate_kgps = 0.0
        self.Injection_Temperature_K = 0.0
        self.Injection_Enthalpy_kJpkg = 0.0
        self.Injection_PressureLimit_MPa = 0.0
        self.Production_Pressure_MPa = 0.0
        self.Production_Rate_kgps = 0.0
        self.Production_Temperature_K = 0.0
        self.Production_Enthalpy_kJpkg = 0.0
        self.Economics_Capital_USD = 0.0
        self.Economics_Upkeep_USD = 0.0
        self.Economics_QuakeRisk_USD = 0.0
        self.Economics_PowerProfit_USD = 0.0
        self.Economics_BulkSales_USD = 0.0
        self.LCOE_Bulk_USDpMWh = 0.0
        self.LCOE_Thermal_USDpMWh = 0.0
        self.Drilling_TotalLength_m = 0.0
        self.Drilling_BottomTemperature_K = 0.0
        self.Drilling_BottomPressure_MPa = 0.0
        self.Drilling_WellDiameter_m = 0.0
        self.Drilling_CasingWeight_kgpm = 0.0
        self.Drilling_CementVolume_m3 = 0.0
        self.Stimulation_SandMass_kg = 0.0
        self.Stimulation_WaterVolume_m3 = 0.0
        self.Stimulation_PumpPressure_MPa = 0.0
        self.Circulation_WorkingVolume_m3 = 0.0
        self.Circulation_Recovery_pct = 0.0
        self.Circulation_MakeupWater_m3 = 0.0
        self.Circulation_PumpPressure_MPa = 0.0
        self.Stress_s3_MPa = 0.0
        self.Stress_s2_MPa = 0.0
        self.Stress_s1_MPa = 0.0
        self.Stress_pp_MPa = 0.0
        self.Seismicity_MaxQuake_Mw = 0.0
        self.Seismicity_NumQuake_ea = 0.0
        self.Fractures_ShearStim_ea = 0.0
        self.Fractures_HydroStim_ea = 0.0
        self.Fractures_Hydroprop_ea = 0.0

class reservoir:
    def re_init(self):
        #pick a strategy
        if self.strategy == 'ALL':
            pick = np.random.uniform(0,1)
            if pick < 0.333:
                self.strategy = 'AGS'
            elif pick < 0.666:
                self.strategy = 'EGS'
            else:
                self.strategy = 'CGS'
            self.setup.Strategy_Design_type = self.strategy
        elif not(self.strategy in ['AGS','EGS','CGS','O&G']):
            self.strategy = 'EGS'
            self.setup.Strategy_Design_type = 'EGS'
            
        #clusters and intervals
        if self.strategy == 'AGS':
            self.w_intervals = 1
            self.perf_clusters = 1
        elif self.strategy == 'CGS':
            self.w_intervals = max([1,int(self.setup.Stimulation_Clusters_clusters)])
            self.perf_clusters = 1
        else: # 'EGS' 'O&G'
            self.w_intervals = 1
            self.perf_clusters = max([1,int(self.setup.Stimulation_Clusters_clusters)]) 
        
        #well depth
        if self.setup.Strategy_Target_type.upper() == 'TEMP':
            self.ResDepth = (self.setup.Well_TargetTemp_K-self.AmbTempK)/self.ResGradient
            self.setup.Well_TargetDepth_m = self.ResDepth
        else: #'Deep'
            self.ResDepth = self.setup.Well_TargetDepth_m
            self.setup.Well_TargetTemp_K = self.ResDepth*self.ResGradient+self.AmbTempK #K
        
        #standard parameters
        self.ResG = self.ResE/(2.0*(1.0+self.Resv))
        self.BH_T = self.ResDepth*self.ResGradient + self.AmbTempK #K
        self.BH_P = self.PoreRho*g*self.ResDepth + self.AmbPres #Pa
        self.s1 = self.ResRho*g*self.ResDepth #Pa
        self.s2 = self.Ks2*(self.s1-self.BH_P)+self.BH_P #Pa
        self.s3 = self.Ks3*(self.s1-self.BH_P)+self.BH_P #Pa
        self.stress = cauchy()
        self.stress.set_sigG_from_Principal(self.s3, self.s2, self.s1, self.s3Azn, self.s3Dip)
        str_dip = self.stress.get_conjugates(plots=False)
        self.fStr = np.asarray([[str_dip[0,0],self.var],
                                [str_dip[1,0],self.var],
                                [str_dip[1,0],self.var]])
                                # [str_dip[2,0],self.var],
                                # [0,0],[0,0]],dtype=float) #m
        self.fDip = np.asarray([[str_dip[0,1],self.var],
                                [str_dip[1,1],self.var],
                                [str_dip[2,1],self.var]])
                                # [str_dip[2,1],self.var],
                                # [0,0],[0,0]],dtype=float) #m
        self.r_perf = 0.2*self.w_spacing #m
        self.Qinj = self.setup.Circulation_InjectionRate_m3ps/self.w_intervals
        self.Qstim = self.setup.Stimulation_InjectionRate_m3ps/self.w_intervals
        self.Vinj = self.Qinj*self.LifeSpan
        
        #injection pressure limit for circulation
        if self.strategy == 'EGS':
            self.Plim = np.min([0.0,self.Plim])
            self.pfinal_max = self.s3 + self.Plim #Pa
        elif self.strategy in ['CGS','XGS']:
            self.pfinal_max = 140e6 + self.BH_P # + self.p_whp
            self.Plim = self.pfinal_max-self.s3
        else: #['AGS']
            self.pfinal_max = self.s3 + self.Plim #Pa
           
        #stimulation volume - this frac model is clunky, but the goal is heat not hydrofracs
        if self.stimtype in ['RAD','R','RADI','RADIUS']:
            stim_trad = self.setup.Stimulation_TargetRadius_m
            stim_mult = np.max([0.9,stim_trad/self.r_perf])
            stim_step = np.log(stim_mult/0.5787)/0.1823
            stim_pres = self.dPi*stim_step
            stim_aper = 8.0*stim_pres*(1.0-self.Resv**2.0)*stim_trad/(np.pi*self.ResE)
            stim_volu = (4.0/3.0)*np.pi*stim_trad**2.0*0.5*stim_aper
            self.Vstim = stim_volu*self.setup.Stimulation_Clusters_clusters
            self.setup.Stimulation_TargetVolume_m3 = self.Vstim*self.w_intervals
        else: # ['V','VOL','VOLUME']
            self.Vstim = self.setup.Stimulation_TargetVolume_m3/self.w_intervals
            stim_volu = np.max([1e-9,self.Vstim/self.setup.Stimulation_Clusters_clusters])
            stim_vol0 = (4.0/3.0)*np.pi*self.r_perf**2.0*0.5*8.0*2.0*self.dPi*(1.0-self.Resv**2.0)*self.r_perf/(np.pi*self.ResE)
            stim_step = 1.6384*np.log(stim_volu/stim_vol0)+1.4738
            self.setup.Stimulation_TargetRadius_m = self.r_perf*0.5787*np.exp(0.1823*stim_step)
            
    def fetch(self,s):
        #--- displayed variables ---
        self.setup = s
        self.strategy = s.Strategy_Design_type.upper()
        self.stimtype = s.Strategy_Stim_type.upper()
        self.size = s.Domain_Size_m
        self.ResGradient = s.Rock_ThermalGradient_Kpm
        self.ResRho = s.Rock_Density_kgpm3
        self.ResKt = s.Rock_ThermalConductivity_WpmK
        self.ResSv = s.Rock_HeatCapacity_kJpm3K
        self.AmbTempK = s.Air_Temperature_K
        self.AmbPres = s.Air_Pressure_Pa
        self.ResE = s.Rock_YoungsModulus_Pa
        self.Resv = s.Rock_PoissonRatio_ratio
        self.Ks3 = s.StressMin_Coefficient_ratio
        self.Ks2 = s.StressMax_Coefficient_ratio
        self.s3Azn = s.StressMin_Azimuth_rad
        self.s3Dip = s.StressMin_Dip_rad
        self.s3AznVar = s.StressMin_Variance_rad
        self.s3DipVar = s.StressMin_Variance_rad
        self.gamma = np.asarray([s.Fracture_SlipLengthMin_ratio,
                                 s.Fracture_SlipLengthNom_ratio,
                                 s.Fracture_SlipLengthMax_ratio])
        self.a = np.asarray([s.Fracture_DilationSlipMin_ratio,
                             s.Fracture_DilationSlipNom_ratio,
                             s.Fracture_DilationSlipMax_ratio])
        self.N = np.asarray([s.Fracture_HydraulicDilationMin_ratio,
                             s.Fracture_HydraulicDilationNom_ratio,
                             s.Fracture_HydraulicDilationMax_ratio])
        self.alpha = np.asarray([s.Fracture_CompressibilityMin_1pPa,
                                 s.Fracture_CompressibilityNom_1pPa,
                                 s.Fracture_CompressibilityMax_1pPa])
        self.prop_alpha = np.asarray([s.Proppant_CompressibilityMin_1pPa,
                                      s.Proppant_CompressibilityNom_1pPa,
                                      s.Proppant_CompressibilityMax_1pPa])
        self.bh = np.asarray([s.Fracture_InitialHydraulicApertureMin_m,
                              s.Fracture_InitialHydraulicApertureNom_m,
                              s.Fracture_InitialHydraulicApertureMax_m])
        self.bh_bound = s.Domain_BoundaryAperture_m
        self.f_roughness = np.asarray([s.Fracture_TortuosityMin_ratio,
                                       s.Fracture_TortuosityNom_ratio,
                                       s.Fracture_TortuosityMax_ratio])
        self.w_count = int(s.Well_ProducerCount_wells)
        self.w_spacing = s.Well_Spacing_m
        self.w_length = s.Well_DeviatedLength_m
        self.w_azimuth = s.Well_Azimuth_rad
        self.w_dip = s.Well_Dip_rad
        self.w_proportion = s.Well_ProducerProportion_ratio
        self.w_phase = s.Well_RotationPhase_rad
        self.w_toe = s.Well_RotationToe_rad
        self.w_skew = s.Well_RotationSkew_rad
        # self.w_intervals = max([1,int(s.Well_InjectionIntervals_intervals)])
        self.ra = s.Well_HydraulicRadius_m
        self.rb = s.Well_CasingRadius_m
        self.rc = s.Well_BoreRadius_m
        self.rgh = s.Well_Roughness_metric
        self.CemKt = s.Cement_ThermalConductivity_WpmK
        self.CemSv = s.Cement_HeatCapacity_kJpm3K
        self.GenEfficiency = s.Power_GeneralEfficiency_ratio
        self.LifeSpan = s.Power_LifeSpan_s
        self.TimeSteps = int(s.Domain_TimeSteps_steps)
        self.p_whp = s.Circulation_Backpressure_Pa
        self.TinjK = s.Circulation_InjectionTemperature_K
        self.H_ConvCoef = s.Power_ConvectionCoefficient_kWpm2K
        self.PoreRho = s.Fluid_NominalDensity_kgpm3 #!!! this ought to be computed
        self.Poremu = s.Fluid_Viscosity_Pas #!!! this ought to be computed
        self.Porek = s.Rock_Permeability_m2 #!!! this isn't used yet
        self.kf = np.asarray([s.Proppant_InitialPermeabilityMin_m2,
                              s.Proppant_InitialPermeabilityNom_m2,
                              s.Proppant_InitialPermeabilityMax_m2])
        self.fNum = np.asarray([int(s.JointSets_Set1_fractures),
                                int(s.JointSets_Set2_fractures),
                                int(s.JointSets_Set3_fractures)])
                                # int(s.JointSets_Set3_fractures),
                                # int(s.JointSets_Set4_fractures),
                                # int(s.JointSets_Set5_fractures)])
        self.fDia = np.asarray([[s.JointSets_DiameterMin_m,s.JointSets_DiameterMax_m],
                                [s.JointSets_DiameterMin_m,s.JointSets_DiameterMax_m],
                                [s.JointSets_DiameterMin_m,s.JointSets_DiameterMax_m]])
                                # [s.JointSets_DiameterMin_m,s.JointSets_DiameterMax_m],
                                # [s.JointSets_DiameterMin_m,s.JointSets_DiameterMax_m],
                                # [s.JointSets_DiameterMin_m,s.JointSets_DiameterMax_m]])
        self.var = s.JointSets_Variance_rad
        # self.perf_clusters = int(s.Stimulation_Clusters_clusters)
        self.perf_dia = s.Stimulation_PerfDiameter_m
        self.perf_per_cluster = s.Stimulation_PerfPerCluster_perfs
        self.sand = s.Stimulation_ProppantConcentration_m3pm3
        self.dPi = s.Domain_PressureIncrement_Pa
        self.stim_limit = int(s.Domain_StimulationLimit_steps)
        self.Plim = s.Circulation_AllowableOverpressure_Pa
        self.phi = np.asarray([s.Shear_FrictionAngleMin_rad,
                               s.Shear_FrictionAngleNom_rad,
                               s.Shear_FrictionAngleMax_rad])
        self.mcc = np.asarray([s.Shear_CohesionMin_Pa,
                               s.Shear_CohesionNom_Pa,
                               s.Shear_CohesionMax_Pa]) #Pa
        self.hfmcc = s.Tension_Cohesion_Pa
        self.hfphi = s.Tension_FrictionAngle_rad
        self.Kic = s.Tension_Toughness_Pasqrtm
        self.sales_kWh = s.Economics_ElectricitySales_USDpkWh              
        self.drill_m = s.Economics_DrillingCost_USDpm
        self.pad_fixed = s.Economics_AccessCost_USD
        self.plant_kWe = s.Economics_EquipmentCost_USDpkW
        self.explore_m = s.Economics_ExplorationCost_USDpm
        self.oper_kWh = s.Economics_Maintenance_USDpkWh
        self.quake_coef = s.Economics_Seismic_USDpMw
        self.quake_exp = s.Economics_Seismic_exp
        
        #--- hidden variables ---
        self.dT0 = 10.0 #K
        self.bval = 1.0 #Gutenberg-Richter magnitude scaling
        self.gradient = True
        
        #--- calculated parameters
        self.re_init()

    def __init__(self):
        self.fetch(setup())
        
#surface object
class surf:
    def __init__(self,x0=0.0, y0=0.0, z0=0.0, dia=1.0, stk=0.0*deg, dip=90.0*deg,
                 ty='fracture', rock = reservoir(),
                 mcc = -1, phi = -1):
        #*** base parameters ***
        #node number of center point
        self.ci = -1
        #geometry
        self.c0 = np.asarray([x0, y0, z0])
        self.dia = dia
        self.str = stk
        self.dip = dip
        self.typ = typ(ty)
        #shear strength
        self.phi = -1.0 #rad
        self.mcc = -1.0 #Pa
        #stress state
        self.sn = 5.0*MPa
        self.En = 50.0*GPa
        self.vn = 0.30
        self.Pc = 0.0*MPa
        self.tau = 0.0*MPa
        #stimulation information
        self.stim = 0
        self.Pmax = 0.0*MPa
        self.Pcen = 0.0*MPa
        self.Mws = [] #[-10.0] #maximum magnitude seismic event tracker
        self.arup = 1.0 #rupture area available for seismicity
        self.hydroprop = False
        self.prop_load = 0.0 #m3 #absolute proppant volume
        self.prop_alpha = norm_trunc(1,rock.prop_alpha[1],rock.prop_alpha[1],rock.prop_alpha[0],rock.prop_alpha[2])[0] #proppant compressibility modulus
        self.roughness = np.random.uniform(rock.f_roughness[0],rock.f_roughness[2],(1))[0] #open flow roughness
        # self.kf = rock.kf #proppant pack permeability
        self.kf = np.random.uniform(rock.kf[0],rock.kf[2],(1))[0] #proppant permeability
        #scaling
        self.u_N = -1.0
        self.u_alpha = -1.0
        self.u_a = -1.0
        #self.u_b = -1.0
        self.u_gamma = -1.0
        #self.u_n1 = -1.0
        #hydraulic geometry
        self.bh = -1.0
        
        #*** stochastic sampled parameters ***
        self.u_gamma = lognorm_trunc(1,np.log10(rock.gamma[1]),np.log10(rock.gamma[1]),np.log10(rock.gamma[0]),np.log10(rock.gamma[2]))[0]
        self.u_a = norm_trunc(1,rock.a[1],0.5*rock.a[1],rock.a[0],rock.a[2])[0]
        self.u_N = contact_trunc(1,0.15,rock.N[1],0.5*rock.N[1],0.5,0.5*rock.N[1]/np.pi,rock.N[0],rock.N[2])[0]
        self.u_alpha = norm_trunc(1,rock.alpha[1],0.5*rock.alpha[1],rock.alpha[0],rock.alpha[2])[0]
        self.bh = norm_trunc(1,rock.bh[1],0.5*rock.bh[1],rock.bh[0],rock.bh[2])[0] #!!! would be nice to replace this with a physics based estimate
        if phi < 0:
            self.phi = np.random.uniform(rock.phi[0],rock.phi[2],(1))[0]
        else:
            self.phi = phi
        if mcc < 0:
            self.mcc = np.random.uniform(rock.mcc[0],rock.mcc[2],(1))[0]
        else:
            self.mcc = mcc
        #stress state
        self.Pc, self.sn, self.tau = rock.stress.Pc_frac(self.str, self.dip, self.phi, self.mcc)
        #apertures
        self.bd = self.bh/self.u_N
        self.bd0 = self.bd
        self.bd0p = 0.0
        self.vol = (4.0/3.0)*pi*0.25*self.dia**2.0*0.5*self.bd
        self.arup = 0.25*np.pi*self.dia**2.0
    #adjust fracture cohesion to prevent runaway stimulation at specified conditions
    def check_integrity(self,rock=reservoir(),pres=0.0):
        #originally assigned values
        Pc0 = self.Pc
        mcc0 = self.mcc
        phi0 = self.phi
        #shear stability limit
        mcc_c = self.tau - (self.sn - pres)*np.tan(self.phi)
        #tensile stability limit
        mcc_t = pres - self.sn
        #update fracture cohesion to ensure stability at the input conditions
        self.mcc = np.max([self.mcc,mcc_c,mcc_t])
        #recompute critical conditions
        self.Pc, self.sn, self.tau = rock.stress.Pc_frac(self.str, self.dip, self.phi, self.mcc)
        #output message
        if self.Pc > Pc0:
            print('alert: critical fracture at Tau = %.2e, sn = %.2e, and Pp %.2e having phi = %.2e rad, mcc = %.2e Pa, Pc = %.2e Pa adjusted to phi = %.2e rad, mcc = %.2e Pa, Pc = %.2e Pa'
                  %(self.tau,self.sn,pres,phi0,mcc0,Pc0,self.phi,self.mcc,self.Pc))
    #set fracture cohesion to critical value at specified conditions
    def make_critical(self,rock=reservoir(),pres=0.0):
        #shear stability limit
        mcc_c = self.tau - (self.sn - pres)*np.tan(self.phi)
        #tensile stability limit
        mcc_t = pres - self.sn
        #update fracture cohesion to enforce critical stability at the input conditions
        self.mcc = np.max([mcc_c,mcc_t])
        #recompute critical conditions
        self.Pc, self.sn, self.tau = rock.stress.Pc_frac(self.str, self.dip, self.phi, self.mcc)
        #output message
        print('         hydrofrac critical cohesion calculated as %.2e Pa to obtain Pc = %.2e Pa' %(self.mcc,self.Pc))
        #error case
        if self.Pc <= self.sn:
            print('         *** error: solution gives closed supercritical shear instead of hydrofracture')

#file-based dfn
class dfn:
    def __init__(self):
        self.data = [[]]
        self.names = []
        
    #create a file containing the fracture geometry from this network
    def export(self,filename='dfn_exp.csv',fracs=[surf()]):
        for i in range(0,len(fracs)):
            out = [['Number_int',i]]
            out += [['Easting_m',fracs[i].c0[0]]]
            out += [['Northing_m',fracs[i].c0[1]]]
            out += [['Height_m',fracs[i].c0[2]]]
            out += [['Diameter_m',fracs[i].dia]]
            out += [['Strike_rad',fracs[i].str]]
            out += [['Dip_rad',fracs[i].dip]]
            out += [['Type_type',fracs[i].typ]]
            out += [['FrictionAngle_rad',fracs[i].phi]]
            out += [['Cohesion_Pa',fracs[i].mcc]]
            out += [['ProppantLoad_m3',fracs[i].prop_load]]
            out += [['ProppantCompressibility_1pPa',fracs[i].prop_alpha]]
            out += [['Tortuosity_ratio',fracs[i].roughness]]
            out += [['ProppantPermeability_m2',fracs[i].kf]]
            out += [['HydraulicDilation_ratio',fracs[i].u_N]]
            out += [['Compressibility_1pPa',fracs[i].u_alpha]]
            out += [['DilationSlip_ratio',fracs[i].u_a]]
            out += [['SlipLength_ratio',fracs[i].u_gamma]]
            out += [['InitialHydraulicAperture_m',fracs[i].bh]]
            if i < 1:
                save_csv(out=out,filename=filename,append=False)
            else:
                save_csv(out=out,filename=filename,append=True)
            
    #load fracture information and add to network
    def intake(self,filename='dfn.csv',fracs=[surf()],rock=reservoir()):
        self.names, self.data = load_csv(filename)
        if True:
            for i in range(0,len(self.data)):
                fracs += [surf(self.data['Easting_m'][i],
                                    self.data['Northing_m'][i],
                                    self.data['Height_m'][i],
                                    self.data['Diameter_m'][i],
                                    self.data['Strike_rad'][i],
                                    self.data['Dip_rad'][i],
                                    typ(int(self.data['Type_type'][i]+0.001)),
                                    rock)]
                if 'FrictionAngle_rad'            in self.names: fracs[-1].phi = self.data['FrictionAngle_rad'][i]
                if 'Cohesion_Pa'                  in self.names: fracs[-1].mcc = self.data['Cohesion_Pa'][i]
                if 'ProppantLoad_m3'              in self.names: fracs[-1].prop_load = self.data['ProppantLoad_m3'][i]
                if 'ProppantCompressibility_1pPa' in self.names: fracs[-1].prop_alpha = self.data['ProppantCompressibility_1pPa'][i]
                if 'Tortuosity_ratio'             in self.names: fracs[-1].roughness = self.data['Tortuosity_ratio'][i]
                if 'ProppantPermeability_m2'      in self.names: fracs[-1].kf = self.data['ProppantPermeability_m2'][i]
                if 'HydraulicDilation_ratio'      in self.names: fracs[-1].u_N = self.data['HydraulicDilation_ratio'][i]
                if 'Compressibility_1pPa'         in self.names: fracs[-1].u_alpha = self.data['Compressibility_1pPa'][i]
                if 'DilationSlip_ratio'           in self.names: fracs[-1].u_a = self.data['DilationSlip_ratio'][i]
                if 'SlipLength_ratio'             in self.names: fracs[-1].u_gamma = self.data['SlipLength_ratio'][i]
                if 'InitialHydraulicAperture_m'   in self.names: fracs[-1].bh = self.data['InitialHydraulicAperture_m'][i]
        else:
            print('%s improperly formatted' %(filename))        

#line objects
class line:
    def __init__(self,x0=0.0,y0=0.0,z0=0.0,length=1.0,azn=0.0*deg,dip=0.0*deg,w_type='pipe',
                 ra=0.0254*3.0,rb=0.0254*3.5,rc=0.0254*3.5,rough=80.0,pID=0):
        #position geometry
        self.c0 = np.asarray([x0, y0, z0]) #origin
        self.leg = length #length
        self.azn = azn #axis azimuth north
        self.dip = dip #axis dip from horizontal
        self.typ = typ(w_type) #type of well
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

#node list object
class nodes:
    #initialization
    def __init__(self):
        self.r0 = np.asarray([np.inf,np.inf,np.inf])
        self.all = np.asarray([self.r0])
        self.tol = 0.0005
        self.num = len(self.all)
        self.p = np.zeros(self.num,dtype=float)
        self.T = np.zeros(self.num,dtype=float)
        self.h = np.zeros(self.num,dtype=float)
        self.Tr = np.zeros(self.num,dtype=float)
#        self.f_id = [[-1]]
    #add a node
    def add(self,c=np.asarray([0.0,0.0,0.0])): #,f_id=-1):
        #round to within tolerance
        if not(np.isinf(c[0])):
            c = np.rint(c/self.tol)*self.tol
        #check for duplicate existing node
        ck_1 = np.asarray([np.isin(self.all[:,0],c[0]),np.isin(self.all[:,1],c[1]),np.isin(self.all[:,2],c[2])])
        ck_2 = ck_1[0,:]*ck_1[1,:]*ck_1[2,:]
        ck_i = np.where(ck_2 == 1)
        #yes duplicate -> return index of existing node
        if len(ck_i[0]) > 0:
#            if f_id != -1:
#                self.f_id[ck_i[0][0]] += [f_id]
            return False, ck_i[0][0]
        #no duplicate -> add node -> return index of new node
        else: 
            self.all = np.concatenate((self.all,np.asarray([c])),axis=0)
            self.num = len(self.all)
            self.p = np.concatenate((self.p,np.asarray([0.0])),axis=0)
            self.T = np.concatenate((self.T,np.asarray([0.0])),axis=0)
            self.h = np.concatenate((self.h,np.asarray([0.0])),axis=0)
            self.Tr = np.concatenate((self.Tr,np.asarray([0.0])),axis=0)
#            self.f_id += [[f_id]]
            return True, len(self.all)-1

#pipe list object
class pipes:
    #initialization
    def __init__(self):
        self.num = 0
        self.n0 = [] #source node
        self.n1 = [] #target node
        self.L = [] #length
        self.W = [] #width
        self.typ = [] #property source type
        self.fID = [] #property source index
        self.pID = [] #parent index
        self.K = [] #flow solver coefficient
        self.n = [] #flow solver exponent
        self.Dh = [] #hydraulic aperture/diameter
        self.Dh_max = [] #hydraulic aperture/diameter limit
        self.frict = [] #hydraulic roughness
        self.hydrofraced = False #tracker for fracture initiation
        self.R0 = [] #thermal radius
        self.q = [] #flowrate
    #add a pipe
    def add(self, n0, n1, length, width, featTyp, featID, Dh=1.0, Dh_max=1.0, frict=1.0, pID = 0):
        self.n0 += [n0]
        self.n1 += [n1]
        self.L += [length]
        self.W += [width]
        self.typ += [featTyp]
        self.fID += [featID]
        self.pID += [pID]
        self.K += [1.0]
        self.n += [1.0]
        self.num = len(self.n0)
        self.Dh += [Dh]
        self.Dh_max += [Dh_max]
        self.frict += [frict]
        self.R0 += [0.0]
        self.hydrofraced = False
        self.q += [0.0]
    #set hydraulic aperture limits based on minimum pressure drop at target flow rate
    def Dh_limit(self,Q,dP,rho=980.0,g=9.81,mu=0.9*cP): #,k=0.1*mD):
        for i in range(0,self.num):
            if (int(self.typ[i]) in [typ('injector'),typ('producer'),typ('pipe'),typ('perfcluster'),typ('screen')]):
                Lscaled = self.L[i]*mu/(0.9*cP)
                a_max = (10.7e-4*Lscaled*rho*g*Q/(dP*self.frict[i]**1.852))**(1.0/4.87)
                self.Dh_max[i] = a_max
            elif (int(self.typ[i]) in [typ('boundary'),typ('fracture'),typ('propped'),typ('choke')]):
                bh_max = (12.0e-4*mu*Q*self.L[i]/(dP*self.W[i]))**(1.0/3.0)
                self.Dh_max[i] = bh_max
            elif (int(self.typ[i]) in [typ('perf')]):
                dp_max = (4.0e-4*mu*rho*Q**2.0/(0.9*cP*dP*self.frict[i]**2.0*0.825**2.0))**(0.25)
                self.Dh_max[i] = dp_max
            elif (int(self.typ[i]) in [typ('darcy')]):
                # t_max = 1.0e-4*Q*mu*self.L[i]/(k*self.W[i]*dP)
                t_max = 1.0e-4*Q*mu*self.L[i]/(self.frict[i]*self.W[i]*dP)
                self.Dh_max[i] = t_max
            elif (int(self.typ[i]) in [typ('shunt')]):
                s_max = 1e-10
                self.Dh_max[i] = s_max
            else:
                print( 'error: undefined type of conduit')
                exit
    #get coordinates
    def coords(self,nodes):
        self.c0 = []
        self.c1 = []
        for p in range(0,self.num):
            self.c0 += [nodes.all[self.n0[p]]]
            self.c1 += [nodes.all[self.n1[p]]]
        
#model object, functions, and data as object
class core:
    def __init__(self): #,node=[],pipe=[],fracs=[],wells=[],hydfs=[],bound=[],geo3D=[]): #@@@ are these still used?
        #domain information
        self.rock = reservoir()
        self.nodes = nodes()
        self.pipes = pipes()
        self.fracs = []
        self.wells = []
        self.hydfs = []
        self.bound = []
        self.faces = []
        #intersections tracker
        self.trakr = [] #index of fractures in chain
        #stimulation solver
        self.iters = 0
        self.pin = np.random.randint(100000000,999999999,1)[0]
        self.Qi = []
        self.Pi = []
        self.Vi = []
        self.ti = []
        self.Mi = []
        self.progress = []
        #flow solver
        self.H = [] #boundary pressure head array, m
        self.Q = [] #boundary flow rate array, m3/s
        self.q = [] #calculated pipe flow rates
        self.v5 = [] #calculated inlet specific volume, m3/kg
        # self.i_p = [] #constant flow well pressures, Pa
        # self.i_q = [] #constant flow well rates, m3/s
        self.p_p = [] #constant pressure well pressures, Pa
        self.p_q = [] #constant pressure well rates, m3/s
        self.b_p = [] #boundary pressure, Pa
        self.b_q = [] #boundary rates, m3/s
        #heat solver
        self.Tb = [] #boundary temperatures
        self.R0 = [] #thermal radius
        self.Rt = [] #thermal radius over time
        self.ms = [] #pipe mass flow rates
        self.Et = [] #energy in rock over time
        self.Qt = [] #heat flow from rock over time
        self.Tt = [] #heat flow from rock over time
        self.ht = [] #enthalpy over time
        self.ts = [] #time stamps
        self.w_h = [] #wellhead enthalpy
        self.w_m = [] #wellhead mass flow
        self.p_E = [] #production energy
        self.b_h = [] #boundary enthalpy
        self.b_m = [] #boundary mass flow
        self.b_E = [] #boundary energy
        self.i_E = [] #injection energy
        self.i_mm = [] #mixed injection mass flow rate
        self.p_mm = [] #mixed produced mass flow rate
        self.p_hm = [] #mixed produced enthalpy
        #power solver
        self.Qout = [] #pumping power out
        self.Pout = [] #net power out
        self.Tout = [] #bulk power out
        self.Hout = [] #thermal power out
        # self.dhout = [] #heat extraction
        #validation stuff
        self.v_Rs = []
        self.v_ts = []
        self.v_Ps = []
        self.v_ws = []
        self.v_Vs = []
        self.v_Pn = []
        #economics
        self.eNPV = 0.0 #net present value
        self.eC = 0.0 #capital cost
        self.eQ = 0.0 #quake cost
        self.eG = 0.0 #generation cost
        self.eP = 0.0 #net power sales
        self.eB = 0.0 #bulk power sales
        self.CpM_elec = 0.0 #cost per net MWe-h
        self.CpM_prod = 0.0 #cost per produced MWe-h
        self.CpM_ther = 0.0 #cost per MWt-h
        self.drill_length = 0.0 #total drilled length
        #key outputs
        self.payoff = payoff()
        
    # Propped fracture property estimation with geomechanics
    def hydromech(self, f_id, fix=False): #, pp=-666.0):
        #initialize stim
        stim = False
        
        #check if fracture will be stimulated
        if (self.faces[f_id].Pmax >= self.faces[f_id].Pc) and (not(fix)):
            self.faces[f_id].stim += 1
            stim = True
            print( '   + fracture stimulated: %i' %(f_id))

        #pressures for analysis
        e_max = self.faces[f_id].sn - self.faces[f_id].Pmax
        
        #original bd0
        bd0 = np.max([1.0e-9,self.faces[f_id].bd0])

        #stimulation enabled
        if stim and not(fix):
            #*** shear ***
            #don't let shear be zero
            if self.faces[f_id].tau < 1.0:
                self.faces[f_id].tau = 1.0
            #maximum moment (shear stress method)
            M0max = self.faces[f_id].tau*self.faces[f_id].arup**(3.0/2.0)
            Mwmax = (np.log10(M0max)-9.1)/1.5
            #mw sample from G-R
            mw = exponential_trunc(1,bval=self.rock.bval,Mmax=Mwmax,Mwin=1.0,prob=0.1)[0]
            #convert to moment magnitude (Mo)
            mo = 10.0**(mw * 1.5 + 9.1)
            # rupture length
            Lr = (4*((mo/self.faces[f_id].tau)**(2/3))/np.pi)**0.5
            # intermediate variables
            ds = mo / ((0.25*np.pi*Lr**2.0) * self.rock.ResG)
            # d0 = 0.5 * self.faces[f_id].u_gamma * (self.faces[f_id].dia ** self.faces[f_id].u_n1)
            d0 = 0.5 * self.faces[f_id].u_gamma * self.faces[f_id].dia
            # bd0 = self.faces[f_id].u_a * (d0 ** self.faces[f_id].u_b)
            bd0 = self.faces[f_id].u_a * d0
            d1 = d0 + ds
            # bd1 = self.faces[f_id].u_a * (d1 ** self.faces[f_id].u_b)
            bd1 = self.faces[f_id].u_a * d1
            dbd = bd1- bd0
            #record stimulation magnitude and correct for seismic overestimation
            self.faces[f_id].Mws += [mw]
            #add to zero-stress dilatant aperture
            bd0 = self.faces[f_id].bd0 + dbd
            #override rupture length in tensile fractures to avoid overpredicting tensile seismicity
            if (e_max < 0.0):
                Lr = self.faces[f_id].dia
                
            #*** growth ***
            #grow fracture by larger of 20% fracture size or 5% domain size
            add_dia = np.max([0.2*self.faces[f_id].dia,0.05*self.rock.size])
            #deduct event's rupture area from residual rupture area
            self.faces[f_id].arup = np.max([self.faces[f_id].arup - 0.25*np.pi*Lr**2.0,  0.1])           
            #update rupture area
            self.faces[f_id].arup = self.faces[f_id].arup + 0.25*np.pi*((self.faces[f_id].dia + add_dia)**2.0 - self.faces[f_id].dia**2.0)
            #update fracture size
            self.faces[f_id].dia += add_dia
            
        #fracture parameters
        f_radius = 0.5*self.faces[f_id].dia
        #propped closure with ideal loose spherical packing (Allen, 1985; Frings et al., 2011)
        prop_load = np.max([0.0,self.faces[f_id].prop_load])
        bd0p = (prop_load/0.64)/(0.25*np.pi*self.faces[f_id].dia**2.0)
        bd0p = np.max([1e-9,bd0p])
        #closure pressure onto a proppant pack
        e_crit = (bd0p*np.pi*self.rock.ResE)/(-8.0*(1.0-self.rock.Resv**2.0)*f_radius)
        e_crit = np.min([-1e-3,e_crit])
        #hydropropped fracture with proppant pillars
        if (e_max < e_crit):
            #maximum aperture from Sneddon's (PKN) penny fracture aperture
            bdt = (-8.0*e_max*(1.0-self.rock.Resv**2.0)*f_radius)/(pi*self.rock.ResE) #m
            #channel width ratios by filled volume to total volume
            wp_wt = bd0p/bdt
            wo_wt = 1.0 - wp_wt
            #flow through open channels
            qo = wo_wt*(self.faces[f_id].roughness * bdt)**3.0
            #flow through propped area
            qp = wp_wt*(12.0 * self.faces[f_id].kf * bdt)
            #flow through shear channels
            qs = (self.faces[f_id].u_N * bd0)**3.0
            #total dilated aperture
            bd = bdt + bd0
            #total hydraulic aperture
            bh = (qo + qp + qs)**(1.0/3.0)
            #note that fracture is hydropropped
            self.faces[f_id].hydroprop = True
        #proppant propped fracture 
        elif (e_max < 0.0):
            #flow through propped area
            qp = 12.0 * self.faces[f_id].kf * bd0p*np.exp(-self.faces[f_id].prop_alpha*(e_max-e_crit))
            #flow through shear channels
            qs = (self.faces[f_id].u_N * bd0)**3.0
            #total dilated aperture
            bd = bd0p*np.exp(-self.faces[f_id].prop_alpha*(e_max-e_crit)) + bd0
            #total hydraulic aperture
            bh = (qp + qs)**(1.0/3.0)
            #note that fracture is not hydropropped
            self.faces[f_id].hydroprop = False
        #proppant propped fracture 
        else:
            #flow through propped area
            qp = 12.0 * self.faces[f_id].kf * bd0p*np.exp(-self.faces[f_id].prop_alpha*(e_max-e_crit))
            #flow through shear channels
            qs = (self.faces[f_id].u_N * bd0 * np.exp(-self.faces[f_id].u_alpha*e_max))**3.0
            #total dilated aperture
            bd = (bd0p*np.exp(-self.faces[f_id].prop_alpha*(e_max-e_crit)) +
                  bd0*np.exp(-self.faces[f_id].u_alpha*e_max))
            #total hydraulic aperture
            bh = (qp + qs)**(1.0/3.0)
            #note that fracture is not hydropropped
            self.faces[f_id].hydroprop = False
            
        #override for boundary fractures
        if (int(self.faces[f_id].typ) in [typ('boundary')]):
            bh = self.rock.bh_bound
        
        #overrides for bad math
        if not(np.isfinite(bh)):
            bh = 1e-9
        elif bh < 1e-9:
            bh = 1e-9
        
        #volume
        vol = (4.0/3.0)*pi*0.25*self.faces[f_id].dia**2.0*0.5*bd
        
        #update proppant loading
        dvol = vol - self.faces[f_id].vol
        self.faces[f_id].prop_load = self.faces[f_id].prop_load + np.max([0.0,dvol])*self.rock.sand

        #update fracture properties
        self.faces[f_id].bd = bd
        self.faces[f_id].bh = bh
        self.faces[f_id].vol = vol
        self.faces[f_id].bd0 = bd0
        self.faces[f_id].bd0p = bd0p 
        
        #return true if stimulated
        return stim
    
    def pre_save(self):
        self.payoff.Drilling_BottomTemperature_K = self.rock.BH_T
        self.payoff.Drilling_BottomPressure_MPa = self.rock.BH_P
        self.payoff.Drilling_WellDiameter_m = 2.0*self.rock.rc
        self.payoff.Drilling_CasingWeight_kgpm = 7850*np.pi*(self.rock.rb**2 - self.rock.ra**2)
        self.payoff.Drilling_CementVolume_m3 = (self.payoff.Drilling_TotalLength_m * np.pi *
                                                (self.rock.rc**2 - self.rock.rb**2))
        self.payoff.Stimulation_SandMass_kg = self.rock.Vstim*self.rock.sand*1538.0
        self.payoff.Stimulation_WaterVolume_m3 = self.rock.Vstim
        self.payoff.Stimulation_PumpPressure_MPa = 0.0 #!!!
        self.payoff.Circulation_WorkingVolume_m3 = 0.0 #!!!
        if self.payoff.Injection_Rate_kgps != 0:
            self.payoff.Circulation_Recovery_pct = 100*self.payoff.Production_Rate_kgps/self.payoff.Injection_Rate_kgps
        else:
            self.payoff.Circulation_Recovery_pct = 100.0
        self.payoff.Circulation_MakeupWater_m3 = ((self.payoff.Injection_Rate_kgps-self.payoff.Production_Rate_kgps)*
                                                  self.rock.LifeSpan)
        self.payoff.Circulation_PumpPressure_MPa = self.rock.p_whp/MPa
        self.payoff.Stress_s3_MPa = self.rock.s3/MPa
        self.payoff.Stress_s2_MPa = self.rock.s2/MPa
        self.payoff.Stress_s1_MPa = self.rock.s1/MPa
        self.payoff.Stress_pp_MPa = self.rock.BH_P/MPa
        self.payoff.Injection_PressureLimit_MPa = self.rock.pfinal_max/MPa
        num_s = 0
        num_t = 0
        num_h = 0
        for face in self.faces:
            #face.stim = 0
            if face.hydroprop:
                    num_h += 1
            elif typ(face.typ) in ['propped']:
                num_t += 1
            elif face.stim > 0:
                num_s += 1
        self.payoff.Fractures_ShearStim_ea = num_s
        self.payoff.Fractures_HydroStim_ea = num_t
        self.payoff.Fractures_Hydroprop_ea = num_h
    
    def setup_payoff(self,fname='setup_payoff.csv',pin='',aux=[],append=True):
        self.pre_save()
        out = []
        out += [['pin',pin]]
        #inputs
        for var in vars(self.rock.setup):
            out += [[var,getattr(self.rock.setup,var)]]
        #key outputs
        for var in vars(self.payoff):
            out += [[var,getattr(self.payoff,var)]]
        #timeseries
        if len(self.ts)>0:
            pkg_ts = '%.3e' %(self.ts[0])
            try:
                pkg_Ts = '%.3e' %(water.T_from_Ph(P=self.rock.p_whp,h=self.p_hm[0]))
            except:
                pkg_Ts = '0.000'
            pkg_hs = '%.3e' %(self.p_hm[0])
            pkg_Ps = '%.3e' %(self.Pout[0])
            for i in range(1,self.rock.TimeSteps-1):
                pkg_ts += ';%.3e' %(self.ts[i])
                try:
                    pkg_Ts += ';%.3e' %(water.T_from_Ph(h=self.p_hm[i],P=self.rock.p_whp))
                except:
                    pkg_Ts += ';0.000'
                pkg_hs += ';%.3e' %(self.p_hm[i])
                pkg_Ps += ';%.3e' %(self.Pout[i])
        else:
            pkg_ts = '0.000'
            pkg_Ts = '0.000'
            pkg_hs = '0.000'
            pkg_Ps = '0.000'
        out += [['Series_Time_s',pkg_ts]]
        out += [['Series_Temperature_s',pkg_Ts]]
        out += [['Series_Enthalpy_s',pkg_hs]]
        out += [['Series_Power_s',pkg_Ps]]
        #aux
        if aux:
            out += aux
        #to file
        save_csv(out=out,filename=fname,append=append)
    
#     def save(self,fname='input_output.txt',pin='',aux=[],printwells=0,time=True):
#         out = []
        
#         out += [['pin',pin]]
        
#         #input parameters
#         r = self.rock
#         out += [['size',r.size]]
#         out += [['gradient',r.gradient]]
#         out += [['ResDepth',r.ResDepth]]
#         out += [['ResGradient',r.ResGradient]]
#         out += [['ResRho',r.ResRho]]
#         out += [['ResKt',r.ResKt]]
#         out += [['ResSv',r.ResSv]]
#         out += [['AmbTempK',r.AmbTempK]]
#         out += [['AmbPres',r.AmbPres]]
#         out += [['ResE',r.ResE]]
#         out += [['Resv',r.Resv]]
#         out += [['ResG',r.ResG]]
#         out += [['Ks3',r.Ks3]]
#         out += [['Ks2',r.Ks2]]
#         out += [['s3Azn',r.s3Azn]]
#         out += [['s3AznVar',r.s3AznVar]]
#         out += [['s3Dip',r.s3Dip]]
#         out += [['s3DipVar',r.s3DipVar]]
#         for i in range(0,len(r.fNum)):
#             out += [['fNum%i' %(i),r.fNum[i]]]
#             out += [['fDia_min%i' %(i),r.fDia[i][0]]]
#             out += [['fDia_max%i' %(i),r.fDia[i][1]]]
#             out += [['fStr_nom%i' %(i),r.fStr[i][0]]]
#             out += [['fStr_var%i' %(i),r.fStr[i][1]]]
#             out += [['fDip_nom%i' %(i),r.fDip[i][0]]]
#             out += [['fDip_var%i' %(i),r.fDip[i][1]]]
#         for i in range(0,3):
#             out += [['alpha%i' %(i),r.alpha[i]]]
#         for i in range(0,3):
#             out += [['prop_alpha%i' %(i),r.prop_alpha[i]]]
#         for i in range(0,3):
#             out += [['gamma%i' %(i),r.gamma[i]]]
#         for i in range(0,3):
#             out += [['a%i' %(i),r.a[i]]]
#         for i in range(0,3):
#             out += [['N%i' %(i),r.N[i]]]
#         for i in range(0,3):
#             out += [['bh%i' %(i),r.bh[i]]]
#         out += [['bh_bound',r.bh_bound]]
#         for i in range(0,3):
#             out += [['f_roughness%i' %(i),r.f_roughness[i]]]
#         out += [['w_count',r.w_count]]
#         out += [['w_spacing',r.w_spacing]]
#         out += [['w_length',r.w_length]]
#         out += [['w_azimuth',r.w_azimuth]]
#         out += [['w_dip',r.w_dip]]
#         out += [['w_proportion',r.w_proportion]]
#         out += [['w_phase',r.w_phase]]
#         out += [['w_toe',r.w_toe]]
#         out += [['w_skew',r.w_skew]]
#         out += [['w_intervals',r.w_intervals]]
#         out += [['ra',r.ra]]
#         out += [['rb',r.rb]]
#         out += [['rc',r.rc]]
#         out += [['rgh',r.rgh]]
#         out += [['CemKt',r.CemKt]]
#         out += [['CemSv',r.CemSv]]
#         out += [['GenEfficiency',r.GenEfficiency]]
#         out += [['LifeSpan',r.LifeSpan]]
#         out += [['TimeSteps',r.TimeSteps]]
#         out += [['p_whp',r.p_whp]]
#         out += [['TinjK',r.TinjK]]
#         out += [['H_ConvCoef',r.H_ConvCoef]]
#         out += [['dT0',r.dT0]]
#         out += [['PoreRho',r.PoreRho]]
#         out += [['Poremu',r.Poremu]]
#         out += [['Porek',r.Porek]]
#         for i in range(0,3):
#             out += [['kf%i' %(i),r.kf[i]]]
#         #out += [['kf',r.kf]]
#         out += [['BH_T',r.BH_T]]
#         out += [['BH_P',r.BH_P]]
#         out += [['s1',r.s1]]
#         out += [['s2',r.s2]]
#         out += [['s3',r.s3]]
#         out += [['perf_clusters',r.perf_clusters]]
#         out += [['perf_dia',r.perf_dia]]
#         out += [['perf_per_cluster',r.perf_per_cluster]]
#         out += [['r_perf',r.r_perf]]
#         out += [['sand',r.sand]]
#         out += [['dPi',r.dPi]]
#         out += [['stim_limit',r.stim_limit]]
#         out += [['Qinj',r.Qinj]]
#         out += [['Vinj',r.Vinj]]
#         out += [['Qstim',r.Qstim]]
#         out += [['Vstim',r.Vstim]]
#         out += [['pfinal_max',r.pfinal_max]]
#         out += [['bval',r.bval]]
#         for i in range(0,3):
#             out += [['phi%i' %(i),r.phi[i]]]
#         for i in range(0,3):
#             out += [['mcc%i' %(i),r.mcc[i]]]
#         out += [['hfmcc',r.hfmcc]]
#         out += [['hfphi',r.hfphi]]
#         out += [['sales_kWh',r.sales_kWh]]
#         out += [['drill_m',r.drill_m]]
#         out += [['pad_fixed',r.pad_fixed]]
#         out += [['plant_kWe',r.plant_kWe]]
#         out += [['explore_m',r.explore_m]]
#         out += [['oper_kWh',r.oper_kWh]]
#         out += [['quake_coef',r.quake_coef]]
#         out += [['quake_exp',r.quake_exp]]
#         out += [['NPV',self.eNPV]]
#         out += [['CapitalCost',self.eC]]
#         out += [['QuakeRisk',self.eQ]]
#         out += [['GenerationCost',self.eG]]
#         out += [['NetSales',self.eP]]
#         out += [['MaxSales',self.eB]]
#         out += [['CpM_net',self.CpM_elec]]
#         out += [['CpM_gro',self.CpM_prod]]
#         out += [['CpM_the',self.CpM_ther]]
#         out += [['drill_length',self.drill_length]]
        
#         #auxillary printed inputs & outputs
#         if aux:
#             out += aux
            
#         #compile well types
#         wtyp = np.zeros(len(self.wells))
#         for i in range(0,len(wtyp)):
#             wtyp[i] = self.wells[i].typ
        
#         #total injection rate (+ in)
#         qinj = 0.0
#         # key = np.where(np.asarray(self.i_q)>0.0)[0]
#         # for i in key:
#         #     if wtyp[i] < 0:
#         #         qinj += self.i_q[i]
#         key = np.where(np.asarray(self.p_q)>0.0)[0]
#         for i in key:
#             # if wtyp[i] < 0:
#             if (int(wtyp[i]) in [typ('injector')]):
#                 qinj += self.p_q[i]
#         out += [['qinj',qinj]]
        
#         #total production rate (-out)
#         qpro = 0.0
#         # key = np.where(np.asarray(self.i_q)<0.0)[0]
#         # for i in key:
#         #     if wtyp[i] < 0:
#         #         qpro += self.i_q[i]
#         key = np.where(np.asarray(self.p_q)<0.0)[0]
#         for i in key:
#             # if wtyp[i] < 0:
#             if (int(wtyp[i]) in [typ('producer')]):
#                 qpro += self.p_q[i]
#         out += [['qpro',qpro]]
        
#         #total leakoff rate (-out)
#         qoff = 0.0
#         key = np.where(np.asarray(self.b_q)<0.0)[0]
#         for i in key:
#             qoff += self.b_q[i]
#         out += [['qleak',qoff]]
        
#         #total boundary uptake (+in)
#         qup = 0.0
#         key = np.where(np.asarray(self.b_q)>0.0)[0]
#         for i in key:
#             qup += self.b_q[i]
#         out += [['qgain',qup]]
        
#         #recovery
#         qrec = 0.0
#         if qinj > 0:
#             qrec = -qpro/qinj
#         else:
#             qrec = 1.0
#         out += [['recovery',qrec]]
        
#         #largest quake
#         quake = -10.0
#         for i in range(0,len(self.faces)):
#             if self.faces[i].Mws:
#                 quake = np.max([quake,np.max(self.faces[i].Mws)])
#         out += [['max_quake',quake]]
        
#         #injection pressure
#         pinj = np.max(self.nodes.p) - r.BH_P
#         out += [['pinj',pinj]]
        
#         #injection enthalpy
#         T5 = self.rock.TinjK # K
#         # if pinj > 100.0:
#         #     pinj = 100.0
#         try:
#             hinj = water.h_from_PT(P=pinj,T=T5)
#         except:
#             hinj = 0.0
#         out += [['hinj',hinj]]
        
#         #injection specific volume
#         out += [['v5',self.v5]]
        
#         #injection intercepts
#         ixint = 0
#         for i in range(0,self.pipes.num):
#             if (int(self.pipes.typ[i]) in [typ('injector'),typ('perfcluster')]):
#                 ixint += 1
#         ixint = ixint - self.rock.w_intervals
#         out += [['ixint',ixint]]
        
#         #production intercepts
#         pxint = 0
#         for i in range(0,self.pipes.num):
#             if (int(self.pipes.typ[i]) in [typ('producer'),typ('screen')]):
#                 pxint += 1
#         pxint = pxint - self.rock.w_count
#         out += [['pxint',pxint]]
        
#         #stimulated fractures
#         hfstim = 0
#         nfstim = 0
#         for i in range(0,len(self.faces)):
#             if (self.faces[i].stim > 0):
#                 if (int(self.faces[i].typ) in [typ('propped')]):
#                     hfstim += 1
#                 elif (int(self.faces[i].typ) in [typ('fracture')]):
#                     nfstim += 1
#         out += [['hfstim',hfstim]]
#         out += [['nfstim',nfstim]]

#         #injection mass flow rate
#         out += [['minj',self.i_mm]]
        
#         #production mass flow rate
#         out += [['mpro',self.p_mm]]
        
#         #production enthalpy over time
#         for t in range(0,len(self.ts)-1):
#             out += [['hpro:%.3f' %(self.ts[t]/yr),self.p_hm[t]]]
        
#         #flash power
#         Pout = []
#         if len(self.Pout) < 1:
#             Pout = np.full(len(self.ts),np.nan)
#         else:
#             Pout = self.Pout
#         for t in range(0,len(self.ts)-1):
#             out += [['Pout:%.3f' %(self.ts[t]/yr),Pout[t]]]
        
#         #per well values
#         if (printwells != 0) and (self.w_m.any()):
#             # #per well volume flow rate
#             # dummy = np.zeros(20,dtype=float)
#             # for i in range(0,len(self.p_q)):
#             #     if i > 20:
#             #         print('warning: number of wells exceeds placeholder of 20 so not all values will be printed')
#             #         break
#             #     dummy[i] = self.p_q[i]
#             # for i in range(0,len(dummy)):
#             #     out += [['q%i' %(i),dummy[i]]]
    
#             #collect well temp (T), volume rate (q), mass rate (m), and enthalpy (h)
#             w_h = []
#             w_T = []
#             w_m = []
#             w_q = []
#             for w in range(0,len(self.wells)):
#                 #coordinates
#                 source = self.wells[w].c0
#                 #find index of duplicate
#                 ck, i = self.nodes.add(source)
#                 #record temperature
#                 w_T += [self.Tt[:,i]]
#                 #record enthalpy
#                 w_h += [self.ht[:,i]]
#                 #record mass flow rate
#                 i_pipe = np.where(np.asarray(self.pipes.n0) == i)[0][0]
#                 w_m += [self.q[i_pipe]/self.v5]
#                 #record volume flow rate
#                 w_q += [self.q[i_pipe]]
#             w_h = np.asarray(w_h)
#             w_T = np.asarray(w_T)
#             w_m = np.asarray(w_m)
#             w_q = np.asarray(w_q)
            
#             #per well mass flow rate #!!!
#             max_w = printwells
#             w_q_s = np.zeros(max_w,dtype=float)
#             w_m_s = np.zeros(max_w,dtype=float)
#             w_T_s = np.zeros((max_w,r.TimeSteps-1),dtype=float)
#             w_h_s = np.zeros((max_w,r.TimeSteps-1),dtype=float)
#             if len(self.wells) > max_w:
#                 print('warning: number of wells exceeds specified output number of %i so not all values will be printed' %(max_w))
#             #use array math to store values in fixed-width arrays
#             if self.wells: #error handling for case with no wells
#                 n = np.min([len(self.wells),max_w])
#                 w_q_s[:n] = w_q_s[:n] + w_q[:n]
#                 w_m_s[:n] = w_m_s[:n] + w_m[:n]
#                 w_T_s[:n] = w_T_s[:n,:] + w_T[:n,:-2]
#                 w_h_s[:n] = w_h_s[:n,:] + w_h[:n,:-2]
#             #rearrange into single-line format
#             for w in range(0,max_w):
#                 #volume flow rate, m3/s
#                 out += [['q%i' %(w),w_q_s[w]]]
#                 #mass flow rate, kg/s
#                 out += [['m%i' %(w),w_m_s[w]]]
#                 if time:
#                     for t in range(0,r.TimeSteps-1):
#                         #temperature over time, K
#                         out += [['T%i:%.3f' %(w,self.ts[t]/yr), w_T_s[w,t]]]
#                     for t in range(0,r.TimeSteps-1):    
#                         #enthalpy over time, h
#                         out += [['h%i:%.3f' %(w,self.ts[t]/yr), w_h_s[w,t]]]
        
#         #output to file
#         # out = zip(*out)
#         out = list(map(list,zip(*out)))
#         head = out[0][0]
#         for i in range(1,len(out[0])):
#             head = head + ',' + out[0][i]
#         data = '%i' %(out[1][0])
#         for i in range(1,len(out[1])):          
#             if not out[1][i]:
#                 data = data + ',0.0' # ',nan'
#             else:
#                 data = data + ',%.5e' %(out[1][i])
#         try:
#             with open(fname,'r') as f:
#                 test = f.readline()
#             f.close()
#             if test != '':
#                 with open(fname,'a') as f:
#                     f.write(data + '\n')
#                 f.close()
#             else:
#                 with open(fname,'a') as f:
#                     f.write(head + '\n')
#                     f.write(data + '\n')
#                 f.close()
#         except:
#             with open(fname,'a') as f:
#                 f.write(head + '\n')
#                 f.write(data + '\n')
#             f.close()
# #        
# #        return out
        
    def build_vtk(self,fname='default',vtype=[1,1,1,1,1,1]):
        #******   scaling       ******
        r = 0.002*self.rock.size
        
        #******   paint wells   ******
        if vtype[0]:
            w_obj = [] #fractures
            w_col = [] #fractures colors
            w_lab = [] #fractures color labels
            w_lab = ['Well_Number','Well_Type','Inner_Radius','Roughness','Outer_Radius']
            w_0 = []
            w_1 = []
            w_2 = []
            w_3 = []
            w_4 = []
            #nodex = np.asarray(self.nodes)
            for i in range(0,len(self.wells)): #skip boundary node at np.inf
                #add colors
                w_0 += [i]
                w_1 += [self.wells[i].typ]
                w_2 += [self.wells[i].ra]
                w_3 += [self.wells[i].rgh]
                w_4 += [self.wells[i].rc]
                #add geometry
                azn = self.wells[i].azn
                dip = self.wells[i].dip
                leg = self.wells[i].leg
                vAxi = np.asarray([math.sin(azn)*math.cos(-dip), math.cos(azn)*math.cos(-dip), math.sin(-dip)])
                c0 = self.wells[i].c0
                c1 = c0 + vAxi*leg
                w_obj += [sg.cylObj(x0=c0, x1=c1, r=1.5*r)]
            #vtk file
            w_col = [w_0,w_1,w_2,w_3,w_4]
            sg.writeVtk(w_obj, w_col, w_lab, vtkFile=(fname + '_wells.vtk'))
        
        #******   paint fractures   ******
        if vtype[1] and len(self.faces)>6:
            f_obj = [] #fractures
            f_col = [] #fractures colors
            f_lab = [] #fractures color labels
            f_lab = ['Face_Number','Node_Number','Type','Sn_MPa','Pc_MPa','Tau_MPa']
            f_0 = []
            f_1 = []
            f_2 = []
            f_3 = []
            f_4 = []
            f_5 = []
            #nodex = np.asarray(self.nodes)
            if len(self.faces) > 6:
                for i in range(6,len(self.faces)): #skip boundary node at np.inf
                    #add colors
                    f_0 += [i]
                    f_1 += [self.faces[i].ci]
                    f_2 += [self.faces[i].typ]
                    f_3 += [self.faces[i].sn/MPa]
                    f_4 += [self.faces[i].Pc/MPa]
                    f_5 += [self.faces[i].tau/MPa]
                    #add geometry
                    f_obj += [HF(r=0.5*self.faces[i].dia, x0=self.faces[i].c0, strikeRad=self.faces[i].str, dipRad=self.faces[i].dip, h=0.01*r)]
                #vtk file
                f_col = [f_0,f_1,f_2,f_3,f_4,f_5]
                sg.writeVtk(f_obj, f_col, f_lab, vtkFile=(fname + '_fracs.vtk'))
        
        #******   paint flowing fractures   ******
        if vtype[2] and len(self.faces)>6:
            q_obj = [] #fractures
            q_col = [] #fractures colors
            q_lab = [] #fractures color labels
            q_lab = ['Face_Number','Node_Number','Type','Mechanical_mm','Dilation_mm','Hydraulic_mm','Sn_MPa','Pcen_MPa','Pc_MPa','stim','Pmax_MPa','Tau_MPa','Mwmax','Propped_mm','Prop_m3','Conductivity_mDft']
            q_0 = []
            q_1 = []
            q_2 = []
            q_3 = []
            q_4 = []
            q_5 = []
            q_6 = []
            q_7 = []
            q_8 = []
            q_9 = []
            q_10 = []
            q_11 = []
            q_12 = []
            q_13 = []
            q_14 = []
            q_15 = []
            #nodex = np.asarray(self.nodes)
            for i in range(6,len(self.faces)): #skip boundary node at np.inf
                if self.faces[i].ci >= 0:
                    #add colors
                    q_0 += [i]
                    q_1 += [self.faces[i].ci]
                    q_2 += [self.faces[i].typ]
                    q_3 += [self.faces[i].bd*1000]
                    q_4 += [self.faces[i].bd0*1000]
                    q_5 += [self.faces[i].bh*1000]
                    q_6 += [self.faces[i].sn/MPa]
                    q_7 += [self.faces[i].Pcen]
                    q_8 += [self.faces[i].Pc/MPa]
                    q_9 += [self.faces[i].stim]
                    q_10 += [self.faces[i].Pmax]
                    q_11 += [self.faces[i].tau/MPa]
                    if len(self.faces[i].Mws) > 0:
                        q_12 += [np.max(self.faces[i].Mws)]
                    else:
                        q_12 += [-10.0]
                    q_13 +=  [self.faces[i].bd0p*1000]
                    q_14 +=  [self.faces[i].prop_load]
                    q_15 +=  [3.28084*1.01324e15*(self.faces[i].bh**3.0)/12.0]
                    #add geometry
                    q_obj += [HF(r=0.5*self.faces[i].dia, x0=self.faces[i].c0, strikeRad=self.faces[i].str, dipRad=self.faces[i].dip, h=0.02*r)]
            #vtk file
            q_col = [q_0,q_1,q_2,q_3,q_4,q_5,q_6,q_7,q_8,q_9,q_10,q_11,q_12,q_13,q_14,q_15]
            sg.writeVtk(q_obj, q_col, q_lab, vtkFile=(fname + '_fnets.vtk'))
        
        #******   paint nodes   ******
        if vtype[3]:
            n_obj = [] #nodes
            n_col = [] #nodes colors
            n_lab = [] #nodes color labels
            n_lab = ['Node_Number','Node_Pressure_MPa','Node_Temperature_K','Node_Enthalpy_kJ/kg']
            n_0 = []
            n_1 = []
            n_2 = []
            n_3 = []
            #nodex = np.asarray(self.nodes)
            for i in range(1,self.nodes.num): #skip boundary node at np.inf
                #add colors
                n_0 += [i]
                n_1 += [self.nodes.p[i]/MPa]
                n_2 += [self.nodes.T[i]]
                n_3 += [self.nodes.h[i]]
                #add geometry
                n_obj += [sg.cylObj(x0=self.nodes.all[i]+np.asarray([0.0,0.0,-r]), x1=self.nodes.all[i]+np.asarray([0.0,0.0,r]), r=r)]
            #vtk file
            n_col = [n_0,n_1,n_2,n_3]
            sg.writeVtk(n_obj, n_col, n_lab, vtkFile=(fname + '_nodes.vtk'))

        #******   paint pipes   ******
        if vtype[4]:
            self.pipes.coords(self.nodes)
            p_obj = [] #pipes
            p_col = [] #pipes colors
            p_lab = [] #pipes color labels
            p_lab = ['Pipe_Number','Type','Pipe_Flow_Rate_m3_s','Height_m','Length_m','Hydraulic_Aperture_mm','Max_Aperture_mm','Friction','ThermRadius_m','Parent_ID']
            p_0 = []
            p_1 = []
            p_2 = []
            p_3 = []
            p_4 = []
            p_5 = []
            p_6 = []
            p_7 = []
            p_8 = []
            p_9 = []
            # qs = np.asarray(self.q)
            # if not qs.any():
            #     qs = np.zeros(self.pipes.num)
            # qs = np.abs(qs)
            for i in range(0,self.pipes.num):
                #add geometry
                x0 = self.nodes.all[self.pipes.n0[i]]
                x1 = self.nodes.all[self.pipes.n1[i]]
                #don't include boundary node
                if not(np.isinf(x0[0]) or np.isinf(x1[0])):
                    p_obj += [sg.cylObj(x0=x0, x1=x1, r=0.666*r)]            
                    #add colors
                    p_0 += [i]
                    p_1 += [self.pipes.typ[i]]
                    p_2 += [self.pipes.q[i]]
                    p_3 += [self.pipes.W[i]]
                    p_4 += [self.pipes.L[i]]
                    p_5 += [self.pipes.Dh[i]*1000]
                    p_6 += [self.pipes.Dh_max[i]*1000]
                    p_7 += [self.pipes.frict[i]]
                    p_8 += [self.pipes.R0[i]]
                    p_9 += [self.pipes.pID[i]]
            #vtk file
            p_col = [p_0,p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9]
            sg.writeVtk(p_obj, p_col, p_lab, vtkFile=(fname + '_flow.vtk'))

        #******   paint boundaries   ******
        if vtype[5]:
            f_obj = [] #fractures
            f_col = [] #fractures colors
            f_lab = [] #fractures color labels
            f_lab = ['Face_Number','Node_Number','Type','Sn_MPa','Pc_MPa','Tau_MPa']
            f_0 = []
            f_1 = []
            f_2 = []
            f_3 = []
            f_4 = []
            f_5 = []
            #nodex = np.asarray(self.nodes)
            if len(self.faces) >= 6:
                for i in range(0,6): #skip boundary node at np.inf
                    #add colors
                    f_0 += [i]
                    f_1 += [self.faces[i].ci]
                    f_2 += [self.faces[i].typ]
                    f_3 += [self.faces[i].sn/MPa]
                    f_4 += [self.faces[i].Pc/MPa]
                    f_5 += [self.faces[i].tau/MPa]
                    #add geometry
                    f_obj += [HF(r=0.5*self.faces[i].dia, x0=self.faces[i].c0, strikeRad=self.faces[i].str, dipRad=self.faces[i].dip, h=0.01*r)]
                #vtk file
                f_col = [f_0,f_1,f_2,f_3,f_4,f_5]
                sg.writeVtk(f_obj, f_col, f_lab, vtkFile=(fname + '_bounds.vtk'))

    def build_pts(self,spacing=25.0,fname='test_gridx'):
        print( '*** constructing temperature grid ***')
        #structured grid of datapoints
        fname = fname + '_therm.vtk'
        size = self.rock.size
        num = int(2.0*size/spacing)+1
        label = 'temp_K'
        ns = [num,num,num]
        o0 = [-size,-size,-size]
        ss = [spacing,spacing,spacing]
        
        #initialize data to initial rock temperature
        #data = np.ones((num,num,num),dtype=float)*self.rock.BH_T
        data = np.ones((num,num,num),dtype=float)*self.rock.BH_T
        if self.rock.gradient:
            for i in range(0,num):
                data[:,:,i] = data[:,:,i] + (size - spacing*i) * self.rock.ResGradient
        
        #seek temperature drawdown
        for i in range(0,self.pipes.num):
            print( 'pipe %i' %(i))
            #collect fracture parameters
            x0 = self.nodes.all[self.pipes.n0[i]]
            x1 = self.nodes.all[self.pipes.n1[i]]
            T0 = self.nodes.T[self.pipes.n0[i]]
            T1 = self.nodes.T[self.pipes.n1[i]]
            c0 = 0.5*(x0+x1)
            r0 = 0.5*np.linalg.norm(x1-x0)
            #R0 = self.R0[i]
            R0 = self.pipes.R0[i]
            Tr0 = self.nodes.Tr[self.pipes.n0[i]]
            Tr1 = self.nodes.Tr[self.pipes.n1[i]]
            #pipes and wells
            if (int(self.pipes.typ[i]) in [typ('injector'),typ('producer'),typ('pipe'),typ('perfcluster'),typ('screen')]): #pipe, Hazen-Williams
                #well info
                dip = self.wells[self.pipes.fID[i]].dip #dip
                azn = self.wells[self.pipes.fID[i]].azn #azimuth
                vAxi = np.asarray([math.sin(azn)*math.cos(-dip),math.cos(azn)*math.cos(-dip),math.sin(-dip)]) #axial vector
                rc = self.wells[self.pipes.fID[i]].rc #azimuth
                #cycle thruogh all points
                for x in range(0,len(data)):
                    for y in range(0,len(data[0])):
                        for z in range(0,len(data[0,0])):
                            #point coordinates
                            xPt = np.asarray([o0[0]+x*spacing, o0[1]+y*spacing, o0[2]+z*spacing])
                            #spherical radius vector
                            pi = xPt-c0
                            #lengthwise distance along well
                            li = np.linalg.norm(np.dot(pi,vAxi))
                            #radial distance from well
                            ri = ((np.linalg.norm(pi))**2.0 - li**2.0)**0.5
                            #if within length and thermal radius
                            if (li < r0) and (ri < R0) and (ri > rc):
                                #subtract delta T at point based on distance of point versus thermal radius
                                data[x,y,z] = data[x,y,z] + (0.5*(T1+T0-Tr1-Tr0))*(1.0-(np.log(rc/ri)/np.log(rc/R0)))
                            elif (li < r0) and (ri <= rc):
                                data[x,y,z] = data[x,y,z] + (0.5*(T1+T0-Tr1-Tr0))
                
            #fractures and planes
            elif (int(self.pipes.typ[i]) in [typ('fracture'),typ('propped'),typ('choke')]): #(int())Y[i][2] == 1: #fracture, effective cubic law
                #fracture info
                dip = self.faces[self.pipes.fID[i]].dip
                azn = self.faces[self.pipes.fID[i]].str
                vNor = np.asarray([math.sin(azn+90.0*deg)*math.sin(dip),math.cos(azn+90.0*deg)*math.sin(dip),math.cos(dip)])
                vLeg = (x1-c0)/np.linalg.norm(x1-c0)
                vWid = np.cross(vNor,vLeg) 
                #cycle thruogh all points
                for x in range(0,len(data)):
                    for y in range(0,len(data[0])):
                        for z in range(0,len(data[0,0])):
                            #point coordinates
                            xPt = np.asarray([o0[0]+x*spacing, o0[1]+y*spacing, o0[2]+z*spacing])
                            #spherical radius vector
                            pi = xPt-c0
                            #normal distance from fracture
                            ni = np.linalg.norm(np.dot(pi,vNor))
                            #lengthwise distance from fracture
                            li = np.linalg.norm(np.dot(pi,vLeg))
                            #widthwise distance from fracture
                            wi = np.linalg.norm(np.dot(pi,vWid))
                            #if within length, width, normal
                            if (ni <= R0) and (li <= r0) and (wi <= (0.5*self.pipes.W[i])):
                                #subtract delta T at point based on distance of point versus thermal radius
                                data[x,y,z] = data[x,y,z] + (0.5*(T1+T0)-self.rock.BH_T)*(1.0-ni/R0)

        head = '# vtk DataFile Version 2.0\n'
        head += 'pointcloud\n'
        head += 'ASCII\n'
        head += 'DATASET STRUCTURED_POINTS\n'
        head += 'DIMENSIONS %i %i %i\n' %(ns[0],ns[1],ns[2])
        head += 'ORIGIN %f %f %f\n' %(o0[0],o0[1],o0[2])
        head += 'SPACING %f %f %f\n' %(ss[0],ss[1],ss[2])
        head += 'POINT_DATA %i\n' %(ns[0]*ns[1]*ns[2])
        head += 'SCALARS ' + label + ' float 1\n'
        head += 'LOOKUP_TABLE default'
        
        print( head)
        
        try:
            with open(fname,'r') as f:
                test = f.readline()
            f.close()
            if test != '':
                print( 'file already exists')
                f.close()
            else:
                with open(fname,'a') as f:
                    f.write(head + '\n')
                    out = ''
                    for k in range(0,len(data[0,0,:])):
                        for j in range(0,len(data[0,:,0])):
                            for i in range(0,len(data[:,0,0])):
                                out += '%e' %(data[i,j,k]) + '\n'
                    f.write(out)
                f.close()
        except:
            with open(fname,'a') as f:
                f.write(head + '\n')
                out = ''
                for k in range(0,len(data[0,0,:])):
                    for j in range(0,len(data[0,:,0])):
                        for i in range(0,len(data[:,0,0])):
                            out += '%e' %(data[i,j,k]) + '\n'
                f.write(out)
            f.close()


    def re_init(self):
        #clear prior data
        self.nodes = []
        self.nodes = nodes()
        self.pipes = []
        self.pipes = pipes()
        self.faces = [] #all surfaces
        self.faces = self.bound + self.fracs + self.hydfs
        for i in range(0,len(self.faces)):
            self.faces[i].ci = -1
        self.trakr = []
        self.H = []
        self.Q = []
        self.nodes.p = self.nodes.p * 0.0
        self.nodes.T = self.nodes.T * 0.0
        self.nodes.h = self.nodes.h * 0.0
        self.q = []
        return self

    def set_bcs(self,p_bound=0.0*MPa,q_well=[],p_well=[],auto_farfield=True):
        #working variables
        rho = self.rock.PoreRho #kg/m3
        #input & output boundaries
        self.H = []
        self.Q = []
        #outer boundary (always zero index node)
        self.H += [[0, p_bound/(rho*g)]]
        #well boundaries (lists with empty elements)
        #q_well = [None,  None, 0.02,  None]
        #p_well = [None, 3.0e6, None, 1.0e6]
        #default = p_bound
        for w in range(0,len(self.wells)):
            #identify boundary node
            #coordinates
            source = self.wells[w].c0
            #find index of duplicate
            ck, i = self.nodes.add(source)
            if not(ck): #yes duplicate
                #prioritize flow boundary conditions
                if q_well[w] != None:
                    self.Q += [[i, -q_well[w]]]
                #pressure boundary conditions
                elif p_well[w] != None:
                    self.H += [[i, p_well[w]/(rho*g)]]
                #default to outer boundary condition if no boundary condition is explicitly stated
                elif auto_farfield:
                    self.H += [[i, p_bound/(rho*g)]]
                    print( 'warning: a well was assigned the far-field pressure boundary condition')
                else:
                    #do nothing
                    pass
            else:
                print( 'error: flow boundary point not identified')

    def therm_bcs(self,T_bound=0.0,T_inlet=[]):
        #outer boundary (always zero index node)
        self.Tb += [[0, T_bound]]
        
        #boundary conditions from wells
        n = 0
        for w in range(0,len(self.wells)):
            #temperature boundary condition
            if (self.wells[w].typ == typ('injector')):
                #coordinates
                source = self.wells[w].c0
                
                #find index of duplicate
                ck, i = self.nodes.add(source)
                if not(ck): #yes duplicate
                    if len(T_inlet)>1:
                        self.Tb += [[i, T_inlet[n]]]
                        n += 1
                    else:
                        self.Tb += [[i, T_inlet[0]]]
                else:
                    print( 'error: pressure temperature point not identified')
                
    def add_flowpath(self, source, target, length, width, featTyp, featID, Dh=1.0, Dh_max=1.0, frict=1.0, pID = 0):
#        #ignore if source == target
#        if list(source) == list(target):
#            return -1, -1
        #source node
        n_s, so_n = self.nodes.add(source)#,f_id=featID)
        #target node
        n_t, ta_n = self.nodes.add(target)#,f_id=featID)
        #check if source is same as target
        if so_n == ta_n:
            return -1, -1
        #check if the reversed node set already exists (i.e., don't create pipes forward and backward between same nodes)
        if not(n_s) and not(n_t):
            to_ck = np.where(np.asarray(self.pipes.n0) == ta_n)[0]
            if len(to_ck) > 0:
                for i in to_ck:
                    if self.pipes.n1[i] == so_n:
                        return -1, -1
        #add pipe
        self.pipes.add(so_n, ta_n, length, width, featTyp, featID, Dh, Dh_max, frict, pID)
        return so_n, ta_n
        
    #intersections of a line with a plane
    def x_well_all_faces(self,plot=True,sourceID=0,targetID=[],offset=[]): #, path_type=0, aperture=0.22, roughness=80.0): #[x0,y0,zo,len,azn,dip]
        #scaled visual offset
        if not(offset):
            offset = np.max([5.0*self.nodes.tol,0.010*self.rock.size])*np.asarray([0.58,0.58,0.58])
            
        #working array for finding and logging intersection points
        x_well = [] #intercept coord
        o_frac = [] #surface origin
        r_frac = [] #index of well
        i_frac = [] #index of frac
        
        #line location
        c0 = self.wells[sourceID].c0 #line origin
        leg = self.wells[sourceID].leg #line length
        azn = self.wells[sourceID].azn #line azimuth
        dip = self.wells[sourceID].dip #line dip
        dia = 2.0*self.wells[sourceID].ra #line inner diameter
        lty = self.wells[sourceID].typ #line type
        vAxi = np.asarray([math.sin(azn)*math.cos(-dip),math.cos(azn)*math.cos(-dip),math.sin(-dip)])
        cm = c0 + 0.5*leg*vAxi #line midpoint
        c1 = c0 + leg*vAxi #line endpoint
        
        #skip intersections for sealed pipe
        if (int(lty) in [typ('pipe')]):
            #add endpoint to endpoint
            self.add_flowpath(c0,
                             c1,
                             leg,
                             dia,
                             lty,
                             sourceID,
                             Dh=dia,
                             Dh_max=dia,
                             frict=self.wells[sourceID].rgh,
                             pID=self.wells[sourceID].pID)
            return False
        #for all target faces (except boundaries)
        for targetID in range(6,len(self.faces)): #!!! changed with gradient
            #planar object parameters
            t0 = self.faces[targetID].c0 #face origin
            rad = 0.5*self.faces[targetID].dia #face diameter
            azn = self.faces[targetID].str #face strike
            dip = self.faces[targetID].dip #face dip
            fty = self.faces[targetID].typ #face type
            vNor = np.asarray([math.sin(azn+90.0*deg)*math.sin(dip),math.cos(azn+90.0*deg)*math.sin(dip),math.cos(dip)])
            #infinite plane intersection point
            if np.dot(vNor,vAxi) != 0: #not parallel
                x_test = c0 + vAxi*(np.dot(vNor,t0) - np.dot(vNor,c0))/(np.dot(vNor,vAxi))
            
                # test for intersect within plane and line extents
                if (np.linalg.norm(cm-x_test) < (0.5*leg)) and (np.linalg.norm(t0-x_test) < (rad)):
                    x_well += [x_test]
                    o_frac += [t0]
                    r_frac += [rad] 
                    i_frac += [targetID]
#                    fs_i += [targetID]
                    self.trakr += [[-1,targetID]]
        #in case of no intersections
        #add_flowpath(self, source, target, length, width, featTyp, featID, tol = 0.001):
        if len(x_well) == 0:
            #add endpoint to endpoint
            self.add_flowpath(c0,
                             c1,
                             leg,
                             dia,
                             lty,
                             sourceID,
                             Dh=dia,
                             Dh_max=dia,
                             frict=self.wells[sourceID].rgh,
                             pID=self.wells[sourceID].pID)
            #return False
            return False
        #in case of intersections
        else:
            #convert to array
            x_well = np.asarray(x_well)
            
            #sort intersections by distance from origin point
            rs = []
            rs = np.linalg.norm(x_well-c0,axis=1)
            a = rs.argsort()
            #first element well (a live end)
            self.add_flowpath(c0,
                             x_well[a[0]] + offset,
                             rs[a[0]],
                             dia,
                             lty,
                             sourceID,
                             Dh=dia,
                             Dh_max=dia,
                             frict=self.wells[sourceID].rgh,
                             pID=self.wells[sourceID].pID)
            #intersection points
            i = 0
            for i in range(0,len(rs)-1):
                #cased well = perforate
                if (int(lty) in [typ('perfcluster')]):
                    #perforation 
                    self.add_flowpath(x_well[a[i]] + offset,
                                     x_well[a[i]] + 0.5*offset,
                                     self.wells[sourceID].rc-self.wells[sourceID].ra,
                                     self.rock.perf_dia,
                                     typ('perf'),
                                     #sourceID,
                                     i_frac[a[i]],
                                     Dh=self.rock.perf_dia,
                                     Dh_max=self.rock.perf_dia,
                                     frict=self.rock.perf_per_cluster,
                                     pID=self.wells[sourceID].pID)
                    #choke
                    self.add_flowpath(x_well[a[i]] + 0.5*offset,
                                     x_well[a[i]],
                                     3.0*self.wells[sourceID].rc,
                                     math.pi*self.wells[sourceID].rc,
                                     typ('choke'),
                                     i_frac[a[i]],
                                     Dh=self.faces[i_frac[a[i]]].bh,
                                     Dh_max=self.faces[i_frac[a[i]]].bh,
                                     frict=self.faces[i_frac[a[i]]].roughness,
                                     pID=self.wells[sourceID].pID)
                #uncased = no perforation
                else:
                    #choke (circumference of well * 3.0 * diameter = near well flow channel area dimensions, otherwise properties of the fracture)
                    self.add_flowpath(x_well[a[i]] + offset,
                                     x_well[a[i]],
                                     3.0*self.wells[sourceID].rc,
                                     math.pi*self.wells[sourceID].rc,
                                     typ('choke'),
                                     i_frac[a[i]],
                                     Dh=self.faces[i_frac[a[i]]].bh,
                                     Dh_max=self.faces[i_frac[a[i]]].bh,
                                     frict=self.faces[i_frac[a[i]]].roughness,
                                     pID=self.wells[sourceID].pID)
                #fracture (use intercept to center length, but fix width to y at 1/2 cirle radius)
                self.add_flowpath(x_well[a[i]],
                                 o_frac[a[i]],
                                 np.linalg.norm(o_frac[a[i]]-(x_well[a[i]])), #+offset*0.5)),
                                 0.866*r_frac[a[i]],
                                 fty,
                                 i_frac[a[i]],
                                 Dh=self.faces[i_frac[a[i]]].bh,
                                 Dh_max=self.faces[i_frac[a[i]]].bh,
                                 frict=self.faces[i_frac[a[i]]].roughness,
                                 pID=self.wells[sourceID].pID)
                #well continued (+1.0 z offset to prevent non-real links from fracture to well without a choke)
                self.add_flowpath(x_well[a[i]] + offset,
                                  x_well[a[i+1]] + offset,
                                  rs[a[i+1]]-rs[a[i]],
                                  dia,
                                  lty,
                                  sourceID,
                                  Dh=dia,
                                  Dh_max=dia,
                                  frict=self.wells[sourceID].rgh,
                                  pID=self.wells[sourceID].pID)
                #store fracture centerpoint node number
                ck, cki = self.nodes.add(o_frac[a[i]])
                self.faces[i_frac[a[i]]].ci = cki
                
            #last segment choke
            #cased well = perforate
            if (int(lty) in [typ('perfcluster')]):
                #perforation 
                self.add_flowpath(x_well[a[-1]] + offset,
                                 x_well[a[-1]] + 0.5*offset,
                                 self.wells[sourceID].rc-self.wells[sourceID].ra,
                                 self.rock.perf_dia,
                                 typ('perf'),
                                 #sourceID,
                                 i_frac[a[i]],
                                 Dh=self.rock.perf_dia,
                                 Dh_max=self.rock.perf_dia,
                                 frict=self.rock.perf_per_cluster,
                                 pID=self.wells[sourceID].pID)
                #choke
                self.add_flowpath(x_well[a[-1]] + 0.5*offset,
                                 x_well[a[-1]],
                                 3.0*self.wells[sourceID].rc,
                                 math.pi*self.wells[sourceID].rc,
                                 typ('choke'),
                                 i_frac[a[-1]],
                                 Dh=self.faces[i_frac[a[-1]]].bh,
                                 Dh_max=self.faces[i_frac[a[-1]]].bh,
                                 frict=self.faces[i_frac[a[-1]]].roughness,
                                 pID=self.wells[sourceID].pID)
            #uncased = no perforation
            else:
                #choke
                self.add_flowpath(x_well[a[-1]] + offset,
                                 x_well[a[-1]],
                                 3.0*self.wells[sourceID].rc,
                                 math.pi*self.wells[sourceID].rc,
                                 typ('choke'),
                                 i_frac[a[-1]],
                                 Dh=self.faces[i_frac[a[-1]]].bh,
                                 Dh_max=self.faces[i_frac[a[-1]]].bh,
                                 frict=self.faces[i_frac[a[-1]]].roughness,
                                 pID=self.wells[sourceID].pID)
            #last segment fracture
            self.add_flowpath(x_well[a[-1]],# + offset*0.5,
                             o_frac[a[-1]],
                             np.linalg.norm(o_frac[a[-1]]-(x_well[a[-1]])), #+offset*0.5)),
                             0.866*r_frac[a[-1]],
                             fty,
                             i_frac[a[i]],
                             Dh=self.faces[i_frac[a[-1]]].bh,
                             Dh_max=self.faces[i_frac[a[-1]]].bh,
                             frict=self.faces[i_frac[a[-1]]].roughness,
                             pID=self.wells[sourceID].pID)
            #dead end segment
            self.add_flowpath(x_well[a[-1]] + offset,
                              c1,
                              np.linalg.norm(x_well[a[-1]]-c1), #,axis=1),
                              dia,
                              lty,
                              sourceID,
                              Dh=dia,
                              Dh_max=dia,
                              frict=self.wells[sourceID].rgh,
                              pID=self.wells[sourceID].pID)
            #store fracture centerpoint node number
            ck, cki = self.nodes.add(o_frac[a[i]])
            self.faces[i_frac[a[-1]]].ci = cki
            return True
    
    #intersections of a plane with a plane
    def x_frac_face(self,plot=True,sourceID=0,targetID=1): #[x0,y0,z0,dia,azn,dip]
        #plane 1
        dia1 = 0.5*self.faces[sourceID].dia
        dip = self.faces[sourceID].dip
        azn = self.faces[sourceID].str
        vNor1 = np.asarray([math.sin(azn+90.0*deg)*math.sin(dip),math.cos(azn+90.0*deg)*math.sin(dip),math.cos(dip)])
        c01 = self.faces[sourceID].c0
        f1_t = self.faces[sourceID].typ
        
        #plane 2
        dia2 = 0.5*self.faces[targetID].dia
        dip = self.faces[targetID].dip
        azn = self.faces[targetID].str
        vNor2 = np.asarray([math.sin(azn+90.0*deg)*math.sin(dip),math.cos(azn+90.0*deg)*math.sin(dip),math.cos(dip)])
        c02 = self.faces[targetID].c0
        f2_t = self.faces[targetID].typ
        
        #intersection vector
        vInt = []
        vInt = np.cross(vNor1,vNor2)
        
        #if not parallel
        if np.dot(vNor1,vNor2) < 0.999999:
            #intersection vector origin point
            zero = np.argmax(np.abs(np.asarray(vInt)),axis=0)
            d1 = -1*(vNor1[0]*c01[0] + vNor1[1]*c01[1] + vNor1[2]*c01[2])
            d2 = -1*(vNor2[0]*c02[0] + vNor2[1]*c02[1] + vNor2[2]*c02[2])
            vN1 = np.delete(vNor1, zero, axis=0)
            vN2 = np.delete(vNor2, zero, axis=0)
            cInt = (np.asarray([np.linalg.det([[vN1[1],vN2[1]],[d1,d2]]),
                   np.linalg.det([[d1,d2],[vN1[0],vN2[0]]])]) 
                   / np.linalg.det([[vN1[0],vN2[0]],[vN1[1],vN2[1]]]))
            cInt = np.insert(cInt,zero,0.0,axis=0)

            #endpoints - plane 1 intersection
            c = (cInt[0]-c01[0])**2.0 + (cInt[1]-c01[1])**2.0 + (cInt[2]-c01[2])**2.0 - dia1**2.0
            b = 2.0*( (cInt[0]-c01[0])*vInt[0] + (cInt[1]-c01[1])*vInt[1] + (cInt[2]-c01[2])*vInt[2] )
            a = vInt[0]**2.0 + vInt[1]**2.0 + vInt[2]**2.0
            if (b**2.0-4.0*a*c) >= 0.0:
                l1 = (-b+(b**2.0-4.0*a*c)**0.5)/(2.0*a)
                l2 = (-b-(b**2.0-4.0*a*c)**0.5)/(2.0*a)
                f1a = cInt + l1 * vInt
                f1b = cInt + l2 * vInt
            else:
                f1a = np.asarray([np.nan,np.nan,np.nan])
                f1b = np.asarray([np.nan,np.nan,np.nan])            
            
            #endpoints - plane 2 intersection
            c = (cInt[0]-c02[0])**2.0 + (cInt[1]-c02[1])**2.0 + (cInt[2]-c02[2])**2.0 - dia2**2.0
            b = 2.0*( (cInt[0]-c02[0])*vInt[0] + (cInt[1]-c02[1])*vInt[1] + (cInt[2]-c02[2])*vInt[2] )
            a = vInt[0]**2.0 + vInt[1]**2.0 + vInt[2]**2.0
            if (b**2.0-4.0*a*c) >= 0.0:
                l1 = (-b+(b**2.0-4.0*a*c)**0.5)/(2.0*a)
                l2 = (-b-(b**2.0-4.0*a*c)**0.5)/(2.0*a)
                f2a = cInt + l1 * vInt
                f2b = cInt + l2 * vInt
            else:
                f2a = np.asarray([np.nan,np.nan,np.nan])
                f2b = np.asarray([np.nan,np.nan,np.nan])
            
            #midpoint
            xInt = np.asarray([f1a,f1b,f2a,f2b])
            xInt = np.unique(xInt,axis=0)
            slot = len(xInt)-1
            xMid = []
            for i in range(0,slot+1):
                if (np.linalg.norm(xInt[slot-i]-c01) >= 1.01*(dia1)) or (np.linalg.norm(xInt[slot-i]-c02) >= 1.01*(dia2)) or (np.sum(np.isnan(xInt)) > 0):
                    xInt = np.delete(xInt,(slot-i),axis=0)
            if len(xInt) == 2 and np.sum(np.isnan(xInt)) == 0:
                xMid = 0.5*(xInt[0]+xInt[1])
            #add pipes to network            
            if (np.sum(np.isnan(xInt)) == 0) and len(xInt) == 2:
                #normal fracture-fracture connection
                if (f1_t != typ('boundary')) and (f2_t != typ('boundary')):
                    #source-center to intersection midpoint
                    #add_flowpath(self, source, target, length, width, featTyp, featID):
                    self.add_flowpath(c01,
                                     xMid,
                                     np.linalg.norm(xMid-c01),
                                     np.linalg.norm(xInt[1]-xInt[0]),
                                     f1_t,
                                     sourceID,
                                     Dh=self.faces[sourceID].bh,
                                     Dh_max=self.faces[sourceID].bh,
                                     frict=self.faces[sourceID].roughness)
                    #intersection midpoint to target-center
                    p_1, p_2 = self.add_flowpath(xMid,
                                     c02,
                                     np.linalg.norm(xMid-c02),
                                     np.linalg.norm(xInt[1]-xInt[0]),
                                     f2_t,
                                     targetID,
                                     Dh=self.faces[targetID].bh,
                                     Dh_max=self.faces[targetID].bh,
                                     frict=self.faces[targetID].roughness)
                    #store fracture centerpoint node number
                    if p_2 >= 0:
                        self.faces[targetID].ci = p_2
                    
                    #update tracker
                    self.trakr += [[sourceID,targetID]] 
                #fracture-boundary connection (boundary type = -3)
                elif (f2_t == typ('boundary')) and (f1_t != typ('boundary')):
                    #source-center to intersection midpoint
                    self.add_flowpath(c01,
                                     xMid,
                                     np.linalg.norm(xMid-c01),
                                     np.linalg.norm(xInt[1]-xInt[0]),
                                     f1_t,
                                     sourceID,
                                     Dh=self.faces[sourceID].bh,
                                     Dh_max=self.faces[sourceID].bh,
                                     frict=self.faces[sourceID].roughness)
                    #intersection midpoint to far-field
                    p_1, p_2 = self.add_flowpath(xMid,
                                     self.nodes.r0,
                                     100.0*dia2,
                                     np.linalg.norm(xInt[1]-xInt[0]),
                                     f2_t,
                                     targetID,
                                     Dh=self.faces[targetID].bh,
                                     Dh_max=self.faces[targetID].bh,
                                     frict=self.faces[targetID].roughness)
                    #store fracture centerpoint node number
                    if p_2 >= 0:
                        self.faces[targetID].ci = p_2
                
    # ********************************************************************
    # domain creation
    # ********************************************************************
    def gen_domain(self,plot=True):
        print( '*** domain boundaries module ***')
        #clear old data
        self.bound = []
        #working variables
        size = self.rock.size
        #create boundary faces for analysis
        domb3D = []
        domb3D += [surf(-size,0.0,0.0,4.0*size,00.0*deg,90.0*deg,'boundary',self.rock)]
        domb3D += [surf(size,0.0,0.0,4.0*size,00.0*deg,90.0*deg,'boundary',self.rock)]
        domb3D += [surf(0.0,-size,0.0,4.0*size,90.0*deg,90.0*deg,'boundary',self.rock)]
        domb3D += [surf(0.0,size,0.0,4.0*size,90.0*deg,90.0*deg,'boundary',self.rock)]
        domb3D += [surf(0.0,0.0,-size,4.0*size,00.0*deg,00.0*deg,'boundary',self.rock)]
        domb3D += [surf(0.0,0.0,size,4.0*size,00.0*deg,00.0*deg,'boundary',self.rock)]
        #add to model domain
        self.bound = domb3D        

    # ********************************************************************
    # fracture network creation
    # ********************************************************************
    def gen_fixfrac(self,clear=True,
                     c0 = [0.0,0.0,0.0],
                     dia = 500.0,
                     azn = 80.0*deg,
                     dip = 90.0*deg):
        # print( '-> manual fracture placement module')
        #clear old data
        if clear:
            self.fracs = []
        
        #place fracture
        c0 = np.asarray(c0)
        
        #compile list of fractures
        frac3D = [surf(c0[0],c0[1],c0[2],dia,azn,dip,'fracture',self.rock)]
        
        # print( 'dia = %.1f, azn = %.1f, dip = %.1f' %(dia, azn, dip))

        #add to model domain
        self.fracs += frac3D
    
    # ********************************************************************
    # fracture network creation
    # ********************************************************************
    def gen_natfracs(self,clear=True,
                     f_num = 40,
                     f_dia = [200.0,900.0],
                     f_azn = [79.0*deg,8.0*deg],
                     f_dip = [90.0*deg,12.5*deg]):
        print( '-> dfn seed module')
        #clear old data
        if clear:
            self.fracs = []
        #size of domain for reservoir
        size = self.rock.size 
        #working variables
        frac3D = []
        #populate fractures
        for n in range(0,f_num):
            #Fracture parameters
            #dia = np.random.uniform(f_dia[0],f_dia[1])
            logmu = 0.5*(np.log10(f_dia[0])+np.log10(f_dia[1]))
            dia = lognorm_trunc(1,logmu,logmu,np.log10(f_dia[0]),np.log10(f_dia[1]))[0] #!!!
            azn = np.random.uniform(f_azn[0],f_azn[1])
            dip = np.random.uniform(f_dip[0],f_dip[1])
            #Build geometry
            x = np.random.uniform(-size,size)
            y = np.random.uniform(-size,size)
            z = np.random.uniform(-size,size)
            c0 = np.asarray([x,y,z])
            #compile list of fractures
            frac3D += [surf(c0[0],c0[1],c0[2],dia,azn,dip,'fracture',self.rock)]
            
        #add to model domain
        self.fracs += frac3D
        
    # ************************************************************************
    # deterministic fracture placement
    # ************************************************************************
    def gen_from_dfn(self,clear=False,filename=''):
        if filename == '':
            filename = self.rock.setup.Strategy_DiscreteFractures_path
        if filename != '':
            print( '*** read-in dfn: %s ***' %(filename))
            #clear old joint sets
            if clear:
                self.fracs = []
            #generate and parameterize fractures from dfn
            dfn().intake(filename=filename,fracs=self.fracs,rock=self.rock)
    
    # ************************************************************************
    # stochastic fracture placement
    # ************************************************************************
    def gen_joint_sets(self,clear=True):
        print( '*** joint set module ***')
        #clear old joint sets
        if clear: 
            self.fracs = []
        #generate fractures from sets
        for i in range(0,len(self.rock.fNum)):
            self.gen_natfracs(False,
                         self.rock.fNum[i],
                         self.rock.fDia[i],
                         self.rock.fStr[i],
                         self.rock.fDip[i])            

    # ************************************************************************
    # stimulation - static
    # ************************************************************************
    def gen_stimfracs(self,plot=True,
                      target=0,
                      stages=1,
                      perfs=1,
                      f_dia = [1266.0,107.0],
                      f_azn = [65.0*deg,4.0*deg],
                      f_dip = [57.5*deg,3.75*deg],
                      clear = True):
        #print( '*** stimulation module ***')
        #clear old data
        if clear == True:
            self.hydfs = []
        #hydraulic fractures
        num = stages*perfs
        #stimulate target well
        spa = self.wells[target].leg/(num+1)
        c0 = self.wells[target].c0
        azn = self.wells[target].azn
        dip = self.wells[target].dip
        vAxi = np.asarray([math.sin(azn)*math.cos(-dip),math.cos(azn)*math.cos(-dip),math.sin(-dip)])
        for n in range(0,num):
            #Fracture parameters
            leg = np.random.normal(0,1)*f_dia[1] + f_dia[0]
            azn = np.random.normal(0,1)*f_azn[1] + f_azn[0]
            dip = np.random.normal(0,1)*f_dip[1] + f_dip[0]
            c0 = c0 + spa*vAxi
            
            #add to model domain
            self.hydfs += [surf(c0[0],c0[1],c0[2],leg,azn,dip,'propped',self.rock, mcc = self.rock.hfmcc, phi = self.rock.hfphi)]
            self.hydfs[-1].bd0 = 0.0
    
    # ************************************************************************
    # well placement
    # ************************************************************************
    def gen_wells(self,clear=True,wells=[],style='EGS'):
        print( '*** well placement module ***')
        #clear old data
        if clear:
            self.wells = []
        #manual placement
        if len(wells) > 0:
            self.wells = wells
        #EGS style placement
        elif (style == 'EGS'):
            n = len(self.wells)
            wells = self.wells
            #center
            i0 = np.asarray([0.0, 0.0, 0.0])
            #ref axes (parallel injector and 90 horizontal to the right)
            azn = self.rock.w_azimuth
            dip = self.rock.w_dip
            vInj = np.asarray([math.sin(azn)*math.cos(-dip),
                     math.cos(azn)*math.cos(-dip),math.sin(-dip)])
            vRht = np.asarray([math.sin(azn+pi/2),math.cos(azn+pi/2),0.0])
            vNor = np.asarray([math.sin(azn)*math.sin(-dip),
                     math.cos(azn)*math.sin(-dip),math.cos(-dip)])
            #offset center
            p0 = i0 + vRht*self.rock.w_spacing
            #toe in production well
            toe = self.rock.w_toe
            vPro = sg.rotatePoints([vInj],vNor,toe)[0]
            #skew in production well
            skew = self.rock.w_skew
            vPro = sg.rotatePoints([vPro],vRht,skew)[0]
            #phase of production wells
            p0s = []
            vPros = []
            phase = self.rock.w_phase
            num = self.rock.w_count
            for i in range(0,num):
                p0s += [sg.rotatePoints([p0],vInj,(i*2.0*pi/num + phase))[0]]
                vPros += [sg.rotatePoints([vPro],vInj,(i*2.0*pi/num + phase))[0]]
            #lengths of wells
            length = self.rock.w_length
            proportion = self.rock.w_proportion
            iLen = length*proportion
            pLen = length
            cLen = pLen - iLen
            #injection well reservoir segments
            i1s = [] #startpoints
            i2s = [] #endpoints
            seg = self.rock.w_intervals
            leg = iLen/(seg+seg-1)
            i1s += [i0 - 0.5*vInj*iLen]
            i2s += [i1s[0] + leg*vInj]
            for i in range(1,seg):
                i1s += [i2s[i-1] + leg*vInj]
                i2s += [i1s[i] + leg*vInj]
            #injection well segments to surface
            ibs = i1s[0] - 0.5*cLen*vInj
            iss = np.asarray([ibs[0],ibs[1],self.rock.ResDepth])
            #production well segments
            p1s = [] #startpoints
            p2s = [] #endpoints
            for i in range(0,num):
                p1s += [p0s[i] - 0.5*vPros[i]*pLen]
                p2s += [p0s[i] + 0.5*vPros[i]*pLen]    
            #production well segments to surface
            pss = []
            for i in range(0,num):
                pss += [np.asarray([p1s[i][0],p1s[i][1],self.rock.ResDepth+0.1])]
            #wellheads
            h = np.asarray([0.0,0.0,self.rock.w_spacing])
            azn, dip, dis = azn_dip(iss+h,iss)
            wells += [line(iss[0]+h[0],iss[1]+h[1],iss[2]+h[2],dis,azn,dip,'injector',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
            for i in range(0,num):
                azn, dip, dis = azn_dip(pss[i]+h,pss[i])
                wells += [line(pss[i][0]+h[0],pss[i][1]+h[1],pss[i][2]+h[2],dis,azn,dip,'producer',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+i+n)]
            #place surface segments
            azn, dip, dis = azn_dip(iss,ibs)
            wells += [line(iss[0],iss[1],iss[2],dis,azn,dip,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
            for i in range(0,num):
                azn, dip, dis = azn_dip(pss[i],p1s[i])
                wells += [line(pss[i][0],pss[i][1],pss[i][2],dis,azn,dip,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+i+n)]
            #place injector base segment
            azn, dip, dis = azn_dip(ibs,i1s[0])
            wells += [line(ibs[0],ibs[1],ibs[2],dis,azn,dip,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
            #place injection wells
            for i in range(0,seg-1):
                azn, dip, dis = azn_dip(i1s[0],i2s[0])
                wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],dis,azn,dip,'perfcluster',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
                azn, dip, dis2 = azn_dip(i2s[0],i1s[1])
                wells += [line(i2s[i][0],i2s[i][1],i2s[i][2],dis2,azn,dip,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
            azn, dip, dis = azn_dip(i1s[-1],i2s[-1])
            wells += [line(i1s[-1][0],i1s[-1][1],i1s[-1][2],dis,azn,dip,'perfcluster',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
            #place production wells
            for i in range(0,num):
                azn, dip, dis = azn_dip(p1s[i],p2s[i])
                wells += [line(p1s[i][0],p1s[i][1],p1s[i][2],dis,azn,dip,'screen',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+i+n)]
            #add to model domain
            self.wells = wells
        #EGS line placement with wine rack
        elif (style == 'O&G'):
            n = len(self.wells)
            wells = self.wells
            #center
            i0 = np.asarray([0.0, 0.0, 0.0])
            #ref axes (parallel injector and 90 horizontal to the right)
            azn = self.rock.w_azimuth
            dip = self.rock.w_dip
            vInj = np.asarray([math.sin(azn)*math.cos(-dip),
                     math.cos(azn)*math.cos(-dip),math.sin(-dip)])
            vRht = np.asarray([math.sin(azn+pi/2),math.cos(azn+pi/2),0.0])
            vNor = np.asarray([math.sin(azn)*math.sin(-dip),
                     math.cos(azn)*math.sin(-dip),math.cos(-dip)])
            #offset center
            p0 = i0 + vRht*self.rock.w_spacing
            #toe in production well
            toe = self.rock.w_toe
            vPro = sg.rotatePoints([vInj],vNor,toe)[0]
            #skew in production well
            skew = self.rock.w_skew
            vPro = sg.rotatePoints([vPro],vRht,skew)[0]
            vPro = vInj
            #populate rack of production and injection wells
            i0s = []
            p0s = []
            phase = self.rock.w_phase
            num = self.rock.w_count
            for i in range(0,num):
                if i == 0: 
                    #first production well includes a mirrored well
                    p0s += [sg.rotatePoints([p0],vInj,(pi - phase))[0]]
                    i0s += [i0]
                else: 
                    #place additional producers
                    p0s += [p0s[-1] + 2.0*vRht*self.rock.w_spacing*math.cos(phase)]
                #don't place injector in second step
                if i > 1:
                    i0s += [i0s[-1] + 2.0*vRht*self.rock.w_spacing*math.cos(phase)]
            #recenter
            c0 = np.zeros(3)
            if num < 2:
                c0 = c0 - 0.5*vRht*self.rock.w_spacing*math.cos(phase)
            else:
                c0 = c0 + (num-2)*vRht*self.rock.w_spacing*math.cos(phase)
            for i in range(0,len(i0s)):
                i0s[i] += -1.0*c0
            for i in range(0,len(p0s)):
                p0s[i] += -1.0*c0  
            #lengths of wells
            length = self.rock.w_length
            proportion = self.rock.w_proportion
            iLen = length*proportion
            pLen = length
            cLen = pLen - iLen
            seg = self.rock.w_intervals
            leg = iLen/(seg+seg-1)
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
                iss = np.asarray([ibs[0],ibs[1],self.rock.ResDepth])
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
                pss += [np.asarray([p1s[i][0],p1s[i][1],self.rock.ResDepth+0.1])]
            #place wellheads
            h = np.asarray([0.0,0.0,self.rock.w_spacing])
            iid = np.asarray(range(0,len(i0s)))+n
            pid = np.asarray(range(0,len(p0s)))+n+len(i0s)
            for i in range(0,len(i0s)):
                iss = isss[i]
                azn, dip, dis = azn_dip(iss+h,iss)
                wells += [line(iss[0]+h[0],iss[1]+h[1],iss[2]+h[2],dis,azn,dip,'injector',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,iid[i])]
            for i in range(0,len(p0s)):
                azn, dip, dis = azn_dip(pss[i]+h,pss[i])
                wells += [line(pss[i][0]+h[0],pss[i][1]+h[1],pss[i][2]+h[2],dis,azn,dip,'producer',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,pid[i])]
            #place surface segments
            for i in range(0,len(i0s)):
                iss = isss[i]
                ibs = ibss[i]
                azn, dip, dis = azn_dip(iss,ibs)
                wells += [line(iss[0],iss[1],iss[2],dis,azn,dip,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,iid[i])]
            for i in range(0,len(p0s)):
                azn, dip, dis = azn_dip(pss[i],p1s[i])
                wells += [line(pss[i][0],pss[i][1],pss[i][2],dis,azn,dip,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,pid[i])]
            #place injector base segment
            for i in range(0,len(i0s)):
                ibs = ibss[i]
                i1s = i1ss[i]
                i2s = i2ss[i]
                azn, dip, dis = azn_dip(ibs,i1s[0])
                wells += [line(ibs[0],ibs[1],ibs[2],dis,azn,dip,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,iid[i])]
                #place injection internal segments
                for i in range(0,seg-1):
                    azn, dip, dis = azn_dip(i1s[0],i2s[0])
                    wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],dis,azn,dip,'perfcluster',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,iid[i])]
                    azn, dip, dis2 = azn_dip(i2s[0],i1s[1])
                    wells += [line(i2s[i][0],i2s[i][1],i2s[i][2],dis2,azn,dip,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,iid[i])]
                azn, dip, dis = azn_dip(i1s[-1],i2s[-1])
                wells += [line(i1s[-1][0],i1s[-1][1],i1s[-1][2],dis,azn,dip,'perfcluster',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,iid[i])]
            #place production wells
            for i in range(0,len(p0s)):
                azn, dip, dis = azn_dip(p1s[i],p2s[i])
                wells += [line(p1s[i][0],p1s[i][1],p1s[i][2],dis,azn,dip,'screen',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,pid[i])]
            #add to model domain
            self.wells = wells
        elif style == 'AGS': #closed loop geothermal
            n = len(self.wells)
            wells = self.wells
            segle = self.rock.ResDepth/20.0 #segment length, m
            separ = self.rock.w_spacing #well surface separation, m
            branc = int(self.rock.w_count/2) #single-sided number of branches (mirrored about center well) 
            #input overrides for economics
            self.rock.w_proportion = 1.0
            #general geometry
            azn = self.rock.w_azimuth
            dip = self.rock.w_dip
            vAxi = np.asarray([math.sin(azn)*math.cos(-dip),math.cos(azn)*math.cos(-dip),math.sin(-dip)])
            vert_z0 = np.sin(dip)*self.rock.w_length
            orig_d0 = self.rock.ResDepth-np.sin(dip)*0.5*self.rock.w_length
            hori_r0 = self.rock.w_length*np.cos(dip)
            offs_z0 = self.rock.w_spacing/np.cos(dip)
            #place injector and producer wells as the first two segments
            inje_c0 = np.asarray([-0.5*hori_r0*np.sin(azn),-0.5*hori_r0*np.cos(azn),orig_d0])
            prod_c0 = np.asarray([-(0.5*hori_r0-separ)*np.sin(azn),-(0.5*hori_r0-separ)*np.cos(azn),orig_d0])
            wells += [line(inje_c0[0],inje_c0[1],inje_c0[2],segle,0.0*deg,90.0*deg,'injector',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
            wells += [line(prod_c0[0],prod_c0[1],prod_c0[2],segle,0.0*deg,90.0*deg,'producer',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
            #segment information: i = injector; v = vertical; l = length; n = count
            ivl = self.rock.ResDepth - vert_z0 - segle
            ivn = int(ivl/segle) + 1
            ihl = self.rock.w_length
            ihn = int(ihl/segle) + 1
            pvl = self.rock.ResDepth - vert_z0 - offs_z0 + separ*np.tan(dip) - segle
            pvn = int(pvl/segle) + 1
            phl = self.rock.w_length - separ/np.cos(dip)
            phn = int(phl/segle) + 1
            #place injector segments
            oi = inje_c0 + np.asarray([0,0,-segle])
            ol = np.zeros(3,dtype=float)
            for b in range(0,ivn):
                wells += [line(oi[0],oi[1],oi[2], ivl/ivn, 0.0*deg,90.0*deg, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
                oi = oi + np.asarray([0.0,0.0,-ivl/ivn])
            ol = ol + oi
            for b in range(0,ihn):
                wells += [line(oi[0],oi[1],oi[2], ihl/ihn, azn, dip, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
                oi = oi + vAxi*(ihl/ihn)
            #place producer segments
            op = prod_c0 + np.asarray([0,0,-segle])
            ou = np.zeros(3,dtype=float)
            for b in range(0,pvn):
                wells += [line(op[0],op[1],op[2], pvl/pvn, 0.0*deg,90.0*deg, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
                op = op + np.asarray([0.0,0.0,-pvl/pvn])
            ou = ou + op
            for b in range(0,phn):
                wells += [line(op[0],op[1],op[2], phl/phn, azn, dip, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
                op = op + vAxi*(phl/phn)
            #place bottom hole connector
            wells += [line(oi[0],oi[1],oi[2], offs_z0, 0.0*deg, -90.0*deg, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
            #transverse vector
            vNor = np.cross(vAxi,np.asarray([0,0,1]))
            vNor = vNor/np.linalg.norm(vNor)
            #add branches
            for v in range(0,branc):
                #add transverse elements
                wells += [line(ou[0],ou[1],ou[2], separ*(v+1), azn+90.0*deg,0.0*deg,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
                wells += [line(ou[0],ou[1],ou[2], separ*(v+1), azn-90.0*deg,0.0*deg,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
                wells += [line(ol[0],ol[1],ol[2], separ*(v+1), azn+90.0*deg,0.0*deg,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
                wells += [line(ol[0],ol[1],ol[2], separ*(v+1), azn-90.0*deg,0.0*deg,'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
                #place injector segments
                oi = np.copy(ol) + vNor*separ*(v+1)
                for b in range(0,ihn):
                    wells += [line(oi[0],oi[1],oi[2], ihl/ihn, azn, dip, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
                    oi = oi + vAxi*(ihl/ihn)
                #place producer segments
                op = np.copy(ou) + vNor*separ*(v+1)
                for b in range(0,phn):
                    wells += [line(op[0],op[1],op[2], phl/phn, azn, dip, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
                    op = op + vAxi*(phl/phn)
                #place bottom connectors
                wells += [line(oi[0],oi[1],oi[2], offs_z0, 0.0*deg, -90.0*deg, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
                #place injector segments
                oi = np.copy(ol) - vNor*separ*(v+1)
                for b in range(0,ihn):
                    wells += [line(oi[0],oi[1],oi[2], ihl/ihn, azn, dip, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
                    oi = oi + vAxi*(ihl/ihn)
                #place producer segments
                op = np.copy(ou) - vNor*separ*(v+1)
                for b in range(0,phn):
                    wells += [line(op[0],op[1],op[2], phl/phn, azn, dip, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
                    op = op + vAxi*(phl/phn)  
                #place bottom connectors
                wells += [line(oi[0],oi[1],oi[2], offs_z0, 0.0*deg, -90.0*deg, 'pipe',self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+n)]
            #add to model domain
            self.wells = wells
        elif style == 'XGS': #patterned geothermal systems
            n = len(self.wells)
            wells = self.wells
            #center
            i0 = np.asarray([0.0, 0.0, 0.0])
            #ref axes (parallel injector and 90 horizontal to the right)
            azn = self.rock.w_azimuth
            dip = self.rock.w_dip
            vInj = np.asarray([math.sin(azn)*math.cos(-dip),
                     math.cos(azn)*math.cos(-dip),math.sin(-dip)])
            vRht = np.asarray([math.sin(azn+pi/2),math.cos(azn+pi/2),0.0])
            vNor = np.asarray([math.sin(azn)*math.sin(-dip),
                     math.cos(azn)*math.sin(-dip),math.cos(-dip)])
            #offset center
            p0 = i0 + vRht*self.rock.w_spacing
            #toe in production well
            toe = self.rock.w_toe
            vPro = sg.rotatePoints([vInj],vNor,toe)[0]
            #skew in production well
            skew = self.rock.w_skew
            vPro = sg.rotatePoints([vPro],vRht,skew)[0]
            #phase of production wells
            p0s = []
            vPros = []
            phase = self.rock.w_phase
            num = self.rock.w_count
            for i in range(0,num):
                p0s += [sg.rotatePoints([p0],vInj,(i*2.0*pi/num + phase))[0]]
                vPros += [sg.rotatePoints([vPro],vInj,(i*2.0*pi/num + phase))[0]]
            #lengths of wells
            length = self.rock.w_length
            proportion = self.rock.w_proportion
            iLen = length*proportion
            pLen = length
            #injection well segments
            i1s = []
            i2s = []
            seg = self.rock.w_intervals
            leg = iLen/(seg+seg-1)
            i1s += [i0 - 0.5*vInj*iLen]
            i2s += [i1s[0] + leg*vInj]
            for i in range(1,seg):
                i1s += [i2s[i-1] + leg*vInj]
                i2s += [i1s[i] + leg*vInj]
            #length of producers
            p1s = []
            p2s = []
            for j in range(0,num):
                p1s += [p0s[j] - 0.5*vPros[j]*pLen]
                p2s += [p0s[j] + 0.5*vPros[j]*pLen]
            #depth to model origin
            oDepth = self.rock.ResDepth + i2s[-1][2]            
            #place injection wells
            wells = []
            azn, dip, dis = azn_dip(i1s[0],i2s[0])
            vAxi = np.asarray([math.sin(azn)*math.cos(-dip),math.cos(azn)*math.cos(-dip),math.sin(-dip)])
            for i in range(0,seg):
                wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],0.1*leg,azn,dip,'injector',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,i+n)]
            #place production well surface segments
            for j in range(0,num):
                extra = 0.05*self.rock.ResDepth
                wells += [line(p1s[j][0],p1s[j][1],oDepth+extra,extra,0.0,90*deg,'producer',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+i+j+n)]
            #place injection well lower segments
            wells += [line(i1s[0][0],i1s[0][1],oDepth,oDepth-i1s[0][2],0.0,90.0*deg,'pipe',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
            #place perforation clusters
            for i in range(0,seg):
                i1s[i] = i1s[i] + vAxi*0.1*leg
                wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],0.9*leg,azn,dip,'perfcluster',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,i+n)]
            #place production well lower segments
            for j in range(0,num):
                #vertical
                wells += [line(p1s[j][0],p1s[j][1],oDepth,oDepth-p1s[j][2],0.0,90.0*deg,'pipe',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+i+j+n)]
                #inclined
                azn, dip, dis = azn_dip(p1s[j],p2s[j])
                vAxi = np.asarray([math.sin(azn)*math.cos(-dip),math.cos(azn)*math.cos(-dip),math.sin(-dip)])
                wells += [line(p1s[j][0],p1s[j][1],p1s[j][2],pLen,azn,dip,'screen',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+i+j+n)]
            #add to model domain
            self.wells = wells
        else: #CGS well placement (default)
            n = len(self.wells)
            wells = self.wells
            #center
            i0 = np.asarray([0.0, 0.0, 0.0])
            #ref axes (parallel injector and 90 horizontal to the right)
            azn = self.rock.w_azimuth
            dip = self.rock.w_dip
            vInj = np.asarray([math.sin(azn)*math.cos(-dip),
                     math.cos(azn)*math.cos(-dip),math.sin(-dip)])
            vRht = np.asarray([math.sin(azn+pi/2),math.cos(azn+pi/2),0.0])
            vNor = np.asarray([math.sin(azn)*math.sin(-dip),
                     math.cos(azn)*math.sin(-dip),math.cos(-dip)])
            #offset center
            p0 = i0 + vRht*self.rock.w_spacing
            #toe in production well
            toe = self.rock.w_toe
            vPro = sg.rotatePoints([vInj],vNor,toe)[0]
            #skew in production well
            skew = self.rock.w_skew
            vPro = sg.rotatePoints([vPro],vRht,skew)[0]
            #phase of production wells
            p0s = []
            vPros = []
            phase = self.rock.w_phase
            num = self.rock.w_count
            for i in range(0,num):
                p0s += [sg.rotatePoints([p0],vInj,(i*2.0*pi/num + phase))[0]]
                vPros += [sg.rotatePoints([vPro],vInj,(i*2.0*pi/num + phase))[0]]
            #lengths of wells
            length = self.rock.w_length
            proportion = self.rock.w_proportion
            iLen = length*proportion
            pLen = length
            #injection well segments
            i1s = []
            i2s = []
            seg = self.rock.w_intervals
            leg = iLen/(seg+seg-1)
            i1s += [i0 - 0.5*vInj*iLen]
            i2s += [i1s[0] + leg*vInj]
            for i in range(1,seg):
                i1s += [i2s[i-1] + leg*vInj]
                i2s += [i1s[i] + leg*vInj]
            #length of producers
            p1s = []
            p2s = []
            for j in range(0,num):
                p1s += [p0s[j] - 0.5*vPros[j]*pLen]
                p2s += [p0s[j] + 0.5*vPros[j]*pLen]
            #depth to model origin
            oDepth = self.rock.ResDepth + i2s[-1][2]            
            #place injection wells
            wells = []
            azn, dip, dis = azn_dip(i1s[0],i2s[0])
            vAxi = np.asarray([math.sin(azn)*math.cos(-dip),math.cos(azn)*math.cos(-dip),math.sin(-dip)])
            for i in range(0,seg):
                wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],0.1*leg,azn,dip,'injector',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,i+n)]
            #place production well surface segments
            for j in range(0,num):
                extra = 0.05*self.rock.ResDepth
                wells += [line(p1s[j][0],p1s[j][1],oDepth+extra,extra,0.0,90*deg,'producer',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+i+j+n)]
            #place injection well lower segments
            wells += [line(i1s[0][0],i1s[0][1],oDepth,oDepth-i1s[0][2],0.0,90.0*deg,'pipe',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,0+n)]
            #place perforation clusters
            for i in range(0,seg):
                i1s[i] = i1s[i] + vAxi*0.1*leg
                wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],0.9*leg,azn,dip,'perfcluster',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,i+n)]
            #place production well lower segments
            for j in range(0,num):
                #vertical
                wells += [line(p1s[j][0],p1s[j][1],oDepth,oDepth-p1s[j][2],0.0,90.0*deg,'pipe',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+i+j+n)]
                #inclined
                azn, dip, dis = azn_dip(p1s[j],p2s[j])
                vAxi = np.asarray([math.sin(azn)*math.cos(-dip),math.cos(azn)*math.cos(-dip),math.sin(-dip)])
                wells += [line(p1s[j][0],p1s[j][1],p1s[j][2],pLen,azn,dip,'screen',
                               self.rock.ra,self.rock.rb,self.rock.rc,self.rock.rgh,1+i+j+n)]
            #add to model domain
            self.wells = wells

    # ************************************************************************
    # stimulation - add frac
    # ************************************************************************
    def add_frac(self,typ = 'propped',
                      c0 = np.asarray([0.0,0.0,0.0]),
                      dia = [1266.0,107.0],
                      azn = [65.0*deg,4.0*deg],
                      dip = [57.5*deg,3.75*deg]):
        print( '   + placing new frac')
        fleg = np.random.normal(0,1)*dia[1] + dia[0]
        fazn = np.random.normal(0,1)*azn[1] + azn[0]
        fdip = np.random.normal(0,1)*dip[1] + dip[0]
        self.hydfs += [surf(c0[0],c0[1],c0[2],fleg,fazn,fdip,typ,self.rock)]
    
    # ************************************************************************
    # find intersections
    # ************************************************************************
    def gen_pipes(self,plot=True):
#        print( '*** intersections module ***')
        #This code segment:
        #   1. finds intersections & break lines into segments for the flow
        #   2. builds the input deck for the flow model                              
        #   3. be compatible with heat transfer model
        #   4. not produce duplicate segments
        #   5. be as efficient as I can muster
        #   6. reject dead ends if possible
        #   7. add infinite fracture elements to nodes on boundary
#        #re-initialize the geometry - building list of faces
#        self.re_init()
        #calculate intersections for each well
        found = False
        for w in range(0,len(self.wells)):
            found += self.x_well_all_faces(sourceID=w)
        if found > 0:
            #chain intersections without repeating same comparitors
            iters = 0
            maxit = 20
            lock = 0
            hold = 0
            while 1:
                #loop breaker
                iters += 1
                if iters > maxit:
                    print( '-> intersection search stopped after 20 iterations')
                    break
                #focus on chain connecting back to wells
                track = np.asarray(self.trakr)
                hold = len(self.trakr)
                conn = True
                #remove duplicate sources
                s_s = track[lock:,1]
                s_s = np.unique(s_s,axis=0)
                for s in s_s:
                    if s >= 0:
                        #search through all faces
                        for t in range(0,len(self.faces)):
                            #check if already searched
                            ck_1 = np.asarray([np.isin(track[:,0],t),np.isin(track[:,1],s)])
                            ck_2 = ck_1[0,:]*ck_1[1,:]
                            ck_3 = np.sum(ck_2)
                            #evaluate if value is valid
                            if (s != t) and (t >= 0) and (ck_3 == 0) and not(s in track[:,0]): #this pair is fresh, check for intersections
                                self.x_frac_face(plot=plot,sourceID=s,targetID=t)
                                conn = False
                        #update tracker
                        self.trakr += [[s,-1]]
                #lockout repeat searches
                lock = hold
                #break if no connections in this search
                if conn:
                    break
            #report number of iterations used to find intersections
            print( '   x all intersections found using %i iters of %i allowable' %(iters,maxit))
        else:
            print( '   # wells do not intersect any fractures')
            #add phantom paths
            # count = 0
            # for w in range(0,len(self.wells)):
            #     if (int(self.wells[w].typ) in [typ('injector'),typ('producer')]):
            #         count += 1
            #         if count == 1:
            #             c1 = self.wells[w].c0
            #         else:
            #             c0 = self.wells[w].c0
            #             self.add_flowpath(c0,c1,1.0,1e-10,typ('shunt'),w,1e-10,1e-10,80.0,self.wells[w].pID)
        #shunts for solver stability
        for w in range(0,len(self.wells)):
            if (int(self.wells[w].typ) in [typ('injector'),typ('producer')]):
                c1 = self.nodes.all[0]
                c0 = self.wells[w].c0
                self.add_flowpath(c0,c1,1.0,1e-10,typ('shunt'),w,1e-10,1e-10,80.0,self.wells[w].pID)
        self.pipes.coords(self.nodes)

    # ************************************************************************
    #energy generation - single flash steam rankine cycle
    # ************************************************************************
    def get_power(self,detail=False):
        #initialization
        self.Qout = []
        self.Pout = []
        self.Tout = []
        self.Hout = []
        Pump_Power = []
        Net_Power = []
        Bulk_Power = []
        Therm_Power = []
        Tpro = []
        hpro = []
        #injection fluid at surface (5)
        P5 = np.max(self.p_p) - self.rock.BH_P# + self.rock.p_whp #Pa; i_p
        T5 = self.rock.TinjK #K
        h5 = water.h_from_PT(T=T5,P=P5) #kJ/kg
        v5 = water.v_from_PT(T=T5,P=P5) #m3/kg
        #massflows
        mt = -self.p_mm
        mi = self.i_mm
        Vt = (v5*mt)*(1000*60) # L/min
        Vi = (v5*mi)*(1000*60) # L/min
        #for each moment in time
        for t in range(0,len(self.p_hm)-1):
            # Undisturbed Reservoir (r)
            Tr = self.rock.BH_T #K
            Pr = self.rock.BH_P #Pa
            hr = water.h_from_PT(T=Tr,P=Pr)
            vr = water.v_from_PT(T=Tr,P=Pr)
            # Surface Production Well (2)
            P2 = self.rock.p_whp #Pa
            h2 = self.p_hm[t] #kJ/kg
            T2 = water.T_from_Ph(P=P2,h=h2)
            v2 = water.v_from_PT(P=P2,T=T2)
            Tpro += [T2]
            hpro += [h2]
            # Bulk power
            effy = np.max([0.0, 0.078795*np.log(h2) - 0.45651]) #(Zarrouk and Moon, 2014: Geothermics) --- f(electricity out, thermal in)
            Bulk = effy*h2*mt
            Pump = -1e-3*mi*v5*P5/self.rock.GenEfficiency #reinjection pumping only
            # Net power
            Net = Bulk + Pump
            # Record results
            Pump_Power += [Pump]
            Net_Power += [Net]
            Bulk_Power += [Bulk] #[Flash + Binary]
            Therm_Power += [mt*(h2 - h5)]
        #record results
        self.Qout = np.asarray(Pump_Power)
        self.Pout = np.asarray(Net_Power)
        self.Tout = np.asarray(Bulk_Power)
        self.Hout = np.asarray(Therm_Power)
        self.payoff.Injection_Rate_m3ps = v5*mi
        self.payoff.Injection_Rate_kgps = mi
        self.payoff.Production_Pressure_MPa = P2/MPa
        self.payoff.Production_Rate_kgps = mt
        self.payoff.Production_Temperature_K = np.sum(Tpro)/len(Tpro)
        self.payoff.Production_Enthalpy_kJpkg = np.sum(hpro)/len(hpro)
        
        #print( details)
        if detail:
            print( '\n*** Electrical Generation Thermal State Values ***')
            print( ("Inject (5): T= %.2f; P= %.2f; h= %.2f, v= %.6f" %(T5,P5,h5,v5)))
            print( ("Reserv (r,1): T= %.2f; P= %.2f; h= %.2f, v= %.6f" %(Tr,Pr,hr,vr)))
            print( ("Produc (2): T= %.2f; P= %.2f; h= %.2f, v= %.6f" %(T2,P2,h2,v2)))
            print( '*** Power Output Estimation ***')
            print( "Injection Rate = %.2f kg/s = %.2f L/min" %(mi, Vi))
            print( "Production Rate = %.2f kg/s = %.2f L/min" %(mt, Vt))
            print( "Bulk Power at %.2f kW" %(Bulk))
            print( "Pumping Power at %.2f kW" %(Pump))
            print( "Net Power at %.2f kW" %(Net))
            
    # ************************************************************************
    # flow network model
    # ************************************************************************
    def get_flow(self,p_bound=0.0*MPa,q_well=[],p_well=[],reinit=True,useprior=False,Qnom=1.0,auto_farfield=False):
        #reinitialize if the mesh has changed
        if reinit:
            #clear data from prior run
            self.re_init()
            #generate pipes
            self.gen_pipes()
        
        #set boundary conditions (m3/s) (Pa); 10 kg/s ~ 0.01 m3/s for water
        self.set_bcs(p_bound=p_bound,q_well=q_well,p_well=p_well,auto_farfield=auto_farfield)
        
        #get fluid properties
        mu = self.rock.Poremu #cP
        rho = self.rock.PoreRho #kg/m3
        
        #flow solver working variables
        N = self.nodes.num
        Np = self.pipes.num
        H = self.H
        Q = self.Q
            
        #initial guess for nodal pressure head
        if useprior and not(reinit):
            h = self.nodes.p/(rho*g)
        else:
            h = (1.0*MPa/(rho*g))*np.random.rand(N) + p_bound/(rho*g)
        q = np.zeros(N)
        
        #install boundary condition
        for i in range(0,len(H)):
            h[H[i][0]] = H[i][1]
        for i in range(0,len(Q)):
            q[Q[i][0]] = Q[i][1]
        
        #hydraulic resistance equations
        K = np.zeros(Np)
        n = np.zeros(Np)
        
        #stabilizing limiters
        hup = (self.rock.s1+10.0*MPa)/(rho*g)
        hlo = -101.4/(rho*g)
        
        #convergence
        goal = 1e-4 #1.0e-8 #0.0001
        
        #dimension limiters for flow solver stability
        self.pipes.Dh_limit(Qnom,goal*rho*g,rho,g,mu) #,self.rock.kf)
        
        #hydraulic resistance terms
        for i in range(0,Np):
            #working variables
            u = self.pipes.fID[i]
            
            #pipes and wells #modified 10/6/23 - go jags!
            if (int(self.pipes.typ[i]) in [typ('injector'),typ('producer'),typ('pipe'),typ('perfcluster'),typ('screen')]):
                self.pipes.Dh[i] = np.max([np.min([self.pipes.Dh_max[i],2.0*self.wells[u].ra]),1e-12]) # added 1e-12 lower limit to prevent div0
                Lscaled = self.pipes.L[i]*mu/(0.9*cP)
                K[i] = 10.7*Lscaled/(self.pipes.frict[i]**1.852*self.pipes.Dh[i]**4.87) #viscosity scaled Hazen-Williams (is this valid?)
                n[i] = 1.852
            #fractures and planes
            elif (int(self.pipes.typ[i]) in [typ('boundary'),typ('fracture'),typ('propped'),typ('choke')]):
                self.pipes.Dh[i] = np.max([np.min([self.pipes.Dh_max[i],self.faces[u].bh]),1e-12])
                K[i] = (12.0*mu*self.pipes.L[i])/(rho*g*self.pipes.W[i]*self.pipes.Dh[i]**3.0)
                n[i] = 1.0      
            #casing perforations
            elif (int(self.pipes.typ[i]) in [typ('perf')]):
                self.pipes.Dh[i] = np.max([np.min([self.pipes.Dh_max[i],self.rock.perf_dia]),1e-12])
                K[i] = 4.0*mu/(0.9*cP*self.rock.perf_per_cluster**2.0*self.pipes.Dh[i]**4.0*0.825**2.0*g)
                n[i] = 2.0
            #porous media
            elif (int(self.pipes.typ[i]) in [typ('darcy')]):
                self.pipes.Dh[i] = np.max([np.min([self.pipes.Dh_max[i],self.faces[u].bd]),1e-12])
                # K[i] = mu*self.pipes.L[i]/(rho*g*self.pipes.Dh[i]*self.pipes.W[i]*self.rock.kf)
                K[i] = mu*self.pipes.L[i]/(rho*g*self.pipes.Dh[i]*self.pipes.W[i]*self.pipes.frict[i])
                n[i] = 1.0
            elif (int(self.pipes.typ[i]) in [typ('shunt')]):
                K[i] = 1e12
                n[i] = 1.0
            #type not handled
            else:
                print( 'error: undefined type of conduit')
                exit
        #record info
        self.pipes.K = K
        self.pipes.n = n
    
        #iterative Newton-Rhapson solution to solve flow
        iters = 0
        max_iters = 50
        z = h
        while 1:
            #loop breaker
            iters += 1
            if iters > max_iters:
                print( '   > Flow solver converged to <%.2e m head after %i iterations' %(np.max(np.abs(z)),iters-1))
                break
            elif (np.max(np.abs(z)) < goal) and (iters > 3): #np.max(np.abs(z/(h+z))) < goal:
                print( '   > Flow solver converged to <%.2e m head using %i iterations' %(goal,iters-1))
                break
            
            #re-initialize working variables
            F = np.zeros(N)
            F += q
            D = np.zeros((N,N))
            
            #build matrix equations
            for i in range(0,Np):
                #working variables
                n0 = self.pipes.n0[i]
                n1 = self.pipes.n1[i]
                
                #Node flow equations
                R = np.sign(h[n0]-h[n1])*(abs((h[n0]-h[n1]))/K[i])**(1.0/n[i])
                F[n0] += R
                F[n1] += -R
                
                #Jacobian (first derivative of the inflow-outflow equations)
                if abs((h[n0]-h[n1])) == 0:
                    J = 1.0
                else:
                    J = ((1.0/n[i])/K[i])*(abs((h[n0]-h[n1]))/K[i])**(1.0/n[i]-1.0)
                D[n0,n0] += J
                D[n1,n1] += J
                D[n0,n1] += -J
                D[n1,n0] += -J
            
            #remove defined boundary values
            for i in range(0,len(H)):
                F = np.delete(F,H[-1-i][0],0)
                D = np.delete(D,H[-1-i][0],0)
                D = np.delete(D,H[-1-i][0],1)
            
            #solve matrix equations
            z = solve(D[:,:],F[:])
            
            #update pressures
            for i in range(0,len(H)):
                z = np.insert(z,H[i][0],0,axis=0)
            h = h - 0.9*z #0.7 scaling factor seems to give faster convergence by reducing solver overshoot
            
            #apply physical limits
            h[h > hup] = hup
            h[h < hlo] = hlo
        
        #flow rates
        q = np.zeros(Np)
        for i in range(0,Np):
            #working variables
            n0 = self.pipes.n0[i]
            n1 = self.pipes.n1[i]
            #flow rate
            q[i] = np.sign(h[n0]-h[n1])*(abs(h[n0]-h[n1])/K[i])**(1.0/n[i])
            
        #record results in class
        self.nodes.p = h*rho*g
        self.q = q
        self.pipes.q = np.abs(q)

        #collect well rates and pressures
        p_q = np.zeros(len(self.wells),dtype=float)
        p_p = np.zeros(len(self.wells),dtype=float)
        
        for w in range(0,len(self.wells)):
            #coordinates
            source = self.wells[w].c0
            #find index of duplicate
            ck, i = self.nodes.add(source)
            #record pressure
            p_p[w] = self.nodes.p[i]
            #record flow rate
            i_pipe = np.where(np.asarray(self.pipes.n0) == i)[0][0]
            p_q[w] = self.q[i_pipe]
                
        #collect boundary rates and pressures
        b_nodes = [0]
        b_pipes = []
        b_q = []
        b_p = []
        w = b_nodes[0]
        b_pipes = np.where(np.asarray(self.pipes.n1) == w)[0]
        b_p += [self.nodes.p[w]]
        if len(b_pipes) > 0:
            for w in b_pipes:
                b_q += [-self.q[w]]
        b_p = list(np.ones(len(b_q),dtype=float)*b_p[0])
                
        #store key well and boundary flow data
        self.p_p = p_p
        self.p_q = p_q
        self.b_p = b_p
        self.b_q = b_q
    
    # ************************************************************************
    # heat transfer model 
    # - mod 1-28-2021: correct errors
    # - mod 9-13-2022: reduce overshoot in initial timesteps
    # - mod 3-27-2024: replace IAPWS with properties lookup-table
    # ************************************************************************
    def get_heat(self,plot=True,
                 detail=False,
                 lapse=False):
        print( '*** heat flow module ***')
        #****** default parameters ******
        t_n = self.rock.TimeSteps
        t_f = self.rock.LifeSpan
        H = self.rock.H_ConvCoef
        dT0 = self.rock.dT0
        
        #****** boundary parameters ******
        # Surface Injection Well (5)
        P5 = np.max(self.p_p) - self.rock.BH_P# + self.rock.p_whp #Pa; i_p
        self.payoff.Injection_Pressure_MPa = P5/MPa #MPa
        T5 = self.rock.TinjK #K
        h5 = water.h_from_PT(P=P5,T=T5) #kJ/kg
        v5 = water.v_from_PT(P=P5,T=T5) #m3/kg
        self.v5 = v5 #m3/kg
        print( ("Inject (5): T= %.2f; P= %.2f; h= %.2f, v= %.6f" %(T5,P5,h5,v5)))
        self.payoff.Injection_Temperature_K = T5 #K
        self.payoff.Injection_Enthalpy_kJpkg = h5 #kJ/kg
        
        # Undisturbed Reservoir (r)
        Tr = self.rock.BH_T #K
        Pr = self.rock.BH_P #Pa
        hr = water.h_from_PT(P=Pr,T=Tr)
        vr = water.v_from_PT(P=Pr,T=Tr)
        print( ("Reserv (r): T= %.2f; P= %.2f; h= %.2f, v= %.6f" %(Tr,Pr,hr,vr)))
        
        #****** enthalpy function linearization *******
        water.set_Pref(P=Pr)
        
        #****** nodal rock temperature initialization *******
        self.nodes.Tr = (self.rock.ResDepth - self.nodes.all[:,2]) * self.rock.ResGradient + self.rock.AmbTempK
        self.nodes.Tr[0] = Tr

        #****** explicit solver setup ******
        #working variables
        ts = np.linspace(0.0,t_f,t_n+1)
        dt = ts[1]-ts[0]
        #set boundary conditions #K
        self.therm_bcs(T_bound=Tr,T_inlet=[T5])
        #set solver parameters
        N = self.nodes.num
        Np = self.pipes.num
        Tb = self.Tb
        #initial guess for temperatures
        Tn = np.ones(N)*self.nodes.Tr #K
        
        #install boundary condition
        for i in range(0,len(Tb)):
            Tn[Tb[i][0]] = Tb[i][1]
            
        #more working variables
        Lp = np.asarray(self.pipes.L)
        ms = np.asarray(self.q)/v5
        hn = np.zeros(N)
        R0 = np.zeros(Np)
        ResSv = self.rock.ResSv
        ResKt = self.rock.ResKt
        CemKt = self.rock.CemKt
        
        #memory working variables
        ht = np.zeros((t_n+1,N),dtype=float)
        Tt = np.zeros((t_n+1,N),dtype=float)
        Rt = np.zeros((t_n+1,Np),dtype=float)
        Et = np.zeros((t_n+1,Np),dtype=float)
        Qt = np.zeros((t_n+1,Np),dtype=float)
        Er = np.zeros(Np,dtype=float)
        
        #well geometry
        Ris = self.rock.ra
        Ric = self.rock.rb
        Rir = self.rock.rc
        
        #functional form of energy vs thermal radius for specified borehole diameter (per unit dT and dL)
        Ror = np.linspace(1,2000,2000)
        ERor = 2*pi*ResSv*1.0*((1-np.log(Rir)/np.log(Rir/Ror))*(Ror**2-Rir**2)/2
            + (1/(2*np.log(Rir/Ror)))*(Ror**2*(np.log(Ror)-0.5)-Rir**2*(np.log(Rir)-0.5)))
        lnRor = np.log(Ror[:])
        lnERor = np.log(ERor[:])
        ERm, ERb, r_value, p_value, std_err = stats.linregress(lnERor,lnRor)
        
        # #key indecies for convergence criteria (i.e., boundary nodes)
        # key = []
        # for e in self.H:
        #     key += [e[0]]
        #     key += [e[0]+1]
        # for e in self.Q:
        #     key += [e[0]]
        #     key += [e[0]+1]
    
        #convergence criteria
        goal = 0.5
        goalE = 0.03

        #iterate over time
        for t in range(0,len(ts)-1):
            #calculate nodal temperatures by pipe flow and nodal mixing
            iters = 0
            err = np.ones(N)*goal*10
            E0_update_intervals = 5
            max_iters = 20*E0_update_intervals #gives 100 max iterations
            E0e = np.zeros(Np) + Et[0,:]
            dE0e = np.zeros(Np)
            
            #iterate for temperature stability
            while 1:
                #calculate nodal fluid enthalpy #@@@ requires consideration of sensible heating (0<X<1) if not supercritical
                for n in range(0,N):
                    hn[n] = water.h_from_T(T=Tn[n]) #kJ/kg
                # hn = hTP[0]*Tn**3.0 + hTP[1]*Tn**2.0 + hTP[2]*Tn**1.0 + hTP[3] #!!! old method has unresolved issues with overflow and invalid values
                #temperature convergence
                kerr = np.max(np.abs(err))
                #energy convergence
                eerr = np.max(np.abs(dE0e)/(np.abs(E0e)+1.0e-9))
                #loop breaker
                iters += 1
                if iters > max_iters:
                    print( 'Heat solver converged to %e K after %i iterations' %(kerr,iters-1))
                    break
                elif (kerr < goal) and (eerr < goalE):
                    print( 'Heat solver converged to %e K after %i iterations' %(kerr,iters-1))
                    break
                #follow the flow to estimate heating/cooling of fluid per pipe
                Ap = np.zeros(Np,dtype=float)
                Bp = np.zeros(Np,dtype=float)
                Cp = np.zeros(Np,dtype=float)
                Ep = np.zeros(Np,dtype=float)
                hm = np.zeros(N,dtype=float)
                mi = np.zeros(N,dtype=float)
                for p in range(0,Np):
                    #working variables
                    n0 = self.pipes.n0[p]
                    n1 = self.pipes.n1[p]
                    Tr0 = self.nodes.Tr[n0] 
                    Tr1 = self.nodes.Tr[n1] 
                    #get overshoot limit for timestep for each pipe
                    if (iters-1)%E0_update_intervals == 0:
                        #dT = 0.5*(Tn[n1] + Tn[n0]) - Tr #!!! solution before gradient
                        dT = 0.5*(Tn[n1] + Tn[n0]) - 0.5*(Tr0 + Tr1)
                        dT = np.max([np.abs(dT0),np.abs(dT)])
                        a = 1.0/(self.rock.ResSv*np.abs(dT)*self.rock.ResKt*10**-3)
                        b = 1.0/self.rock.H_ConvCoef
                        c = -2.0*np.abs(dT)*dt
                        dE0 = (-b + (b**2.0 - 4.0*a*c)**0.5)/(2.0*a)
                        #get extracted energy - Et, thermal radius - R0, heat flow rate - Qt
                        if (int(self.pipes.typ[p]) in [typ('injector'),typ('producer'),typ('pipe'),typ('perfcluster'),typ('screen'),typ('perf')]):  #pipe, radial heat flow
                            #energy withdraw
                            E0p = dE0*2.0*pi*Rir*Lp[p] #kJ
                            Esp = Et[t,p]
                            Etp = np.max([E0p,Esp])
                            dE0e[p] = Etp - E0e[p]
                            E0e[p] = Etp
                            #initial rock thermal radius
                            R0[p] = np.exp(ERm*(np.log(np.abs(Etp/(Lp[p]*dT))))+ERb) + Rir #m
                            #initial rock energy transfer rates
                            Qt[t,p] = Lp[p]/(1.0/(2.0*pi*Ris*H) + np.log(R0[p]/Rir)/(2.0*pi*ResKt*10**-3) + np.log(Rir/Ric)/(2.0*pi*CemKt*10**-3)) #kJ/K-s
                        elif (int(self.pipes.typ[p]) in [typ('boundary'),typ('fracture'),typ('propped'),typ('darcy'),typ('choke')]): #fracture, plate heat flow
                            #rock energy withdraw
                            E0p = dE0*self.pipes.W[p]*Lp[p] #kJ
                            Esp = Et[t,p]
                            Etp = np.max([E0p,Esp])
                            dE0e[p] = Etp - E0e[p]
                            E0e[p] = Etp
                            #rock thermal radius
                            R0[p] = Etp/(ResSv*self.pipes.W[p]*Lp[p]*dT) #m
                            #rock energy transfer rates
                            Qt[t,p] = (2.0*self.pipes.W[p]*Lp[p])/(1.0/(H) + R0[p]/(ResKt*10**-3)) #kJ/K-s
                        elif (int(self.pipes.typ[p]) in [typ('shunt')]):
                            #null
                            Qt[t,p] = 0.0
                            R0[p] = 1e-10
                            E0e[p] = 0.0
                            dE0e[p] = 0.0
                        else:
                            print( 'error: segment type %s not identified' %(typ(int(self.pipes.typ[p]))))

                    #equilibrium enthalpy
                    heq0 = water.h_from_T(T=Tr0)
                    heq1 = water.h_from_T(T=Tr1)
                    
                    #working variables
                    Kp = 0.0
                    Qp = 0.0
                    
                    #positive flow
                    if ms[p] > 0:
                        #non-equilibrium conduction limited heating
                        Ap[p] = Qt[t,p]*(0.5*(Tr1 + Tr0) - 0.5*(Tn[n1] + Tn[n0]))
                        #equilibrium conduction limited heating
                        Bp[p] = Qt[t,p]*(Tr1 - 0.5*(Tr1 + Tn[n0]))
                        #flow limited cooling to equilibrium
                        Cp[p] = ms[p]*(heq1 - hn[n0])
    
                    #negative flow    
                    else:
                        #non-equilibrium conduction limited heating
                        Ap[p] = Qt[t,p]*(0.5*(Tr0 + Tr1) - 0.5*(Tn[n0] + Tn[n1]))
                        #equilibrium conduction limited heating
                        Bp[p] = Qt[t,p]*(Tr0 - 0.5*(Tr0 + Tn[n1]))
                        #flow limited cooling to equilibrium
                        Cp[p] = -ms[p]*(heq0 - hn[n1])
                        
                    #take maximum of conduction terms because this will drive conduction heat deliverability
                    Kp = np.asarray([Ap[p],Bp[p]])[np.argmax(np.abs([Ap[p],Bp[p]]))]
                    #if flow limits heat extraction from rock, it will go to equilibrium
                    Qp = Cp[p]
                    #get the limiting term
                    KorQ = np.argmin(np.abs([Kp,Qp]))
                    Ep[p] = np.asarray([Kp,Qp])[np.argmin(np.abs([Kp,Qp]))]
                    
                    #conduction limited
                    if KorQ == 0:
                        #positive flow
                        if ms[p] > 0:
                            mi[n1] += ms[p]
                            hm[n1] += Ep[p] + ms[p]*hn[n0]
                        #negative flow
                        else:
                            mi[n0] += -ms[p]
                            hm[n0] += Ep[p] + -ms[p]*hn[n1]
                    #flow limited
                    else:
                        #positive flow
                        if ms[p] > 0:
                            mi[n1] += ms[p]
                            hm[n1] += ms[p]*heq1
                        #negative flow
                        else:
                            mi[n0] += -ms[p]
                            hm[n0] += -ms[p]*heq0
                
                #calculate nodal temperatures
                z = []
                z = np.zeros(N,dtype=float) #K
                hu = np.zeros(N,dtype=float) #kJ/kg
                for n in range(0,N):
                    #mixed inflow enthalpy
                    if mi[n] > 0:
                        hu[n] = hm[n]/mi[n]
                    #no inflow use rock enthalpy
                    else:
                        hu[n] = water.h_from_T(T=self.nodes.Tr[n]) #!!!
                    #calculate temperature at new enthalpy
                    z[n] = water.T_from_h(h=hu[n])
                    
                
                #install boundary condition
                for j in range(0,len(Tb)):
                    z[Tb[j][0]] = Tb[j][1]
                    hu[Tb[j][0]] = water.h_from_T(T=z[Tb[j][0]]) #!!!
                    
                #calculate error
                err = [] #np.zeros(N,dtype=float)
                err = Tn - z
                
                #update Tn
                Tn = z
                
            #store ht
            ht[t] = hu
            #store temperatures
            Tt[t,:] = Tn
            #store thermal radii
            Rt[t,:] = R0
            
            #timelapse 3D
            self.nodes.T = Tn
            self.nodes.h = hn
            if lapse:
                self.build_vtk(fname='t%02d' %(t))
    
            #extracted energy during this time step
            Et[t+1] = Et[t] + np.abs(Ep)*dt
            #rock energy tracker (added or lost)
            Er += Ep*dt
            for i in range(0,Np):
                #working variables
                n0 = self.pipes.n0[i]
                n1 = self.pipes.n1[i]
                Tr0 = self.nodes.Tr[n0] 
                Tr1 = self.nodes.Tr[n1] 
                dT = 0.5*(Tn[n1] + Tn[n0]) - 0.5*(Tr0 + Tr1)
                dT = np.max([np.abs(dT0),np.abs(dT)])

                #thermal radius for next time step
                if (int(self.pipes.typ[i]) in [typ('injector'),typ('producer'),typ('pipe'),typ('perfcluster'),typ('screen'),typ('perf')]): #pipe, radial heat flow
                    if (dT > 0) and (Et[t+1,i] > 0):
                        R0[i] = np.exp(ERm*(np.log(np.abs(Et[t+1,i]/(Lp[i]*dT))))+ERb) + Rir #+2.0*Rir # m
                    #Et[0,i] = np.exp((np.log(R0[i])-ERb)/ERm)*Lp[i]*(Tr-0.5*(Tn[Y[i][1]]+Tn[Y[i][0]])) # kJ
                    Qt[t+1,i] = Lp[i]/(1.0/(2.0*pi*Ris*H) + np.log(R0[i]/Rir)/(2.0*pi*ResKt*10**-3) + np.log(Rir/Ric)/(2.0*pi*CemKt*10**-3)) # kJ/K-s # Note converted Kt in W/m-K to kW/m-K
                    
                elif (int(self.pipes.typ[i]) in [typ('boundary'),typ('fracture'),typ('propped'),typ('darcy'),typ('choke')]): #fracture, plate heat flow
                    if (dT > 0) and (Et[t+1,i] > 0):
                        R0[i] = Et[t+1,i]/(ResSv*self.pipes.W[i]*Lp[i]*dT)
                    #Et[0,i] = ResSv*R0[i]*Y[i][4]*Lp[i]*(Tr-0.5*(Tn[Y[i][1]]+Tn[Y[i][0]])) # kJ
                    Qt[t+1,i] = (2.0*self.pipes.W[i]*Lp[i])/(1.0/(H) + R0[i]/(ResKt*10**-3)) # kJ/K-s # Note converted Kt in W/m-K to kW/m-K

                # elif (int(self.pipes.typ[i]) in [typ('perf')]): #pipe, radial heat flow
                #     R0[i] = 1.0e-6 # m
                #     #Et[0,i] = np.exp((np.log(R0[i])-ERb)/ERm)*Lp[i]*(Tr-0.5*(Tn[Y[i][1]]+Tn[Y[i][0]])) # kJ
                #     Qt[t+1,i] = 0.0 # kJ/K-s
                elif (int(self.pipes.typ[p]) in [typ('shunt')]):
                    #null
                    Qt[t+1,p] = 0.0
                    R0[p] = 1e-10

                else:
                    print( 'error: segment type %s not identified' %(typ(int(self.pipes.typ[i]))))
        
        #store results
        self.Et = Et
        self.Qt = Qt
        self.nodes.T = Tn
        self.nodes.h = hn
        self.ms = ms
        self.Tt = Tt
        self.ht = ht
        self.R0 = R0
        self.Rt = Rt
        
        #updated storage approach
        self.pipes.R0 = R0
        
        #parameters of interest
        #**********************
        #collect well temp, rate, enthalpy
        w_h = []
        w_T = []
        w_m = []
        for w in range(0,len(self.wells)):
            #type
            if (int(self.wells[w].typ) in [typ('injector'),typ('producer')]):
                #coordinates
                source = self.wells[w].c0
                #find index of duplicate
                ck, i = self.nodes.add(source)
                #record temperature
                w_T += [Tt[:,i]]
                #record enthalpy
                w_h += [ht[:,i]]
                #record mass flow rate
                i_pipe = np.where(np.asarray(self.pipes.n0) == i)[0][0]
                w_m += [self.q[i_pipe]/v5]
        w_h = np.asarray(w_h)
        w_T = np.asarray(w_T)
        w_m = np.asarray(w_m)
            
        #identify injectors and producers
        iPro = np.where(w_m<=0.0)[0]
        iInj = np.where(w_m>0.0)[0]
        
        #boundary flows
        b_nodes = [0]
        b_pipes = []
        b_h = [] 
        b_T = []
        b_m = []
        w = b_nodes[0]
        b_pipes = np.where(np.asarray(self.pipes.n1) == w)[0]
        b_h += [ht[:,w]]
        b_T += [Tt[:,w]]
        if len(b_pipes) > 0:
            for w in b_pipes:
                b_m += [self.q[w]/v5]
        
        #total production energy
        p_ET = 0.0
        for i in iPro:
            i = int(i)
            p_ET += np.sum(w_m[i]*w_h[i]) #kJ/s
        p_ET = p_ET*dt/t_f #/len(p_h[w]) #kJ/s avg over all
        
        #mixed produced enthalpy and mass flow rate
        p_mm = []
        p_hm = np.zeros(len(ts))
        p_mm = 0.0
        for i in iPro:
            i = int(i)
            p_mm += w_m[i]
            p_hm += w_h[i]*w_m[i]
        if p_mm < 0: 
            p_hm = p_hm/p_mm
        else:
            p_hm = np.ones(len(ts))*h5
        for i in p_hm:
            if i <= 0.0:
                i = water.h_from_PT(P=self.rock.AmbPres,T=self.rock.AmbTempK)
        self.p_mm = p_mm #mixed produced mass flow rate
        self.p_hm = p_hm #mixed produced enthalpy
                
        #total injection energy
        i_ET = 0.0
        for i in iInj:
            i = int(i)
            i_ET += np.sum(w_m[i]*w_h[i]) #kJ/s
        i_ET = i_ET*dt/t_f #/len(i_h[w])

        #mixed injection mass flow rate
        i_mm = []
        i_mm = 0.0
        for i in iInj:
            i = int(i)
            i_mm += w_m[i]
        self.i_mm = i_mm #mixed produced mass flow rate
        
        #total boundary energy flow
        b_ET_out = 0.0
        b_ET_in = 0.0
        b_ET = 0.0
        b_ni = 0
        b_no = 0
        if len(b_pipes) > 0:
            for w in range(0,len(b_pipes)):
                #outflow
                if b_m[w] > 0:
                    b_ET_out += -np.sum(b_m[w]*b_h[0])
                    b_no += 1
                #inflow
                else:
                    b_ET_in += -np.sum(b_m[w]*b_h[0])
                    b_ni += 1
            if b_no > 0:
                b_ET_out = b_ET_out*dt/t_f #/b_no
            if b_ni > 0:
                b_ET_in = b_ET_in*dt/t_f #/b_ni
            b_ET = b_ET_out + b_ET_in
    
        #rock energy for thermal radius
        #E_tot = (np.sum(Et[-1,:]) - np.sum(Et[0,:]))/t_f #kJ/s
        
        #rock energy change
        E_roc = np.sum(Er)/t_f #kJ/s
        
        #net system energy (less error in model if closer to zero)
        E_net = E_roc + p_ET + i_ET + b_ET #kJ/s
        print( '\nNet System Energy (kJ/s):' )
        print( E_net)
        
        #save key values
        self.ts = ts
        self.w_h = w_h
        self.w_m = w_m
        self.b_h = b_h
        self.b_m = b_m
        self.w_T = w_T
        
        #calculate power output
        self.get_power(detail=detail)
        
        # #thermal energy extraction
        # dht = np.zeros(len(ts),dtype=float)
        # for i in range(0,len(w_m)):
        #     dht += -w_m[i]*w_h[i]
        # self.dhout = dht

        if plot: #plots
            fig = pylab.figure(figsize=(6.0, 8.5), dpi=96, facecolor='w', edgecolor='k',tight_layout=True) # Medium resolution
            font = {'family' : 'serif',   'size'   : 12}
            pylab.rc('font', serif='Arial')
            pylab.rc('font', **font)
            #temperature
            ax3 = fig.add_subplot(211)
            for i in iInj:
                i = int(i)
                ax3.plot(ts[:-1]/yr,w_T[i][:-1],linewidth=1.5)
            for i in iPro:
                i = int(i)
                ax3.plot(ts[:-1]/yr,w_T[i][:-1],linewidth=1.5)
            ax3.set_xlabel('Time (yr)')
            ax3.set_ylabel('Production Temperature (K)')
            #power
            ax2 = fig.add_subplot(212)
            ax2.plot(ts[:-1]/yr,self.Tout[:],linewidth=1.0,color='blue')
            ax2.plot(ts[:-1]/yr,self.Qout[:],linewidth=1.0,color='red')
            ax2.plot(ts[:-1]/yr,self.Pout[:],linewidth=1.5,color='green')
            ax2.set_xlabel('Time (yr)')
            ax2.set_ylabel('Gross-Blue Pump-Cyan Net-Green (kWe)')

    #calculate flow with sequentially-coupled stress-sensitive mechanics
    def solve_hydromech(self,i_key,p_key,tip,reinit=True,fix=False,
                        stim=False,plim=True,tryflowbound=False,visuals=False):
        #call counter
        self.iters += 1
        if len(self.Qi) == 0:
            self.Qi = np.zeros(len(i_key))
        
        #get max pressure on each fracture from all the nodes associated with that fracture
        face_pmax = np.ones(len(self.faces),dtype=float)*self.rock.BH_P
        for i in range(0,self.pipes.num):
            if (int(self.pipes.typ[i]) in [typ('boundary'),typ('fracture'),typ('propped'),typ('choke'),typ('perf')]):
                face_pmax[self.pipes.fID[i]] = np.max([face_pmax[self.pipes.fID[i]], 
                                                       self.nodes.p[self.pipes.n0[i]], 
                                                       self.nodes.p[self.pipes.n1[i]]])
        #update fracture properties using the maximum nodal pressure (note: using center pressure doesn't work)
        Vii = 0.0
        stim_num = 0
        for i in range(0,len(self.faces)):
            #record maximum and center node pressures
            self.faces[i].Pmax = face_pmax[i]
            self.faces[i].Pcen = self.nodes.p[self.faces[i].ci]
            #compute fracture properties
            stim_num += self.hydromech(i,fix)
            #calculate total fracture volume
            if typ(self.faces[i].typ) != 'boundary':
                Vii += (4.0/3.0)*pi*0.25*self.faces[i].dia**2.0*0.5*self.faces[i].bd
                
        #fetch target flow rate
        if stim:
            qi = self.rock.Qstim
        else:
            qi = self.rock.Qinj
        
        #fetch maximum injection pressure
        if plim:
            pmax = self.rock.pfinal_max
        else:
            pmax = np.inf
            
        #set boundary conditions with excess flow check
        q_well = np.full(len(self.wells),None)
        p_well = np.full(len(self.wells),None)
        #set injection boundary conditions
        for i in range(0,len(i_key)):
            # if (self.Qi[i] > 0.95*qi) and (tryflowbound):
            #     q_well[i_key[i]] = qi
            # else:
            #     p_well[i_key[i]] = np.min([tip[i],pmax])
            p_well[i_key[i]] = np.min([tip[i],pmax])

        #set production boundary conditions
        for i in range(0,len(p_key)):
            if not(stim):
                p_well[p_key[i]] = self.rock.BH_P + self.rock.p_whp
            else:
                q_well[p_key[i]] = 0.0
        
        #solve flow
        self.get_flow(p_bound=self.rock.BH_P,p_well=p_well,q_well=q_well,reinit=reinit,Qnom=qi)
        Pii = []
        for i in range(0,len(i_key)):
            source = self.wells[i_key[i]].c0
            ck, j = self.nodes.add(source)
            Pii += [self.nodes.p[j]]
        Pii = np.asarray(Pii)    
        Qii = self.p_q[i_key]
        self.Qi = Qii
        self.Pi = Pii
        self.Vi = Vii
        
        #stabilize flow without recomputing mechanics
        flowbound = False
        if tryflowbound:
            #inform of pre-solution *** debugging
            print('-> flow solution prior to flow boundary')
            print('P\t\tQ')
            print('%.2e\t%.2e' %(self.Pi[0],self.Qi[0]))
            #set boundary conditions with excess flow check
            q_well = np.full(len(self.wells),None)
            p_well = np.full(len(self.wells),None)
            #set injection boundary conditions
            for i in range(0,len(i_key)):
                if (self.Qi[i] > 0.95*qi):
                    flowbound = True
                    q_well[i_key[i]] = qi
                else:
                    p_well[i_key[i]] = np.min([tip[i],pmax])
    
            #set production boundary conditions
            for i in range(0,len(p_key)):
                if not(stim):
                    p_well[p_key[i]] = self.rock.BH_P + self.rock.p_whp
                else:
                    q_well[p_key[i]] = 0.0
            
            #solve flow
            if flowbound:
                self.get_flow(p_bound=self.rock.BH_P,p_well=p_well,q_well=q_well,reinit=reinit,Qnom=qi)
            
        
        #3D visuals of propagation process
        self.progress += [[self.faces]]
        if visuals:
            fname2 = 'stim_%i_%i' %(self.iters,self.pin)
            self.build_vtk(fname2,vtype=[0,0,1,1,1,0])
        
        #TODO: add seismicity tracker
        
        #record injection rates and pressures
        Pii = []
        for i in range(0,len(i_key)):
            source = self.wells[i_key[i]].c0
            ck, j = self.nodes.add(source)
            Pii += [self.nodes.p[j]]
        Pii = np.asarray(Pii)    
        Qii = self.p_q[i_key]
        self.Qi = Qii
        self.Pi = Pii
        self.Vi = Vii
        # print(self.Pi)
        # print(self.Qi)
        return stim_num

    # ************************************************************************
    # stimulation - add frac (Vinj is total volume per stage, sand ratio to frac slurry by volume)
    # ************************************************************************
    def solve(self, visuals=True):
        print( '*** geodt solver module ***')
        #working variables
        iters = 0
        maxit = 40
        Pis = [] #Pa
        Qis = [] #m3/s
        timesteps = []
        
        #get index of key wells (index will match p_q and i_q from flow solver outputs)
        i_div = 0
        i_key = [] #injectors
        p_div = 0
        p_key = [] #producers
        s_div = 0
        s_key = [] #perforated zones
        for i in range(0,len(self.wells)):
            if (int(self.wells[i].typ) in [typ('injector')]):
                i_key += [int(i)]
                i_div += 1
            elif (int(self.wells[i].typ) in [typ('producer')]):
                p_key += [int(i)]
                p_div += 1
            elif (int(self.wells[i].typ) in [typ('perfcluster')]):
                s_key += [int(i)]
                s_div += 1
        i_key = np.asarray(i_key,dtype=int)
        p_key = np.asarray(p_key,dtype=int)
        s_key = np.asarray(s_key,dtype=int)
        
        #volume tracker variables
        vol_ini = 0.0 #m3
        vol_new = 0.0 #m3
        vol_old = 0.0 #m3
        
        #initial fracture parameters and network volume
        self.re_init()
        self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
        for i in range(0,len(self.faces)):
            #stabilize critically weak fractures
            self.faces[i].check_integrity(rock=self.rock,pres=(self.rock.BH_P+0.5*self.rock.dPi))
            #initial conditions for fractures
            self.faces[i].Pmax = self.rock.BH_P
            self.faces[i].Pcen = self.rock.BH_P
            self.hydromech(i)
            if typ(self.faces[i].typ) != 'boundary':
                vol_old += (4.0/3.0)*pi*0.25*self.faces[i].dia**2.0*0.5*self.faces[i].bd
        vol_ini = vol_old
        
        #initalization
        if len(self.Qi) == 0:
            self.Qi = np.zeros(len(i_key))
            
        #### stimulation ************
        if not(self.rock.strategy in ['AGS']):
            print( '\n+> stimulation using %.2e m3/s, %.2e m3 water, %.2e m3/m3 sand' %(self.rock.Qstim,self.rock.Vstim,self.rock.sand))
            
            #intermediate variables
            completed = np.zeros(i_div,dtype=bool) #T/F segment stimulation complete
            tip = np.ones(i_div,dtype=float)*(self.rock.s3 - self.rock.dPi) #trial injection pressure
            
            #stimulation loop -- production wells shut in
            # iters = 0
            # maxit = 40
            # Pis = [] #Pa
            # Qis = [] #m3/s
            # timesteps = []
            while 1:
                #loop counter and breaker
                if iters >= maxit:
                    break
                iters += 1
                print( '\n[%i] rock stimulation step' %(iters))
                
                #update test injection pressures
                for i in range(0,i_div):
                    #set to boundary pressure if interval was completed
                    if completed[i]:
                        tip[i] = self.rock.BH_P
                    #set with pressure incrementer if not completed
                    else:
                        tip[i] += self.rock.dPi
                print( '   + injection pressure increased to %.0f' %(np.max(tip)))
                
                #complete if pressure driven injection rate exceeds stimulation rate (+=in, -=out)
                for i in range(0,i_div):
                    if self.Qi[i] > self.rock.Qstim:
                        completed[i] = True
                if np.sum(completed) == i_div:
                    print( '   % stimulation complete by sustained flow')
                    break
                
                #solve hydromechanics step
                stim_num = self.solve_hydromech(i_key,p_key,tip,reinit=True,fix=False,
                            stim=True,plim=False,tryflowbound=False,visuals=visuals)
                vol_new = self.Vi
    
                #stimulation workflow
                timestep = [] 
                for i in range(0,i_div):
                    #calculate timestep by volume change
                    timestep += [(vol_new-vol_old)/(self.rock.Qstim-self.Qi[i])]
                    
                    #stimulate hydraulic fractures if criteria met
                    if not(completed[i]) and (self.wells[i_key[i]].hydrofrac==False) and (tip[i]>(self.rock.s3+self.rock.hfmcc)):
                        print( '   ~ (%i) interval hydraulic fractures' %(i))
                        self.wells[i_key[i]].hydrofrac = True
                        #seed hydraulic fracture, note that hydrofractures are in perfclusters
                        self.gen_stimfracs(target=s_key[i],perfs=self.rock.perf_clusters,
                                      f_dia = [2.0*self.rock.r_perf,0.0],
                                      f_azn = [self.rock.s3Azn+np.pi/2.0, self.rock.s3AznVar],
                                      f_dip = [np.pi/2.0-self.rock.s3Dip, self.rock.s3DipVar],
                                      clear=False)
                        #ensure stability of the hydraulic fracture (we don't want supercritical hydrofracs)
                        for notch in range(0,self.rock.perf_clusters):
                            fi = len(self.hydfs)-notch-1
                            si = self.hydfs[fi].sn
                            self.hydfs[fi].make_critical(rock=self.rock,pres=(si+self.rock.hfmcc))
                timesteps += [[timestep]]
                
                #complete if injected volume exceeds Vstim
                if (vol_new - vol_ini) > self.rock.Vstim:
                    print( '   - stimulation complete by injected volume')
                    #correct for excess volume
                    if (vol_new > vol_ini):
                        excess = self.rock.Vstim/(vol_new - vol_ini)
                    else:
                        excess = 1.0
                    for i in range(0,len(self.faces)):
                        self.faces[i].prop_load = self.faces[i].prop_load * excess
                    break
                else:
                    print( '   - stimulation volume remaining %.2e m3' %(self.rock.Vstim-vol_new+vol_ini))           
                                
                #update fracture network volume
                vol_old = vol_new
            
        #### circulation ************
        print( '\n+> circulation using %.2e m3/s, %.2e m3 water, %.2e Pa pinj_max' %(self.rock.Qinj,self.rock.Vinj,self.rock.pfinal_max))
        
        #disable sand
        remember_sand = self.rock.sand
        self.rock.sand = 0.0

        #intermediate variables
        completed = np.zeros(i_div,dtype=bool) #T/F segment stimulation complete
        tip = np.ones(i_div,dtype=float)*(self.rock.s3 - self.rock.dPi) #trial injection pressure
        
        #reset fracture parameters
        self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
        for i in range(0,len(self.faces)):
            #initial conditions for fractures
            self.faces[i].Pmax = self.rock.BH_P
            self.faces[i].Pcen = self.rock.BH_P
            self.hydromech(i)
            if typ(self.faces[i].typ) != 'boundary':
                vol_old += (4.0/3.0)*pi*0.25*self.faces[i].dia**2.0*0.5*self.faces[i].bd
        
        #circulation loop -- production wells flowing
        Pis = [] #Pa
        Qis = [] #m3/s
        while 1:
            #loop counter and breaker
            if iters >= maxit:
                break
            iters += 1
            print( '\n[%i] rock stimulation step' %(iters))
            
            #update test injection pressures
            for i in range(0,i_div):
                #set to boundary pressure if interval was completed
                if completed[i]:
                    tip[i] = self.rock.BH_P
                #set with pressure incrementer if not completed
                else:
                    tip[i] += self.rock.dPi
            print( '   + injection pressure increased to %.0f' %(np.max(tip)))

            #complete if pressure driven injection rate exceeds stimulation rate (+=in, -=out)
            for i in range(0,i_div):
                if self.Qi[i] > self.rock.Qinj:
                    completed[i] = True
            if np.sum(completed) == i_div:
                print( '   % full circulation rate successfully acheived')
                break
            
            #solve hydromechanics step
            stim_num = self.solve_hydromech(i_key,p_key,tip,reinit=True,fix=False,
                        stim=False,plim=True,tryflowbound=False,visuals=visuals)
            vol_new = self.Vi
                
            #stimulation workflow
            timestep = [] 
            for i in range(0,i_div):
                #calculate timestep by volume change
                timestep += [(vol_new-vol_old)/(self.rock.Qinj-self.Qi[i])]
                
                #stimulate hydraulic fractures if criteria met
                if not(completed[i]) and (self.wells[i_key[i]].hydrofrac==False) and (tip[i]>(self.rock.s3+self.rock.hfmcc)):
                    print( '   ~ (%i) interval hydraulic fractures' %(i))
                    self.wells[i_key[i]].hydrofrac = True
                    #seed hydraulic fracture, note that hydrofractures are in perfclusters
                    self.gen_stimfracs(target=s_key[i],perfs=self.rock.perf_clusters,
                                  f_dia = [2.0*self.rock.r_perf,0.0],
                                  f_azn = [self.rock.s3Azn+np.pi/2.0, self.rock.s3AznVar],
                                  f_dip = [np.pi/2.0-self.rock.s3Dip, self.rock.s3DipVar],
                                  clear=False)
                    #ensure stability of the hydraulic fracture (we don't want supercritical hydrofracs)
                    for notch in range(0,self.rock.perf_clusters):
                        fi = len(self.hydfs)-notch-1
                        si = self.hydfs[fi].sn
                        self.hydfs[fi].make_critical(rock=self.rock,pres=(si+self.rock.hfmcc))
            timesteps += [[timestep]]
            
            #complete if injected volume exceeds Vstim
            if (vol_new - vol_ini) > self.rock.Vinj:
                print('   - injection fully accomodated by volume change')
                break
            else:
                print('   - circulation volume remaining %.2e m3' %(self.rock.Vinj-vol_new+vol_ini))
            
            #complete if no stimulation at maximum pressure
            if (stim_num == 0) and (np.max(tip) > self.rock.pfinal_max):
                print('   ! stimulation reached maximum pressure limit without propagating fractures')
                break
            
            #update fracture network volume
            vol_old = vol_new
            
        #### hydromechanics ************
        print( '\n+> circulation flow solver with static hydromechanics')
        #working variables
        Qis = [] #array of flows
        Pis = [] #array of pressures
        
        #***** [0] ambient flow (linear)
        print( '\n[0] Ambient pressure')
        #reset nodal pressures
        self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
        #calculate injection rate and pressures
        tip = self.rock.BH_P*np.ones(i_div)
        #1st solve - low flow
        self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                        stim=False,plim=False,tryflowbound=False,visuals=visuals)
        #2nd solve - high flow
        self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                        stim=False,plim=False,tryflowbound=False,visuals=visuals)
        Pis += [self.Pi]
        Qis += [self.Qi]

        #***** [1] flow at critical pressure (pseudo-linear)
        print( '\n[1] Critical pressure')
        #reset nodal pressures
        self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
        #calculate injection rate and pressures
        tip = self.rock.s3*np.ones(i_div) - self.rock.dPi
        #1st solve - low flow
        self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                        stim=False,plim=False,tryflowbound=False,visuals=visuals)
        #2nd solve - high flow
        self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                        stim=False,plim=False,tryflowbound=False,visuals=visuals)
        Pis += [self.Pi]
        Qis += [self.Qi]
        
        #***** [2] flow at over-critical pressure (highly non-linear)
        print( '\n[2] Overcritical pressure')
        #reset nodal pressures
        self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
        #calculate injection rate and pressures
        tip = self.rock.s3*np.ones(i_div) + self.rock.hfmcc + 2*39.0*self.rock.dPi
        #1st solve - low flow
        self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                        stim=False,plim=False,tryflowbound=False,visuals=visuals)
        #2nd solve - high flow
        self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                        stim=False,plim=False,tryflowbound=False,visuals=visuals)
        Pis += [self.Pi]
        Qis += [self.Qi]

        #convert to array
        Pis = np.asarray(Pis)
        Qis = np.asarray(Qis)
        
        #error checking
        print('-> final flow iterations with hydromechanical coupling')
        print(' \tP\t\tQ')
        print('low\t%.2e\t%.2e' %(Pis[0][0],Qis[0][0]))
        print('mid\t%.2e\t%.2e' %(Pis[1][0],Qis[1][0]))
        print('hi \t%.2e\t%.2e' %(Pis[2][0],Qis[2][0]))
        
        #***** [3] pressure limited case
        skipnext = False
        if Pis[2][0] > self.rock.pfinal_max:
            print( '\n[3] Pressure limit upper')
            #reset nodal pressures
            self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
            #calculate injection rate and pressures
            tip = self.rock.s3*np.ones(i_div) + self.rock.hfmcc + 2*39.0*self.rock.dPi
            #1st solve - low flow
            self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                            stim=False,plim=True,tryflowbound=False,visuals=visuals)
            #2nd solve - high flow
            self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                            stim=False,plim=True,tryflowbound=False,visuals=visuals)
            Pis[2] = np.asarray(self.Pi)
            Qis[2] = np.asarray(self.Qi)
            if self.rock.Qinj > np.max(Qis[2]):
                skipnext = True
                
        if not(skipnext):        
            if Pis[1][0] >= Pis[2][0]:
                print( '\n[4] Pressure limit midpoint')
                #reset nodal pressures
                self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
                #calculate injection rate and pressures
                tip = 0.5*(Pis[2] + Pis[0])
                #1st solve - low flow
                self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                                stim=False,plim=True,tryflowbound=False,visuals=visuals)
                #2nd solve - high flow
                self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                                stim=False,plim=True,tryflowbound=False,visuals=visuals)
                Pis[1] = np.asarray(self.Pi)
                Qis[1] = np.asarray(self.Qi)            
            
        
            #***** [5,6,7,8,9] search for best solution using binary search
            for p in range(0,5):
                print('[%i] Hydromechanics binary search' %(p+5))
                #get next test pressures per interval
                a = np.zeros(i_div,dtype=int) + (Qis[1] > self.rock.Qinj)
                Pis[0] = Pis[0]*(a) + Pis[1]*(1-a)
                Pis[2] = Pis[1]*(a) + Pis[2]*(1-a)
                Qis[0] = Qis[0]*(a) + Qis[1]*(1-a)
                Qis[2] = Qis[1]*(a) + Qis[2]*(1-a)
                Pis[1] = 0.5*(Pis[0] + Pis[2])
                #reset nodal pressures
                self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
                #solve for flow
                tip = Pis[1]
                #1st solve - low flow
                self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                                stim=False,plim=True,tryflowbound=False,visuals=visuals)
                #2nd solve - high flow
                self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                                stim=False,plim=True,tryflowbound=False,visuals=visuals)
                Pis[1] = np.asarray(self.Pi)
                Qis[1] = np.asarray(self.Qi)
                #error checking
                print('-> final flow iterations with hydromechanical coupling')
                print(' \tP\t\tQ')
                print('low\t%.2e\t%.2e' %(Pis[0][0],Qis[0][0]))
                print('mid\t%.2e\t%.2e' %(Pis[1][0],Qis[1][0]))
                print('hi \t%.2e\t%.2e' %(Pis[2][0],Qis[2][0]))
        
        #***** [10] flow boundary if flow is excessive, but don't calculate changes in apertures this step
        #reset nodal pressures
        self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
        #solve for flow
        tip = Pis[2]
        #1st solve - low flow
        self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                        stim=False,plim=True,tryflowbound=False,visuals=visuals)
        #2nd solve - high flow
        self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                        stim=False,plim=True,tryflowbound=True,visuals=visuals)
        #error checking
        print('-> flow solution completed')
        print('P\t\tQ')
        print('%.2e\t%.2e' %(self.Pi[0],self.Qi[0]))
        #override to compensate for hydropropped fracture instability
        for i in range(0,i_div):
            source = self.wells[i_key[i]].c0
            ck, j = self.nodes.add(source)
            self.nodes.p[j] = Pis[2][i]
        for w in range(0,len(self.wells)):
            #coordinates
            source = self.wells[w].c0
            #find index of duplicate
            ck, i = self.nodes.add(source)
            #record pressure
            self.p_p[w] = self.nodes.p[i]
            
        #remember sand
        self.rock.sand = remember_sand
    
    # ************************************************************************
    # basic validation of fracture geometry
    # ************************************************************************
    def detournay_visc(self,Q0=0.08):
        #plot initialization
        fig = pylab.figure(figsize=(15.0, 8.5), dpi=96, facecolor='w', edgecolor='k',tight_layout=True) # Medium resolution
        font = {'family' : 'serif',   'size'   : 16}
        pylab.rc('font', serif='Arial')
        pylab.rc('font', **font)
        ax1 = fig.add_subplot(231)
        ax2 = fig.add_subplot(232)
        ax3 = fig.add_subplot(234)
        ax4 = fig.add_subplot(235)
        ax5 = fig.add_subplot(236)
        
        #normalized constants
        E_prime = self.rock.ResE/(1.0-self.rock.Resv**2.0)
        K_prime = 4.0*self.rock.Kic*(2.0/pi)**0.5
        mu_prime = 12.0*self.rock.Poremu
        
        #time
        time = np.linspace(1000.0,200000.0,100,dtype=float)
        rho = np.linspace(0.001,0.999,100,dtype=float)
        drho = rho[1] - rho[0]
        wst = []
        pst = []
        rst = []
        Vt = []
        Pt = []
        
        for t in time:
            #radius
            sca_radi = 0.0
            sca_Kscas = []
            sca_apers = []
            sca_press = []
            for p in range(0,len(rho)):
                #viscosity dominated if K_scaled < 1.0
                K_scaled = K_prime*(t**2.0/(mu_prime**5.0 * Q0**3.0 * E_prime**13.0))**(1.0/18.0)
                
                #scaled aperture (ohmega_bar_m0) and scaled pressure (pi_m0)
                scaled_aper = (((0.6846*70.0**0.5)/3.0  + (13.0*rho[p]-6.0)*(0.07098*4.0*5.0**0.5)/9.0)*(1.0 - rho[p])**(2.0/3.0) + 
                               0.09269*((1.0-rho[p])**0.5*8.0/pi - np.arccos(rho[p])*rho[p]*8.0/pi))
                scaled_pres = 0.3581*(2.479 - 2.0/(3.0*(1.0-rho[p])**(1.0/3.0))) - 0.09269*(np.log(rho[p]/2.0)+1.0)
                
                sca_Kscas += [K_scaled]
                sca_apers += [scaled_aper]
                sca_press += [scaled_pres]
                
                #scaled radius (gamma_m0)
                if p != 0:
                    sca_radi += 0.5*(scaled_aper + sca_apers[p-1])*0.5*(rho[p] + rho[p-1])*drho
            sca_radi = (2.0*pi*sca_radi)**(-1.0/3.0)

            #convert to arrays for math
            sca_Kscas = np.asarray(sca_Kscas)
            sca_apers = np.asarray(sca_apers)
            sca_press = np.asarray(sca_press)
            
            #less scaled aperture (ohmega_m0)
            n_apers = sca_apers*sca_radi
            
            #lets now undo this irritating scaling crap
            Lm = (E_prime*Q0**3.0*t**4.0/mu_prime)**(1.0/9.0)
            em = (mu_prime/(E_prime*t))**(1.0/3.0)
            ws = em*Lm*n_apers
            ps = em*E_prime*sca_press
            rs = Lm*rho*sca_radi
            dr = Lm*drho*sca_radi
            
            #store for plotting
            wst += [ws]
            pst += [ps]
            rst += [rs]
            
            #now for what we actually care about: volume and radius
            vol = 0.0
            pnet = 0.0
            for i in range(1,len(rs)):
                vol += 2.0*pi*(0.5*(rs[i-1]+rs[i]))*(0.5*(ws[i-1]+ws[i]))*dr
                pnet += 2.0*pi*(0.5*(rs[i-1]+rs[i]))*(0.5*(ps[i-1]+ps[i]))*dr
            pnet = pnet/(pi*rs[-1]**2.0)
            Vt += [vol]
            Pt += [pnet]
            
            #plot
            ax1.plot(rs,ps)
            ax2.plot(rs,ws)
        #ax3.plot(time,np.asarray(Vt))
        #ax4.plot(time,np.asarray(Pt))
        Vt = np.asarray(Vt)
        Pt = np.asarray(Pt)
        rst = np.asarray(rst)
        ax3.plot(rst[:,-1],Vt)
        ax4.plot(rst[:,-1],Pt)
        ax5.plot(rst[:,-1],Vt/(pi*rst[:,-1]**2.0))
        
        #labels
        ax1.set_xlabel('Radial Distance (m)')
        ax1.set_ylabel('Pressure (Pa)')
        
        ax2.set_xlabel('Radial Distance (m)')
        ax2.set_ylabel('Aperture (m)')

        ax3.set_xlabel('Radius (m)')
        ax3.set_ylabel('Volume (m3)')
        
        ax4.set_xlabel('Radius (m)')
        ax4.set_ylabel('Avg. Net Pressure (Pa)')
        
        ax5.set_xlabel('Radius (m)')
        ax5.set_ylabel('Avg. Aperture (m)')

    # ************************************************************************
    # optimization objective function (with $!!!), to be normalized to 2021 studies
    # ************************************************************************
    def get_economics(self,
                      detail=False, #print cycle infomation
                      plots=False): #plot results
        #eliminate periods of negative net power production
        keep = (self.Pout > 0.0)
        NPsum = np.sum(self.Pout[keep])
        TPsum = np.sum(self.Tout[keep])
        QPsum = np.sum(self.Qout[keep])
        HPsum = np.sum(self.Hout[keep])
        count = len(self.Pout[keep])
        if count > 0:
            self.payoff.Power_Net_kW = NPsum/len(self.Pout[keep])
            self.payoff.Power_Bulk_kW = TPsum/len(self.Tout[keep])
            self.payoff.Power_Pump_kW = QPsum/len(self.Qout[keep])
            self.payoff.Power_Thermal_kW = HPsum/len(self.Hout[keep])
        else:
            self.payoff.Power_Net_kW = 0
            self.payoff.Power_Bulk_kW = 0
            self.payoff.Power_Pump_kW = 0
            self.payoff.Power_Thermal_kW = 0
        #get system values
        dt = (self.ts[1]-self.ts[0])/yr
        life = self.rock.LifeSpan/yr
        depth = self.rock.ResDepth
        drill_len = 0.0
        for w in range(0,len(self.wells)): #when all the wells are in the model
            drill_len += self.wells[w].leg
        if drill_len < 2*self.rock.ResDepth: #wells were not placed to surface, estimate another way
            lateral = (self.rock.w_length*self.rock.w_count) + (self.rock.w_proportion*self.rock.w_length)
            drill_len = self.rock.ResDepth*(self.rock.w_count+1) + lateral
        Max_Quake = -10.0
        Num_Quake = 0
        for i in range(0,len(self.faces)):
            num = len(self.faces[i].Mws)
            if num > 0:
                Num_Quake += num
                Max_Quake = np.max([Max_Quake,np.max(self.faces[i].Mws)])
        #power profit
        P = self.rock.sales_kWh*NPsum*dt*24.0*365.2425
        #bulk power sales
        B = self.rock.sales_kWh*TPsum*dt*24.0*365.2425
        #capital costs
        C = 0.0
        C += self.rock.drill_m*drill_len
        C += self.rock.pad_fixed
        C += np.max([0.0,self.rock.plant_kWe*TPsum*dt/life])
        C += self.rock.explore_m*depth
        #quake cost
        Q = self.rock.quake_coef * np.exp(Max_Quake*self.rock.quake_exp)
        #generation cost
        G = self.rock.oper_kWh*TPsum*dt*24.0*365.2425
        #net present value
        NPV = P - C - Q - G
        #cost per MWe-h net electic power production
        CpM = 0.0
        if NPsum <= 0:
            CpM = np.inf
        else:
            CpM = (C + Q + G)/(10**-3*NPsum*dt*24.0*365.2425)
        #cost per MWe-h production without pumping losses
        CpM_prod = 0.0
        if TPsum <= 0:
            CpM_prod = np.inf
        else:
            CpM_prod = (C + Q + G)/(10**-3*TPsum*dt*24.0*365.2425)
        #cost per MWt-h bulk extraction
        CpM_ther = 0.0
        if HPsum <= 0:
            CpM_ther = np.inf
        else:
            CpM_ther = (C + Q)/(10**-3*HPsum*dt*24.0*365.2425)
        #detail
        if detail:
            print('\n*** economics module ***')
            print('   sales: $%.0f (%.2f kWh)' %(P,NPsum*dt*24.0*365.2425))
            print('   capital: $%.0f (%.2f m drilled length)' %(C,drill_len))
            print('   seismic: $%.0f (%.2f Mw max quake)' %(Q,Max_Quake))
            print('   net present value (npv): $%.0f' %(NPV))
            print('   cost per megawatt-hour: $%.2f/MWe-h' %(CpM))
        
        self.eNPV = NPV
        self.eC = C
        self.eQ = Q
        self.eG = G
        self.eP = P
        self.eB = B
        self.CpM_elec = CpM #cost per net MWe-h
        self.CpM_prod = CpM_prod #cost per produced MWe-h
        self.CpM_ther = CpM_ther #cost per MWt-h
        self.drill_length = drill_len #calculated drilled length
        
        self.payoff.Economics_NetPresentValue_USD = NPV
        self.payoff.Economics_Capital_USD = C
        self.payoff.Economics_Upkeep_USD = G
        self.payoff.Economics_QuakeRisk_USD = Q
        self.payoff.Economics_PowerProfit_USD = P
        self.payoff.Economics_BulkSales_USD = B
        self.payoff.Seismicity_MaxQuake_Mw = Max_Quake
        self.payoff.Seismicity_NumQuake_ea = Num_Quake
        self.payoff.LCOE_Net_USDpMWh = CpM
        self.payoff.LCOE_Bulk_USDpMWh = CpM_prod
        self.payoff.LCOE_Thermal_USDpMWh = CpM_ther
        self.payoff.Drilling_TotalLength_m = drill_len

        return NPV, P-G, C, Q

    # ************************************************************************
    # calculate flow rate as a function of pressure using GeoDT
    # ************************************************************************
    def rate_vs_pressure(self, tips=[], points=40, visuals=True):
        #working variables
        Pis = [] #Pa
        Qis = [] #m3/s
        
        #disable sand
        remember_sand = self.rock.sand
        self.rock.sand = 0.0
        
        #get index of key wells (index will match p_q and i_q from flow solver outputs)
        i_div = 0
        i_key = [] #injectors
        p_div = 0
        p_key = [] #producers
        s_div = 0
        s_key = [] #perforated zones
        for i in range(0,len(self.wells)):
            if (int(self.wells[i].typ) in [typ('injector')]):
                i_key += [int(i)]
                i_div += 1
            elif (int(self.wells[i].typ) in [typ('producer')]):
                p_key += [int(i)]
                p_div += 1
            elif (int(self.wells[i].typ) in [typ('perfcluster')]):
                s_key += [int(i)]
                s_div += 1
        i_key = np.asarray(i_key,dtype=int)
        p_key = np.asarray(p_key,dtype=int)
        s_key = np.asarray(s_key,dtype=int)
        
        #initalization
        Qis = [] #array of flows
        Pis = [] #array of pressures
        Vis = [] #array of volumes
        Pcl = [] #array of median perf cluster pressure
        Pch = [] #array of median pressure after perf choke
        Psc = [] #array of median screen pressure
        Pba = [] #array pf backpressures
        Ps3 = [] #array of min-stress
        
        #rate from pressure loop
        if tips == []:
            pmax = 70.0*MPa
            pmin = pmax/100.0
            tips = 10.0**np.linspace(np.log10(pmin),np.log10(pmax),points) + self.rock.BH_P + self.rock.p_whp
            #TODO: add critical pressure point
        for i in range(0,points):
            # #critical pressure corner point
            # if (tips[i] > (self.rock.s3 - self.rock.dPi)) and (not(corner)):
            #     corner = True
            #     #reset nodal pressures
            #     self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
            #     #calculate injection rate and pressures
            #     tip = self.rock.s3*np.ones(i_div)
            #     #1st solve - low flow
            #     self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
            #                     stim=False,plim=False,tryflowbound=False,visuals=visuals)
            #     #2nd solve - high flow
            #     self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
            #                     stim=False,plim=False,tryflowbound=False,visuals=visuals)
            #     #get intermediate pressure points
            #     p_clusters = []
            #     p_chokes = []
            #     p_screens = []
            #     for k in range(0,self.pipes.num):
            #         if typ(self.pipes.typ[k]) == 'perfcluster':
            #             p_clusters += [self.nodes.p[self.pipes.n1[k]]]
            #         elif typ(self.pipes.typ[k]) == 'choke':
            #             p_chokes += [self.nodes.p[self.pipes.n1[k]]]    
            #         elif typ(self.pipes.typ[k]) == 'screen':
            #             p_screens += [self.nodes.p[self.pipes.n1[k]]]
            #     #store data
            #     Pis += [self.Pi]
            #     Qis += [self.Qi]
            #     Vis += [self.Vi]
            #     Pcl += [np.median(p_clusters)]
            #     Pch += [np.median(p_chokes)]
            #     Psc += [np.median(p_screens)]
            #     Pba += [self.rock.p_whp]
            #     Ps3 += [self.rock.s3]
                
            #nominal point
            #reset nodal pressures
            self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
            #calculate injection rate and pressures
            tip = tips[i]*np.ones(i_div)
            #1st solve - low flow
            self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                            stim=False,plim=False,tryflowbound=False,visuals=visuals)
            #2nd solve - high flow
            self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
                            stim=False,plim=False,tryflowbound=False,visuals=visuals)
            ####TODO
            '''
            #### pressure corrections
            i_prop = np.zeros(i_div) #index map of wells that are upstream of hydroprop
            for b in range(0,self.pipes.num):
                #start search from fractures and end with wells
                p = self.pipes.num - b
                if typ(self.pipes.typ[p] in ['choke','perf','propped','fracture']):
                    if self.faces[self.pipes.fID[p]].hydroprop:
                        #identify minimum pressure
                        #calculate offset to impose hydroprop
                        #add critical pressure
                elif i_prop[self.pipes.pID[p]] > 0: #check if parent well is upstream of hydroprop
                    pass
                
            
            for f in range(6,len(self.faces)):
                if self.faces[f].hydroprop == True:
                    f_prop[f] self.faces[f].hydroprop = True
                    break
            if hydroprop:
                #identify parent well
            '''
                
            
            #get intermediate pressure points
            p_clusters = []
            p_chokes = []
            p_screens = []
            for k in range(0,self.pipes.num):
                if typ(self.pipes.typ[k]) == 'perfcluster':
                    p_clusters += [self.nodes.p[self.pipes.n1[k]]]
                elif typ(self.pipes.typ[k]) == 'choke':
                    p_chokes += [self.nodes.p[self.pipes.n1[k]]]    
                elif typ(self.pipes.typ[k]) == 'screen':
                    p_screens += [self.nodes.p[self.pipes.n1[k]]]
            #store data
            Pis += [self.Pi]
            Qis += [self.Qi]
            Vis += [self.Vi]
            Pcl += [np.median(p_clusters)]
            Pch += [np.median(p_chokes)]
            Psc += [np.median(p_screens)]
            Pba += [self.rock.p_whp]
            Ps3 += [self.rock.s3]
                
        #convert to array
        Pis = np.asarray(Pis)-self.rock.BH_P
        Qis = np.asarray(Qis)
        Vis = np.asarray(Vis)-Vis[0]
        Pcl = np.asarray(Pcl)-self.rock.BH_P
        Pch = np.asarray(Pch)-self.rock.BH_P
        Psc = np.asarray(Psc)-self.rock.BH_P
        Pba = np.asarray(Pba)
        Ps3 = np.asarray(Ps3)-self.rock.BH_P
        
        # #***** [10] flow boundary if flow is excessive, but don't calculate changes in apertures this step
        # #reset nodal pressures
        # self.nodes.p = self.rock.BH_P + 1.0*MPa*np.random.rand(self.nodes.num)
        # #solve for flow
        # tip = Pis[2]
        # #1st solve - low flow
        # self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
        #                 stim=False,plim=True,tryflowbound=False,visuals=visuals)
        # #2nd solve - high flow
        # self.solve_hydromech(i_key,p_key,tip,reinit=False,fix=True,
        #                 stim=False,plim=True,tryflowbound=True,visuals=visuals)
        # #error checking
        # print('-> flow solution completed')
        # print('P\t\tQ')
        # print('%.2e\t%.2e' %(self.Pi[0],self.Qi[0]))
        # #override to compensate for hydropropped fracture instability
        # for i in range(0,i_div):
        #     source = self.wells[i_key[i]].c0
        #     ck, j = self.nodes.add(source)
        #     self.nodes.p[j] = Pis[2][i]
            
        #remember sand
        self.rock.sand = remember_sand
        
        return Pis, Qis, Vis, Pcl, Pch, Psc, Pba, Ps3

# Visualization tools
class visualization_SP:
    def __init__(self,filename='setup_payoff.csv'):
        #import data and store unaltered copy
        try:
            self.names, self.data = load_csv(filename)
        except:
            geom = core()
            setu = setup()
            geom.rock.fetch(setu)
            geom.setup_payoff('setup_payoff_placeholder.csv',geom.pin,aux=[],append=False)
            self.names, self.data = load_csv('setup_payoff_placeholder.csv')
        # if np.size(self.data) > 1:
        self.orig = np.copy(self.data)
        self.normalize()
    def normalize(self):
        #normalize data
        self.data_n = np.copy(self.data)
        for n in self.names:
            try:
                np.nan_to_num(self.data_n[n],copy=False,nan=0.0)
                #TODO: Add special cases
                # if n == 'Well_RotationSkew_rad':
                #     self.data_n[n] = np.abs(self.data_n[n])
                span = np.max(self.data_n[n])-np.min(self.data_n[n])
                if span > 0:
                    self.data_n[n] = (self.data_n[n]-np.min(self.data_n[n]))/span
                else:
                    self.data_n[n] = 0.0
            except:
                self.data_n[n] = np.zeros(len(self.data[n]))
    #2 axis scatter plot with color per 3rd axis
    def scatter_3D(self,x,y,c,name=['x','y','c','save'],cmap='rainbow_r',
                   xlog=False,ylog=False,vrange=[],view=True,ax=[]):
        if not(ax):
            fig = pylab.figure(figsize=(10.0,5.0),dpi=100)
            ax1 = fig.add_subplot(111)
        else:
            ax1 = ax
        if vrange:
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
    def timeseries(self,ts,ys,name=['t','y','save'],
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
    #data filter
    def multifilter(self,target='pin',compare=None,vs=[0.0,1e12]):
        keep = np.ones(len(self.data),dtype=bool)
        try:
            if compare:
                keep = (self.data[target] < vs[1]*self.data[compare]) * (self.data[target] > vs[0]*self.data[compare])
            else:
                keep = (self.data[target] < vs[1]) * (self.data[target] > vs[0])
        except:
            print('%s not in setup_payoff' %(target))
        self.data = self.data[keep]
        self.data_n = self.data_n[keep]
    #reset
    def reset(self):
        self.data = np.copy(self.orig)
        self.normalize()
    #percentile
    def quantile(self,x,pcts=[],ax=[]):
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
    #percentiles
    def percentiles(self,filename='percentiles'):
        pcts = [0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99,100]
        first = True
        for p in pcts:
            out = [['Quantile','%i' %(p)]]
            out += [['Percentile','%i' %(100-p)]]
            for n in self.names:
                try:
                    pv = '%.3e' %(self.quantile(self.data[n],[p])[1][0])
                except:
                    pv = ''
                out += [[n, pv]]
            if first:
                save_csv(out,filename,False)
                first = False
            else:
                save_csv(out,filename,True)
            
def RUN(setup=setup()):
    geom = core() #create the model object
    geom.rock.fetch(setup) #take in model parameters
    geom.gen_domain() #place model boundaries
    geom.gen_joint_sets() #place natural fractures
    geom.gen_from_dfn(filename=geom.rock.setup.Strategy_DiscreteFractures_path) #place deterministic fractures
    geom.re_init() #model initialization
    geom.gen_wells(True,wells=[],style=geom.rock.strategy) #place wells
    # geom.build_vtk()
    geom.solve(False) #stimulate, circulate, and solve hydromechanics
    geom.get_heat(plot=False,detail=False,lapse=False)
    NPV, P, C, Q = geom.get_economics(detail=True) #calculate economics
    aux = []
    if setup.Strategy_TargetDirectory_path  != '':
        filename = setup.Strategy_TargetDirectory_path + '\\setup_payoff.csv'
    else:
        filename = 'setup_payoff.csv'
    print('************* saving file *************')
    print(filename)
    geom.setup_payoff(filename,geom.pin,aux=aux)
    packup(geom,geom.pin)

# def RUN2(setup=setup()):
if False:
    geom = core() #create the model object
    setting = setup()
    # setting.Stimulation_TargetVolume_m3 = 1.0
    # setting.Circulation_AllowableOverpressure_Pa = 140*MPa
    setting.Strategy_Design_type = "CGS"
    geom.rock.fetch(setting)
    geom.pin = np.random.randint(100000000,999999999)
    geom.gen_domain() #place model boundaries
    geom.gen_joint_sets()#place natural fractures
    geom.re_init() #model initialization
    print(geom.rock.pfinal_max)
    print(geom.rock.s3)
    geom.gen_wells(True,wells=[],style=geom.rock.strategy) #place wells
    geom.solve(True)
    # geom.get_heat(plot=False,detail=False,lapse=False)
    geom.get_heat(plot=True,detail=False,lapse=False)
    plt.savefig('plt_%i.png' %(geom.pin), format='png')
    plt.close()
    NPV, P, C, Q = geom.get_economics(detail=True) #calculate economics
    aux = []
    #geom.save('inputs_results.txt',pin,aux=aux,printwells=0,time=True) #save data
    geom.setup_payoff('setup_payoff.csv',geom.pin,aux=aux)
    packup(geom,geom.pin)


# ****************************************************************************
#### examples
# ****************************************************************************







