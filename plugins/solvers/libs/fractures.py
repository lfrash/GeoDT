# -*- coding: utf-8 -*-
""" 
Functions to create and draw fractures
"""
#initializations
import numpy as np
if __package__ is None or __package__ == '':
    import iogt
    import stats
    import stress
    import wells
    import vtk as sg
    import csv
    from typ import typ
    from units import *
else:
    from . import iogt
    from . import stats
    from . import stress
    from . import wells
    from . import vtk as sg
    from . import csv
    from .typ import typ
    from .units import *
deg=deg
yr=yr
cP=cP
darcy=darcy
mD=mD
g=g

# ********************************************************************
### surface object
# ********************************************************************
class surf:
    def __init__(self,c0=[0,0,0],dia=1.0,stk=0.0*deg,dip=90.0*deg,
                 ty='fracture', s=iogt.setup()):
        #node number of center point
        self.ci = -1
        self.wi = -1
        
        #geometry
        self.c0 = np.asarray(c0,dtype=float)
        self.dia = dia
        self.str = stk
        self.dip = dip
        self.typ = typ(ty)
        
        #properties
        self.En = s.Rock_YoungsModulus_Pa
        self.vn = s.Rock_PoissonRatio_ratio
        self.roughness = stats.loguniform(1,s.Fracture_Tortuosity_ratio,1.0)[0]
        self.alpha = s.Tension_FractureCompressibility_1pPa * (1.0 + stats.norm_trunc(1,0.0,0.05,-0.12,0.12)[0]*s.Strategy_SecondOrderNoise_ratio)
        self.prop_alpha = s.Stimulation_ProppantCompressibility_1pPa * (1.0 + stats.norm_trunc(1,0.0,0.05,-0.12,0.12)[0]*s.Strategy_SecondOrderNoise_ratio)
        self.kf = 10.0**(np.log10(s.Stimulation_ProppantPermeability_m2) * (1.0 + stats.norm_trunc(1,0.0,0.05,-0.12,0.12)[0]*s.Strategy_SecondOrderNoise_ratio))
        self.gamma = 10.0**(np.log10(s.Shear_SlipLength_ratio) * (1.0 + stats.norm_trunc(1,0.0,0.2,-0.4,0.5)[0]*s.Strategy_SecondOrderNoise_ratio))
        self.psi = s.Shear_DilationSlip_ratio * (1.0 + stats.norm_trunc(1,0.0,0.2,-0.46,0.35)[0]*s.Strategy_SecondOrderNoise_ratio) #a
        self.eta = s.Shear_HydraulicDilation_ratio * (1.0 + ((stats.contact_trunc(1,0.15,1.0,0.5,0.5,0.15,0.0,2.0))[0]-1.0)*s.Strategy_SecondOrderNoise_ratio) #n
        self.bd = s.JointSets_InitialDilation_m * (1.0 + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.Strategy_SecondOrderNoise_ratio)
        
        #friction
        if ty in ['propped']:
            self.mu = s.Tension_FrictionCoefficient_ratio * (1.0 + stats.norm_trunc(1,0.0,0.04,-0.08,0.10)[0]*s.Strategy_SecondOrderNoise_ratio)
            self.phi = np.arctan(self.mu)
            self.mcc = s.Tension_Cohesion_Pa * (1.0 + stats.norm_trunc(1,0.0,0.2,-0.5,1.5)[0]*s.Strategy_SecondOrderNoise_ratio)
            # self.bd = 0.0 #TODO: is this needed?
        else:
            self.mu = np.arctan(s.Shear_FrictionCoefficient_ratio * (1.0 + stats.norm_trunc(1,0.0,0.12,-0.28,0.36)[0]*s.Strategy_SecondOrderNoise_ratio))
            self.phi = np.arctan(self.mu)
            self.mcc = s.Shear_Cohesion_Pa * (1.0 + stats.norm_trunc(1,0.0,0.2,-0.5,1.5)[0]*s.Strategy_SecondOrderNoise_ratio)
        
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

# ********************************************************************
### saved discrete fracture network
# ********************************************************************
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
            out += [['HydraulicDilation_ratio',fracs[i].eta]]
            out += [['Compressibility_1pPa',fracs[i].alpha]]
            out += [['DilationSlip_ratio',fracs[i].psi]]
            out += [['SlipLength_ratio',fracs[i].gamma]]
            out += [['InitialHydraulicAperture_m',fracs[i].bh]]
            out += [['Conductivity_m2m',fracs[i].fc]]
            out += [['InitialDilation_m',fracs[i].bd]]
            if i < 1:
                csv.save_csv(out=out,filename=filename,append=False)
            else:
                csv.save_csv(out=out,filename=filename,append=True)
            
    #load fracture information and add to network
    def intake(self,filename='dfn.csv',s=iogt.setup()):
        #load data
        self.names, self.data = csv.load_csv(filename)
        fracs = []
        try:
            for i in range(0,len(self.data)):
                #create fracture
                c0 = np.asarray([self.data['Easting_m'][i],self.data['Northing_m'][i],self.data['Height_m'][i]])
                fracs += [surf(c0,self.data['Diameter_m'][i],self.data['Strike_rad'][i],self.data['Dip_rad'][i],typ(int(self.data['Type_type'][i]+0.001)),s)]
                
                #match fields
                if 'FrictionAngle_rad'            in self.names: fracs[-1].phi = self.data['FrictionAngle_rad'][i]
                if 'Cohesion_Pa'                  in self.names: fracs[-1].mcc = self.data['Cohesion_Pa'][i]
                if 'ProppantLoad_m3'              in self.names: fracs[-1].prop_load = self.data['ProppantLoad_m3'][i]
                if 'ProppantCompressibility_1pPa' in self.names: fracs[-1].prop_alpha = self.data['ProppantCompressibility_1pPa'][i]
                if 'Tortuosity_ratio'             in self.names: fracs[-1].roughness = self.data['Tortuosity_ratio'][i]
                if 'ProppantPermeability_m2'      in self.names: fracs[-1].kf = self.data['ProppantPermeability_m2'][i]
                if 'HydraulicDilation_ratio'      in self.names: fracs[-1].eta = self.data['HydraulicDilation_ratio'][i]
                if 'Compressibility_1pPa'         in self.names: fracs[-1].alpha = self.data['Compressibility_1pPa'][i]
                if 'DilationSlip_ratio'           in self.names: fracs[-1].psi = self.data['DilationSlip_ratio'][i]
                if 'SlipLength_ratio'             in self.names: fracs[-1].gamma = self.data['SlipLength_ratio'][i]
                if 'InitialHydraulicAperture_m'   in self.names: fracs[-1].bh = self.data['InitialHydraulicAperture_m'][i]
                if 'InitialDilation_m'            in self.names: fracs[-1].bd = self.data['InitialDilation_m'][i]
                if 'Conductivity_m2m'             in self.names: fracs[-1].fc = self.data['Conductivity_m2m'][i]
                
                #other stuff
                fracs[-1].mu = np.tan(fracs[-1].phi)
        except:
            print('%s improperly formatted discrete fracture network datafile' %(filename))  
        #result
        return fracs

# ********************************************************************
### generate vtk file for fractures
# ********************************************************************
def vtk(faces,fname='default',scale=0.002*1500):
    #******   scaling       ******
    # r = 0.002*self.rock.size
    r = scale
        
    #******   paint flowing fractures   ******
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
    for i in range(0,len(faces)): #skip boundary node at np.inf
        # if faces[i].ci >= 0:
        #add colors
        q_0 += [i]
        q_1 += [faces[i].ci]
        q_2 += [faces[i].typ]
        q_3 += [faces[i].bd*1000]
        q_4 += [faces[i].bd0*1000]
        q_5 += [faces[i].bh*1000]
        q_6 += [faces[i].sn*1e-6]
        q_7 += [faces[i].Pcen]
        q_8 += [faces[i].Pc*1e-6]
        q_9 += [faces[i].stim]
        q_10 += [faces[i].Pmax]
        q_11 += [faces[i].tau*1e-6]
        if len(faces[i].Mws) > 0:
            q_12 += [np.max(faces[i].Mws)]
        else:
            q_12 += [-10.0]
        q_13 +=  [faces[i].bd0p*1000]
        q_14 +=  [faces[i].prop_load]
        q_15 +=  [3.28084*1.01324e15*(faces[i].bh**3.0)/12.0]
        #add geometry
        q_obj += [sg.HF(r=0.5*faces[i].dia, x0=faces[i].c0, strikeRad=faces[i].str, dipRad=faces[i].dip, h=0.02*r)]
    #vtk file
    q_col = [q_0,q_1,q_2,q_3,q_4,q_5,q_6,q_7,q_8,q_9,q_10,q_11,q_12,q_13,q_14,q_15]
    sg.writeVtk(q_obj, q_col, q_lab, vtkFile=(fname + '.vtk'))

# ********************************************************************
### boundary fractures
# ********************************************************************
def gen_domain(s):
    #working variables
    size = s.Domain_Size_m
    #create boundary faces for analysis
    domb3D = []
    domb3D += [surf([-size,0.0,0.0],4.0*size,00.0*deg,90.0*deg,'boundary',s)]
    domb3D += [surf([size,0.0,0.0],4.0*size,00.0*deg,90.0*deg,'boundary',s)]
    domb3D += [surf([0.0,-size,0.0],4.0*size,90.0*deg,90.0*deg,'boundary',s)]
    domb3D += [surf([0.0,size,0.0],4.0*size,90.0*deg,90.0*deg,'boundary',s)]
    domb3D += [surf([0.0,0.0,-size],4.0*size,00.0*deg,00.0*deg,'boundary',s)]
    domb3D += [surf([0.0,0.0,size],4.0*size,00.0*deg,00.0*deg,'boundary',s)]
    #add to model domain
    return domb3D        

# ********************************************************************
### hydraulic fracture creation
# ********************************************************************
def gen_hydrofracs(s,w):
    hydfs = []
    for i in range(0,len(w)):
        if w[i].type in ['procluster','perfcluster']:
            #hydraulic fractures
            num = s.Well_Clusters_count
            #location information
            spa = w[i].leg/(num+1)
            c0 = w[i].c0
            azn = w[i].azn
            dip = w[i].dip
            vAxi = w[i].vAxi
            #build list of fractures
            for n in range(0,num):
                #fracture parameters
                dia = 2.0*s.Stimulation_Radius_m * (1.0 + stats.norm_trunc(1,0.0,0.2,-0.5,0.5)[0]*s.Strategy_SecondOrderNoise_ratio)
                azn = np.pi/2.0+s.Stress_ShminAzimuth_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad*s.Strategy_SecondOrderNoise_ratio
                dip = np.pi/2.0-s.Stress_ShminDip_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad*s.Strategy_SecondOrderNoise_ratio
                c0 = c0 + spa*vAxi
                #add to model domain
                hydfs += [surf(c0,dia,azn,dip,'propped',s)]
                hydfs[-1].wi = w[i].pID
    #result
    return hydfs

# ********************************************************************
### maximum angle hydraulic fracture creation
# ********************************************************************
def gen_critfrac(s):
    critf = []
    #location information
    c0 = np.zeros(3)
    #normal
    dia = 2.0*s.Stimulation_Radius_m
    azn = np.pi/2.0+s.Stress_ShminAzimuth_rad + 0.0
    dip = np.pi/2.0-s.Stress_ShminDip_rad + 2.0*s.Stress_Variance_rad
    critf += [surf(c0,dia,azn,dip,'propped',s)]
    #strike-slip
    dia = 2.0*s.Stimulation_Radius_m
    azn = np.pi/2.0+s.Stress_ShminAzimuth_rad + 2.0*s.Stress_Variance_rad
    dip = np.pi/2.0-s.Stress_ShminDip_rad + 0.0
    critf += [surf(c0,dia,azn,dip,'propped',s)]
    #skew
    dia = 2.0*s.Stimulation_Radius_m
    azn = np.pi/2.0+s.Stress_ShminAzimuth_rad + 2.0*s.Stress_Variance_rad
    dip = np.pi/2.0-s.Stress_ShminDip_rad + 2.0*s.Stress_Variance_rad
    critf += [surf(c0,dia,azn,dip,'propped',s)]
    #find worst
    hi = np.argmax([critf[0].tau,critf[1].tau,critf[1].tau])
    critf = [critf[hi]]
    #result
    return critf
        
# ********************************************************************
### natural fracture network creation
# ********************************************************************
def gen_shearfracs(s,w):
    #size of domain for reservoir
    size = s.Domain_Size_m
    
    #joint set 1
    frac3D = []
    for n in range(0,s.JointSets_Set1_fractures):
        #fracture parameters
        logmu = 0.5*(np.log10(s.JointSets_DiameterMin_m)+np.log10(s.JointSets_DiameterMax_m))
        dia = stats.lognorm_trunc(1,logmu,logmu,np.log10(s.JointSets_DiameterMin_m),np.log10(s.JointSets_DiameterMax_m))[0]
        azn = s.JointSets_Set1Strike_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad
        dip = s.JointSets_Set1Dip_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad
        #Build geometry
        x = stats.uniform(1,-size,size)[0]
        y = stats.uniform(1,-size,size)[0]
        z = stats.uniform(1,-size,size)[0]
        c0 = np.asarray([x,y,z])
        #compile list of fractures
        frac3D += [surf(c0,dia,azn,dip,'fracture',s)]
        
    #joint set 2
    for n in range(0,s.JointSets_Set2_fractures):
        #fracture parameters
        logmu = 0.5*(np.log10(s.JointSets_DiameterMin_m)+np.log10(s.JointSets_DiameterMax_m))
        dia = stats.lognorm_trunc(1,logmu,logmu,np.log10(s.JointSets_DiameterMin_m),np.log10(s.JointSets_DiameterMax_m))[0]
        azn = s.JointSets_Set2Strike_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad
        dip = s.JointSets_Set2Dip_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad
        #Build geometry
        x = stats.uniform(1,-size,size)[0]
        y = stats.uniform(1,-size,size)[0]
        z = stats.uniform(1,-size,size)[0]
        c0 = np.asarray([x,y,z])
        #compile list of fractures
        frac3D += [surf(c0,dia,azn,dip,'fracture',s)]
        
    #joint set 3
    for n in range(0,s.JointSets_Set3_fractures):
        #fracture parameters
        logmu = 0.5*(np.log10(s.JointSets_DiameterMin_m)+np.log10(s.JointSets_DiameterMax_m))
        dia = stats.lognorm_trunc(1,logmu,logmu,np.log10(s.JointSets_DiameterMin_m),np.log10(s.JointSets_DiameterMax_m))[0]
        azn = s.JointSets_Set3Strike_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad
        dip = s.JointSets_Set3Dip_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad
        #Build geometry
        x = stats.uniform(1,-size,size)[0]
        y = stats.uniform(1,-size,size)[0]
        z = stats.uniform(1,-size,size)[0]
        c0 = np.asarray([x,y,z])
        #compile list of fractures
        frac3D += [surf(c0,dia,azn,dip,'fracture',s)]
        
    #add to model domain
    return frac3D
    
# ************************************************************************
### load fractures from file #TODO: remove cross referenceing to setup?
# ************************************************************************
def gen_from_dfn(s,filename=''):
    #get filepath if none given
    if filename == '':
        filename = s.Strategy_DiscreteFractures_path
    #if a filepath is provided, generate fractures
    if filename != '':
        dfn().intake(filename=filename,s=s)

# ************************************************************************
### seed hydraulic fractures
# ************************************************************************
def gen_seedfracs(s,w):
    hydfs = []
    for i in range(0,len(w)):
        if w[i].type in ['procluster','perfcluster']:
            #hydraulic fractures
            num = s.Well_Clusters_count
            #location information
            spa = w[i].leg/(num+1)
            c0 = w[i].c0
            azn = w[i].azn
            dip = w[i].dip
            vAxi = w[i].vAxi
            #build list of fractures
            for n in range(0,num):
                #fracture parameters
                dia = 2.0*s.Stimulation_SeedRadius_m * (1.0 + stats.norm_trunc(1,0.0,0.2,-0.5,0.5)[0]*s.Strategy_SecondOrderNoise_ratio)
                azn = np.pi/2.0+s.Stress_ShminAzimuth_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad*s.Strategy_SecondOrderNoise_ratio
                dip = np.pi/2.0-s.Stress_ShminDip_rad + stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*s.JointSets_Variance_rad*s.Strategy_SecondOrderNoise_ratio
                c0 = c0 + spa*vAxi
                #add to model domain
                hydfs += [surf(c0,dia,azn,dip,'propped',s)]
                hydfs[-1].wi = w[i].pID
    #result
    return hydfs
# ************************************************************************
### testing
# ************************************************************************
if False:
    #setup
    s=iogt.setup()
    
    #well geometry
    w = wells.gen_wells(s)
    wells.vtk(w,'test',0.005*s.Domain_Size_m)
    
    #fracture geometry
    f = gen_hydrofracs(s,w)
    vtk(f,'test_hf',0.005*s.Domain_Size_m)
    
    #natural fractures
    f = gen_shearfracs(s,w)
    vtk(f,'test_nf',0.005*s.Domain_Size_m)
    
    #critical fracture
    f = gen_critfrac(s)
    vtk(f,'test_critf',0.005*s.Domain_Size_m)