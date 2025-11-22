""" 
Automated well geometries:
    1. Break well into segments per the specified development strategy
    2. Determine depths, lengts, and diameters
    3. Place wells relative to eachother
    4. Output array of line objects that specify well information
"""
#initializations
import numpy as np
if __package__ is None or __package__ == '':
    import iogt
    import vtk as sg
    import csv
    from typ import typ
    from units import *
else:
    from . import iogt
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

#line objects
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
        self.typ = typ(self.type)
        
# ************************************************************************
# create vtk file for wells
# ************************************************************************
def vtk(wells,fname='default',scale=0.002*1500):
    #******   scaling       ******
    # r = 0.002*self.rock.size
    r = scale
    
    #******   paint wells   ******
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
    for i in range(0,len(wells)):
        #add colors
        w_0 += [i]
        w_1 += [wells[i].typ]
        w_2 += [wells[i].ra]
        w_3 += [wells[i].rgh]
        w_4 += [wells[i].rc]
        #add geometry
        azn = wells[i].azn
        dip = wells[i].dip
        leg = wells[i].leg
        vAxi = np.asarray([np.sin(azn)*np.cos(-dip), np.cos(azn)*np.cos(-dip), np.sin(-dip)])
        c0 = wells[i].c0
        c1 = c0 + vAxi*leg
        w_obj += [sg.cylObj(x0=c0, x1=c1, r=1.5*r)]
    #vtk file
    w_col = [w_0,w_1,w_2,w_3,w_4]
    sg.writeVtk(w_obj, w_col, w_lab, vtkFile=(fname + '.vtk'))

# ************************************************************************
# create vtk file for wells
# ************************************************************************
def vtk_zero(setup,wells,fname='default',scale=0.002*1500):
    #******   scaling       ******
    # r = 0.002*self.rock.size
    r = scale
    
    #******   paint wells   ******
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
    for i in range(0,len(wells)):
        #add colors
        w_0 += [i]
        w_1 += [wells[i].typ]
        w_2 += [wells[i].ra]
        w_3 += [wells[i].rgh]
        w_4 += [wells[i].rc]
        #add geometry
        azn = wells[i].azn
        dip = wells[i].dip
        leg = wells[i].leg
        vAxi = np.asarray([np.sin(azn)*np.cos(-dip), np.cos(azn)*np.cos(-dip), np.sin(-dip)])
        c0 = wells[i].c0 - np.asarray([0,0,setup.Well_TargetDepth_m])
        c1 = c0 + vAxi*leg 
        w_obj += [sg.cylObj(x0=c0, x1=c1, r=1.5*r)]
    #vtk file
    w_col = [w_0,w_1,w_2,w_3,w_4]
    sg.writeVtk(w_obj, w_col, w_lab, vtkFile=(fname + '.vtk'))

# ************************************************************************
# create csv file for wells
# ************************************************************************
#create a file containing all well information
def export(filename='well_exp.csv',w=[]):
    for i in range(0,len(w)):
        out = [['num',i]]
        for var in vars(w[i]):
            v = getattr(w[i],var)
            if np.isscalar(v):
                out += [[var,v]]
            elif var in ['c0','c1']:
                out += [['%s_E'%(var),v[0]]]
                out += [['%s_N'%(var),v[1]]]
                out += [['%s_Z'%(var),v[2]]]
        if i < 1:
            csv.save_csv(out=out,filename=filename,append=False)
        else:
            csv.save_csv(out=out,filename=filename,append=True)

# ************************************************************************
# read csv file for wells
# ************************************************************************
#create a file containing all well information
def intake(filename='well_exp.csv'):
    names, data = csv.load_csv(filename)
    w = []
    for i in range(0,len(data)):
        x0 = data['c0_E'][i]
        y0 = data['c0_N'][i]
        z0 = data['c0_Z'][i]
        leg = data['leg'][i]
        azn = data['azn'][i]
        dip = data['dip'][i]
        w_type = data['type'][i]
        ra = data['ra'][i]
        rb = data['rb'][i]
        rc = data['rc'][i]
        rgh = data['rgh'][i]
        pID = int(data['pID'][i]+0.001)
        w += [line(x0,y0,z0,leg,azn,dip,w_type,ra,rb,rc,rgh,pID)]
    return w

# ************************************************************************
# well placement
# ************************************************************************
def gen_wells(s=iogt.setup()):
    #print( '*** well placement module ***')

    #initialization
    wells = []
    n = len(wells)
    #centerpoint, target depth (z0) = 0.0
    i0 = np.asarray([0.0, 0.0, 0.0])
    #ref axes (parallel injector and 90 horizontal to the right)
    azn = s.Well_Azimuth_rad
    dip = s.Well_Dip_rad
    vInj = np.asarray([np.sin(azn)*np.cos(-dip),
             np.cos(azn)*np.cos(-dip),np.sin(-dip)])
    vRht = np.asarray([np.sin(azn+np.pi/2),np.cos(azn+np.pi/2),0.0])
    # vNor = np.asarray([np.sin(azn)*np.sin(-dip),
    #          np.cos(azn)*np.sin(-dip),np.cos(-dip)])
    vNor = np.cross(vRht,vInj)
    vNor = vNor/np.linalg.norm(vNor)
    #offset center
    p0 = i0 + vRht*s.Well_Spacing_m
    #toe in production well
    toe = s.Well_RotationToe_rad
    vPro = sg.rotatePoints([vInj],vNor,toe)[0]
    #skew in production well
    skew = s.Well_RotationSkew_rad
    vPro = sg.rotatePoints([vPro],vRht,skew)[0]
    #lengths of wells
    length = s.Well_DeviatedLength_m
    proportion = s.Well_ProducerProportion_ratio
    iLen = length*proportion
    pLen = length
    cLen = pLen - iLen
    #segments
    depth = s.Well_TargetDepth_m
    segle = depth/s.Well_SegmentCount_units #segment length, m
    num = s.Well_ProducerCount_wells
    phase = -1.0*s.Well_RotationPhase_rad
    seg = s.Well_Intervals_count
    leg = iLen/(seg+seg-1)
    rgh = s.Well_Roughness_metric
    #conductor segment
    h = np.asarray([0.0,0.0,-segle])
    #radii
    pra = s.Well_ProductionRadius_m - s.Well_WallThickness_m
    prb = s.Well_ProductionRadius_m 
    prc = s.Well_ProductionRadius_m + s.Well_AnnularThickness_m
    ira = prc
    irb = ira + s.Well_WallThickness_m
    irc = irb + s.Well_AnnularThickness_m
    sra = irc
    srb = sra + s.Well_WallThickness_m
    src = srb + s.Well_AnnularThickness_m
    
    # ************************************************************************
    # O&G - wine-rack style wells with casing to bottom in all wells
    # ************************************************************************
    if s.Strategy_Design_type.upper() in ['O&G','DCM']:
        #make wells parallel
        vPro = vInj
        #populate rack of production and injection wells
        i0s = []
        p0s = []
        for i in range(0,num):
            if i == 0: 
                #first production well includes a mirrored well
                p0s += [sg.rotatePoints([p0],vInj,(np.pi - phase))[0]]
                i0s += [i0]
            else: 
                #place additional producers
                p0s += [p0s[-1] + 2.0*vRht*s.Well_Spacing_m*np.cos(phase)]
            #don't place injector in second step
            if i > 1:
                i0s += [i0s[-1] + 2.0*vRht*s.Well_Spacing_m*np.cos(phase)]
        #recenter
        c0 = np.zeros(3)
        if num < 2:
            c0 = c0 - 0.5*vRht*s.Well_Spacing_m*np.cos(phase)
        else:
            c0 = c0 + (num-2)*vRht*s.Well_Spacing_m*np.cos(phase)
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
            azn, dip, dis = sg.azn_dip(iss,iss+h)
            wells += [line(iss[0],iss[1],iss[2],dis,azn,dip,'injector',sra,srb,src,rgh,iid[i])]
        for i in range(0,len(p0s)):
            azn, dip, dis = sg.azn_dip(pss[i],pss[i]+h)
            wells += [line(pss[i][0],pss[i][1],pss[i][2],dis,azn,dip,'producer',sra,srb,src,rgh,pid[i])]
        #place surface segments
        for i in range(0,len(i0s)):
            iss = isss[i]
            ibs = ibss[i]
            azn, dip, dis = sg.azn_dip(iss+h,ibs)
            wells += [line(iss[0]+h[0],iss[1]+h[1],iss[2]+h[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,iid[i])]
        for i in range(0,len(p0s)):
            azn, dip, dis = sg.azn_dip(pss[i]+h,p1s[i])
            wells += [line(pss[i][0]+h[0],pss[i][1]+h[1],pss[i][2]+h[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,pid[i])]
        #place injector base segment
        for i in range(0,len(i0s)):
            ibs = ibss[i]
            i1s = i1ss[i]
            i2s = i2ss[i]
            azn, dip, dis = sg.azn_dip(ibs,i1s[0])
            wells += [line(ibs[0],ibs[1],ibs[2],dis,azn,dip,'pipe',pra,prb,prc,rgh,iid[i])]
            #place injection internal segments
            for i in range(0,seg-1):
                azn, dip, dis = sg.azn_dip(i1s[0],i2s[0])
                wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],dis,azn,dip,'perfcluster',pra,prb,prc,rgh,iid[i])]
                azn, dip, dis2 = sg.azn_dip(i2s[0],i1s[1])
                wells += [line(i2s[i][0],i2s[i][1],i2s[i][2],dis2,azn,dip,'pipe',pra,prb,prc,rgh,iid[i])]
            azn, dip, dis = sg.azn_dip(i1s[-1],i2s[-1])
            wells += [line(i1s[-1][0],i1s[-1][1],i1s[-1][2],dis,azn,dip,'perfcluster',pra,prb,prc,rgh,iid[i])]
        #place production wells
        for i in range(0,len(p0s)):
            azn, dip, dis = sg.azn_dip(p1s[i],p2s[i])
            wells += [line(p1s[i][0],p1s[i][1],p1s[i][2],dis,azn,dip,'screen',pra,prb,prc,rgh,pid[i])]
        #stimulate all wells for O&G
        if s.Strategy_Design_type.upper() == 'O&G':
            for w in wells:
                if typ(w.typ) == 'screen':
                    w.type = 'procluster'
                    w.typ = typ('procluster')
            
    # ************************************************************************
    # AGS - closed loop wells with casing
    # ************************************************************************
    elif s.Strategy_Design_type == 'AGS': #closed loop geothermal
        separ = s.Well_Spacing_m #well surface separation, m
        branc = int(s.Well_ProducerCount_wells/2) #single-sided number of branches (mirrored about center well)
        mod = s.Well_ProducerCount_wells%2
        #input overrides for economics
        s.Well_ProducerProportion_ratio = 1.0
        #general geometry
        vAxi = vInj
        vert_z0 = 0.5*np.sin(dip)*length
        hori_r0 = length*np.cos(dip)
        offs_z0 = s.Well_Spacing_m/np.cos(dip)
        #place injector and producer wells as the first two segments
        inje_c0 = np.asarray([-0.5*hori_r0*np.sin(azn),-0.5*hori_r0*np.cos(azn),depth])
        prod_c0 = np.asarray([-(0.5*hori_r0-separ)*np.sin(azn),-(0.5*hori_r0-separ)*np.cos(azn),depth])
        wells += [line(inje_c0[0],inje_c0[1],inje_c0[2],segle,0.0*deg,90.0*deg,'injector',sra,srb,src,rgh,0+n)]
        wells += [line(prod_c0[0],prod_c0[1],prod_c0[2],segle,0.0*deg,90.0*deg,'producer',sra,srb,src,rgh,1+n)]
        #segment information: i = injector; v = vertical; l = length; n = count
        ivl = depth - vert_z0 - segle
        ivn = int(ivl/segle) + 1
        ihl = length
        ihn = int(ihl/segle) + 1
        pvl = depth - vert_z0 - offs_z0 + separ*np.tan(dip) - segle
        pvn = int(pvl/segle) + 1
        phl = length - separ/np.cos(dip)
        phn = int(phl/segle) + 1
        #place injector segments
        oi = inje_c0 + np.asarray([0,0,-segle])
        ol = np.zeros(3,dtype=float)
        for b in range(0,ivn):
            wells += [line(oi[0],oi[1],oi[2], ivl/ivn, 0.0*deg,90.0*deg, 'pipe',ira,irb,irc,rgh,0+n)]
            oi = oi + np.asarray([0.0,0.0,-ivl/ivn])
        ol = ol + oi
        for b in range(0,ihn):
            wells += [line(oi[0],oi[1],oi[2], ihl/ihn, azn, dip, 'pipe',pra,prb,prc,rgh,0+n)]
            oi = oi + vAxi*(ihl/ihn)
        #place producer segments
        op = prod_c0 + np.asarray([0,0,-segle])
        ou = np.zeros(3,dtype=float)
        for b in range(0,pvn):
            wells += [line(op[0],op[1],op[2], pvl/pvn, 0.0*deg,90.0*deg, 'pipe',ira,irb,irc,rgh,1+n)]
            op = op + np.asarray([0.0,0.0,-pvl/pvn])
        ou = ou + op
        for b in range(0,phn):
            wells += [line(op[0],op[1],op[2], phl/phn, azn, dip, 'pipe',pra,prb,prc,rgh,1+n)]
            op = op + vAxi*(phl/phn)
        #place bottom hole connector
        wells += [line(oi[0],oi[1],oi[2], offs_z0, 0.0*deg, -90.0*deg, 'pipe',pra,prb,prc,rgh,1+n)]
        #transverse vector
        vNor = np.cross(vAxi,np.asarray([0,0,1]))
        vNor = vNor/np.linalg.norm(vNor)
        #add branches
        for v in range(0,branc):
            #even branch
            wells += [line(ou[0],ou[1],ou[2], separ*(v+1), azn+90.0*deg,0.0*deg,'pipe',pra,prb,prc,rgh,1+n)]
            wells += [line(ol[0],ol[1],ol[2], separ*(v+1), azn+90.0*deg,0.0*deg,'pipe',pra,prb,prc,rgh,0+n)]
            #place injector segments
            oi = np.copy(ol) + vNor*separ*(v+1)
            for b in range(0,ihn):
                wells += [line(oi[0],oi[1],oi[2], ihl/ihn, azn, dip, 'pipe',pra,prb,prc,rgh,0+n)]
                oi = oi + vAxi*(ihl/ihn)
            #place producer segments
            op = np.copy(ou) + vNor*separ*(v+1)
            for b in range(0,phn):
                wells += [line(op[0],op[1],op[2], phl/phn, azn, dip, 'pipe',pra,prb,prc,rgh,1+n)]
                op = op + vAxi*(phl/phn)
            #place bottom connectors
            wells += [line(oi[0],oi[1],oi[2], offs_z0, 0.0*deg, -90.0*deg, 'pipe',pra,prb,prc,rgh,1+n)]
            #odd branch
            if mod>0:
                wells += [line(ou[0],ou[1],ou[2], separ*(v+1), azn-90.0*deg,0.0*deg,'pipe',pra,prb,prc,rgh,1+n)]
                wells += [line(ol[0],ol[1],ol[2], separ*(v+1), azn-90.0*deg,0.0*deg,'pipe',pra,prb,prc,rgh,0+n)]
                #place injector segments
                oi = np.copy(ol) - vNor*separ*(v+1)
                for b in range(0,ihn):
                    wells += [line(oi[0],oi[1],oi[2], ihl/ihn, azn, dip, 'pipe',pra,prb,prc,rgh,0+n)]
                    oi = oi + vAxi*(ihl/ihn)
                #place producer segments
                op = np.copy(ou) - vNor*separ*(v+1)
                for b in range(0,phn):
                    wells += [line(op[0],op[1],op[2], phl/phn, azn, dip, 'pipe',pra,prb,prc,rgh,1+n)]
                    op = op + vAxi*(phl/phn)  
                #place bottom connectors
                wells += [line(oi[0],oi[1],oi[2], offs_z0, 0.0*deg, -90.0*deg, 'pipe',pra,prb,prc,rgh,1+n)]
        
    # ************************************************************************
    # FGS - zonally isolated injection well with uncased producers #TODO update module to account for well friction to surface in the injector and to have surface at correct elevation
    # ************************************************************************
    elif (s.Strategy_Design_type in ['FGS','KGS']):

        #phase of production wells
        p0s = []
        vPros = []
        for i in range(0,num):
            p0s += [sg.rotatePoints([p0],vInj,(i*2.0*np.pi/num + phase))[0]]
            vPros += [sg.rotatePoints([vPro],vInj,(i*2.0*np.pi/num + phase))[0]]

        #injection well segments
        i1s = []
        i2s = []
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
        oDepth = depth + i2s[-1][2]
        #place injection wells
        wells = []
        azn, dip, dis = sg.azn_dip(i1s[0],i2s[0])
        vAxi = np.asarray([np.sin(azn)*np.cos(-dip),np.cos(azn)*np.cos(-dip),np.sin(-dip)])
        for i in range(0,seg):
            wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],0.1*leg,azn,dip,'injector',sra,srb,src,rgh,i+n)]
        #place production well surface segments
        for j in range(0,num):
            #extra = 0.05*depth
            wells += [line(p1s[j][0],p1s[j][1],oDepth,segle,0.0,90*deg,'producer',sra,srb,src,rgh,1+i+j+n)]
        #place injection well lower segments
        wells += [line(i1s[0][0],i1s[0][1],oDepth,oDepth-i1s[0][2],0.0,90.0*deg,'pipe',ira,irb,irc,rgh,0+n)]
        #place perforation clusters
        for i in range(0,seg):
            i1s[i] = i1s[i] + vAxi*0.1*leg
            wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],0.9*leg,azn,dip,'perfcluster',pra,prb,prc,rgh,i+n)]
        #place production well lower segments
        for j in range(0,num):
            #vertical
            wells += [line(p1s[j][0],p1s[j][1],oDepth+h[2],oDepth+h[2]-p1s[j][2],0.0,90.0*deg,'pipe',ira,irb,irc,rgh,1+i+j+n)]
            #inclined
            azn, dip, dis = sg.azn_dip(p1s[j],p2s[j])
            vAxi = np.asarray([np.sin(azn)*np.cos(-dip),np.cos(azn)*np.cos(-dip),np.sin(-dip)])
            wells += [line(p1s[j][0],p1s[j][1],p1s[j][2],pLen,azn,dip,'screen',pra,prb,prc,rgh,1+i+j+n)]
            
    # ************************************************************************
    # CGS - conventional wells
    # ************************************************************************
    elif (s.Strategy_Design_type in ['CGS']):
        #make wells parallel
        vPro = vInj
        #populate rack of production and injection wells
        i0s = []
        p0s = []
        lRht = 0.5*s.Well_Spacing_m
        lNor = s.Well_Spacing_m
        for i in range(0,num):
            #place production wells
            p0s += [i0 + vRht*s.Well_Spacing_m*i]
            if (s.Well_InjectorCount_wells > 0):
                if i == 0:
                    i0s += [i0 + vRht*lRht + vNor*lNor]
                elif i > 1:
                    i0s += [i0s[-1] + vRht*s.Well_Spacing_m]
        #recenter
        c0 = np.mean(np.asarray(p0s),axis=0)
        for i in range(0,len(i0s)):
            i0s[i] += -1.0*c0
        for i in range(0,len(p0s)):
            p0s[i] += -1.0*c0  
        #injection well segments
        i1s = [] #startpoints
        i2s = [] #endpoints
        iss = [] #surface segments
        for i in range(0,len(i0s)):
            i1s += [i0s[i] - 0.5*vInj*iLen]
            i2s += [i0s[i] + 0.5*vInj*iLen]    
            iss += [np.asarray([i1s[i][0],i1s[i][1],depth+0.1])]
        #production well segments
        p1s = [] #startpoints
        p2s = [] #endpoints
        pss = [] #surface segments
        for i in range(0,len(p0s)):
            p1s += [p0s[i] - 0.5*vPro*pLen]
            p2s += [p0s[i] + 0.5*vPro*pLen]    
            pss += [np.asarray([p1s[i][0],p1s[i][1],depth+0.1])]
        #place wellheads
        iid = np.asarray(range(0,len(i0s)))+n
        pid = np.asarray(range(0,len(p0s)))+n+len(i0s)
        for i in range(0,len(i0s)):
            azn, dip, dis = sg.azn_dip(iss[i],iss[i]+h)
            wells += [line(iss[i][0],iss[i][1],iss[i][2],dis,azn,dip,'injector',sra,srb,src,rgh,iid[i])]
        for i in range(0,len(p0s)):
            azn, dip, dis = sg.azn_dip(pss[i],pss[i]+h)
            wells += [line(pss[i][0],pss[i][1],pss[i][2],dis,azn,dip,'producer',sra,srb,src,rgh,pid[i])]
        #place surface segments
        for i in range(0,len(i0s)):
            azn, dip, dis = sg.azn_dip(iss[i]+h,i1s[i])
            wells += [line(iss[i][0]+h[0],iss[i][1]+h[1],iss[i][2]+h[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,iid[i])]
        for i in range(0,len(p0s)):
            azn, dip, dis = sg.azn_dip(pss[i]+h,p1s[i])
            wells += [line(pss[i][0]+h[0],pss[i][1]+h[1],pss[i][2]+h[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,pid[i])]
        #place injection wells
        for i in range(0,len(i0s)):
            azn, dip, dis = sg.azn_dip(i1s[i],i2s[i])
            wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],dis,azn,dip,'screen',pra,prb,prc,rgh,iid[i])]
        #place production wells
        for i in range(0,len(p0s)):
            azn, dip, dis = sg.azn_dip(p1s[i],p2s[i])
            wells += [line(p1s[i][0],p1s[i][1],p1s[i][2],dis,azn,dip,'screen',pra,prb,prc,rgh,pid[i])]
            
    # ************************************************************************
    # EGS - all wells cased to bottom centered on one injection well
    # ************************************************************************
    else: # (s.Strategy_Design_type in ['EGS','DCM']):
        #phase of production wells
        p0s = []
        vPros = []
        for i in range(0,num):
            p0s += [sg.rotatePoints([p0],vInj,(i*2.0*np.pi/num + phase))[0]]
            vPros += [sg.rotatePoints([vPro],vInj,(i*2.0*np.pi/num + phase))[0]]

        #injection well reservoir segments
        i1s = [] #startpoints
        i2s = [] #endpoints
        i1s += [i0 - 0.5*vInj*iLen]
        i2s += [i1s[0] + leg*vInj]
        for i in range(1,seg):
            i1s += [i2s[i-1] + leg*vInj]
            i2s += [i1s[i] + leg*vInj]
        #injection well segments to surface
        ibs = i1s[0] - 0.5*cLen*vInj
        iss = np.asarray([ibs[0],ibs[1],depth])
        #production well segments
        p1s = [] #startpoints
        p2s = [] #endpoints
        for i in range(0,num):
            p1s += [p0s[i] - 0.5*vPros[i]*pLen]
            p2s += [p0s[i] + 0.5*vPros[i]*pLen]    
        #production well segments to surface
        pss = []
        for i in range(0,num):
            pss += [np.asarray([p1s[i][0],p1s[i][1],depth+0.1])]
        #wellheads
        azn, dip, dis = sg.azn_dip(iss,iss+h)
        wells += [line(iss[0],iss[1],iss[2],dis,azn,dip,'injector',sra,srb,src,rgh,0+n)]
        for i in range(0,num):
            azn, dip, dis = sg.azn_dip(pss[i],pss[i]+h)
            wells += [line(pss[i][0],pss[i][1],pss[i][2],dis,azn,dip,'producer',sra,srb,src,rgh,1+i+n)]
        #place surface segments
        azn, dip, dis = sg.azn_dip(iss+h,ibs)
        wells += [line(iss[0]+h[0],iss[1]+h[1],iss[2]+h[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,0+n)]
        for i in range(0,num):
            azn, dip, dis = sg.azn_dip(pss[i]+h,p1s[i])
            wells += [line(pss[i][0]+h[0],pss[i][1]+h[1],pss[i][2]+h[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,1+i+n)]
        #place injector base segment
        azn, dip, dis = sg.azn_dip(ibs,i1s[0])
        wells += [line(ibs[0],ibs[1],ibs[2],dis,azn,dip,'pipe',pra,prb,prc,rgh,0+n)]
        #place injection wells
        for i in range(0,seg-1):
            azn, dip, dis = sg.azn_dip(i1s[0],i2s[0])
            wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],dis,azn,dip,'perfcluster',pra,prb,prc,rgh,0+n)]
            azn, dip, dis2 = sg.azn_dip(i2s[0],i1s[1])
            wells += [line(i2s[i][0],i2s[i][1],i2s[i][2],dis2,azn,dip,'pipe',pra,prb,prc,rgh,0+n)]
        azn, dip, dis = sg.azn_dip(i1s[-1],i2s[-1])
        wells += [line(i1s[-1][0],i1s[-1][1],i1s[-1][2],dis,azn,dip,'perfcluster',pra,prb,prc,rgh,0+n)]
        #place production wells
        for i in range(0,num):
            azn, dip, dis = sg.azn_dip(p1s[i],p2s[i])
            wells += [line(p1s[i][0],p1s[i][1],p1s[i][2],dis,azn,dip,'screen',pra,prb,prc,rgh,1+i+n)]
            
    # ************************************************************************
    ### offset geometry to zero-surface positioning
    # ************************************************************************
    leg = s.Well_DeviatedLength_m
    azn = s.Well_Azimuth_rad
    dip = s.Well_Dip_rad
    depth = s.Well_TargetDepth_m
    offset = np.asarray([0.0, 0.0, -depth])
    if s.Well_Pattern_count > 1:
        offset += 0.5*(leg+20)*np.asarray([np.sin(azn)*np.cos(-dip), np.cos(azn)*np.cos(-dip), 0.0])
    print(offset)
    for w in wells:
        w.c0 += offset
        w.c1 += offset
    
    # ************************************************************************
    ### return
    # ************************************************************************
    return wells

#testing
if False:
    x = iogt.setup()
    x.Well_DeviatedLength_m = 1500
    x.Well_Dip_rad = 0
    
    x.Strategy_Design_type = 'EGS'
    x.Well_ProducerCount_wells = 3
    w = gen_wells(x)
    vtk_zero(x,w,'EGS3',0.005*x.Domain_Size_m)
    
    x.Strategy_Design_type = 'EGS'
    x.Well_ProducerCount_wells = 1
    x.Well_RotationPhase_rad = np.pi/2
    x.re_init()
    w = gen_wells(x)
    vtk_zero(x,w,'EGS2',0.005*x.Domain_Size_m)
    
    x.Strategy_Design_type = 'AGS'
    x.Well_ProducerCount_wells = 2
    x.Well_RotationPhase_rad = 0.0
    x.re_init()
    w = gen_wells(x)
    vtk_zero(x,w,'AGS2',0.005*x.Domain_Size_m)
    
    x.Strategy_Design_type = 'AGS'
    x.Well_ProducerCount_wells = 11
    x.Well_Dip_rad = 45*deg
    x.Well_Spacing_m = 75
    x.re_init()
    w = gen_wells(x)
    vtk_zero(x,w,'AGS11',0.005*x.Domain_Size_m)
    x.Well_Dip_rad = 0.0*deg
    x.Well_Spacing_m = 175
    
    x.Strategy_Design_type = 'O&G'
    x.Well_ProducerCount_wells = 2
    x.Well_RotationPhase_rad = np.pi/8
    x.re_init()
    w = gen_wells(x)
    vtk_zero(x,w,'O&G2',0.005*x.Domain_Size_m)
    
    x.Strategy_Design_type = 'O&G'
    x.Well_ProducerCount_wells = 4
    x.Well_RotationPhase_rad = np.pi/8
    x.re_init()
    w = gen_wells(x)
    vtk_zero(x,w,'O&G4',0.005*x.Domain_Size_m)

    #upscaling chunks
    for i in w:
        i.c0[0] += -50
        i.c0[2] += -400
        i.c1[0] == -50
        i.c1[2] += -400
    vtk_zero(x,w,'O&G8',0.005*x.Domain_Size_m)
    for i in w:
        i.c0[0] += -50
        i.c0[2] += 400
        i.c1[0] == -50
        i.c1[2] += 400
        i.azn = -i.azn
        i.c1[1] = -i.c1[1]
    vtk_zero(x,w,'O&G12',0.005*x.Domain_Size_m)
    for i in w:
        i.c0[0] += -50
        i.c0[2] += -400
        i.c1[0] == -50
        i.c1[2] += -400
    vtk_zero(x,w,'O&G16',0.005*x.Domain_Size_m)
    
    x.Strategy_Design_type = 'CGS'
    x.Well_ProducerCount_wells = 1
    x.Well_Pattern_count = 0
    x.Well_RotationPhase_rad = np.pi/6
    x.re_init()
    w = gen_wells(x)
    vtk_zero(x,w,'CGS0',0.005*x.Domain_Size_m)
    
    x.Strategy_Design_type = 'CGS'
    x.Well_ProducerCount_wells = 3
    x.Well_Pattern_count = 1
    x.Well_RotationPhase_rad = np.pi/2.2
    sp = x.Well_Spacing_m
    x.Well_Spacing_m = 2000
    x.re_init()
    w = gen_wells(x)
    vtk_zero(x,w,'CGS3',0.005*x.Domain_Size_m)
    
    x.Strategy_Design_type = 'FGS'
    x.Well_ProducerCount_wells = 4
    x.Well_Pattern_count = 1
    x.Well_RotationPhase_rad = 0.0
    x.Well_Spacing_m = sp
    x.re_init()
    w = gen_wells(x)
    vtk_zero(x,w,'FGS4',0.005*x.Domain_Size_m)