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
    import io
    import vtk as sg
    from units import *
else:
    from . import io
    from . import vtk as sg
    from .units import *
deg=deg
yr=yr
cP=cP
darcy=darcy
mD=mD
g=g

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

# ************************************************************************
# well placement #TODO: add ability to remove injection wells
# ************************************************************************
def gen_wells(s=io.setup()):
    print( '*** well placement module ***')
    
    
    
    #dummy
    s.Well_TargetTemp_K = 273.15+325 #K
    s.Well_TargetDepth_m = 6000.0 #m
    s.Well_ProducerCount_wells = 4 #wells
    s.Well_InjectorCount_ratio = 1 #ratio of injection wells to production wells #TODO: needs to be implemented
    s.Well_InjectorCount_wells = 1
    s.Well_Spacing_m = 400.0 #m
    s.Well_DeviatedLength_m = 800.0 #m
    s.Well_ProducerProportion_ratio = 0.6 #m/m
    s.Well_Azimuth_rad = 90.0*deg #rad
    s.Well_Dip_rad = 30.0*deg #rad
    s.Well_RotationPhase_rad = 0.0*deg #rad
    s.Well_RotationToe_rad = 0.0*deg #rad
    s.Well_RotationSkew_rad = 0.0*deg #rad
    s.Well_ScaleFactor_ratio = 1.0 #m TODO: Tool to scale diameter of well
    s.Well_Roughness_metric = 80.0
    
    s.Well_Intervals_count = 1
    s.Well_Clusters_count = 1
    
    s.Well_ProductionRadius_m = 1
    s.Well_IntermediateRadius_m = 1
    s.Well_SurfaceRadius_m = 1
    s.Well_WallThickness_m = 0.0254*0.75 #m #TODO: Implement complete well plan
    s.Well_AnnularThickness_m = 0.0254*1.25 #m #TODO: Implement complete well plan

    
    
    

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
    vNor = np.asarray([np.sin(azn)*np.sin(-dip),
             np.cos(azn)*np.sin(-dip),np.cos(-dip)])
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
    num = s.Well_ProducerCount_wells
    phase = s.Well_RotationPhase_rad
    seg = s.Well_Intervals_count
    leg = iLen/(seg+seg-1)
    rgh = s.Well_Roughness_metric
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
    if (s.Strategy_Design_type == 'O&G'):
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
            iss = np.asarray([ibs[0],ibs[1],s.Well_TargetDepth_m])
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
            pss += [np.asarray([p1s[i][0],p1s[i][1],s.Well_TargetDepth_m+0.1])]
        #place wellheads
        h = np.asarray([0.0,0.0,s.Well_Spacing_m])
        iid = np.asarray(range(0,len(i0s)))+n
        pid = np.asarray(range(0,len(p0s)))+n+len(i0s)
        for i in range(0,len(i0s)):
            iss = isss[i]
            azn, dip, dis = sg.azn_dip(iss+h,iss)
            wells += [line(iss[0]+h[0],iss[1]+h[1],iss[2]+h[2],dis,azn,dip,'injector',sra,srb,src,rgh,iid[i])]
        for i in range(0,len(p0s)):
            azn, dip, dis = sg.azn_dip(pss[i]+h,pss[i])
            wells += [line(pss[i][0]+h[0],pss[i][1]+h[1],pss[i][2]+h[2],dis,azn,dip,'producer',sra,srb,src,rgh,pid[i])]
        #place surface segments
        for i in range(0,len(i0s)):
            iss = isss[i]
            ibs = ibss[i]
            azn, dip, dis = sg.azn_dip(iss,ibs)
            wells += [line(iss[0],iss[1],iss[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,iid[i])]
        for i in range(0,len(p0s)):
            azn, dip, dis = sg.azn_dip(pss[i],p1s[i])
            wells += [line(pss[i][0],pss[i][1],pss[i][2],dis,azn,dip,'pipe',ira,irb,irc,rgh,pid[i])]
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
            
    # ************************************************************************
    # AGS - closed loop wells with casing
    # ************************************************************************
    elif s.Strategy_Design_type == 'AGS': #closed loop geothermal
        segle = s.Well_TargetDepth_m/20.0 #segment length, m
        separ = s.Well_Spacing_m #well surface separation, m
        branc = int(s.Well_ProducerCount_wells/2) #single-sided number of branches (mirrored about center well) 
        #input overrides for economics
        s.Well_ProducerProportion_ratio = 1.0
        #general geometry
        vAxi = vInj
        vert_z0 = np.sin(dip)*length
        orig_d0 = s.Well_TargetDepth_m-np.sin(dip)*0.5*length
        hori_r0 = length*np.cos(dip)
        offs_z0 = s.Well_Spacing_m/np.cos(dip)
        #place injector and producer wells as the first two segments
        inje_c0 = np.asarray([-0.5*hori_r0*np.sin(azn),-0.5*hori_r0*np.cos(azn),orig_d0])
        prod_c0 = np.asarray([-(0.5*hori_r0-separ)*np.sin(azn),-(0.5*hori_r0-separ)*np.cos(azn),orig_d0])
        wells += [line(inje_c0[0],inje_c0[1],inje_c0[2],segle,0.0*deg,90.0*deg,'injector',sra,srb,src,rgh,0+n)]
        wells += [line(prod_c0[0],prod_c0[1],prod_c0[2],segle,0.0*deg,90.0*deg,'producer',sra,srb,src,rgh,1+n)]
        #segment information: i = injector; v = vertical; l = length; n = count
        ivl = s.Well_TargetDepth_m - vert_z0 - segle
        ivn = int(ivl/segle) + 1
        ihl = length
        ihn = int(ihl/segle) + 1
        pvl = s.Well_TargetDepth_m - vert_z0 - offs_z0 + separ*np.tan(dip) - segle
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
            #add transverse elements
            wells += [line(ou[0],ou[1],ou[2], separ*(v+1), azn+90.0*deg,0.0*deg,'pipe',pra,prb,prc,rgh,1+n)]
            wells += [line(ou[0],ou[1],ou[2], separ*(v+1), azn-90.0*deg,0.0*deg,'pipe',pra,prb,prc,rgh,1+n)]
            wells += [line(ol[0],ol[1],ol[2], separ*(v+1), azn+90.0*deg,0.0*deg,'pipe',pra,prb,prc,rgh,0+n)]
            wells += [line(ol[0],ol[1],ol[2], separ*(v+1), azn-90.0*deg,0.0*deg,'pipe',pra,prb,prc,rgh,0+n)]
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
    # FGS - zonally isolated injection well with uncased producers #TODO remove casing
    # ************************************************************************
    elif (s.Strategy_Design_type in ['FGS']): #patterned geothermal systems

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
        oDepth = s.Well_TargetDepth_m + i2s[-1][2]            
        #place injection wells
        wells = []
        azn, dip, dis = sg.azn_dip(i1s[0],i2s[0])
        vAxi = np.asarray([np.sin(azn)*np.cos(-dip),np.cos(azn)*np.cos(-dip),np.sin(-dip)])
        for i in range(0,seg):
            wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],0.1*leg,azn,dip,'injector',sra,srb,src,rgh,i+n)]
        #place production well surface segments
        for j in range(0,num):
            extra = 0.05*s.Well_TargetDepth_m
            wells += [line(p1s[j][0],p1s[j][1],oDepth+extra,extra,0.0,90*deg,'producer',sra,srb,src,rgh,1+i+j+n)]
        #place injection well lower segments
        wells += [line(i1s[0][0],i1s[0][1],oDepth,oDepth-i1s[0][2],0.0,90.0*deg,'pipe',ira,irb,irc,rgh,0+n)]
        #place perforation clusters
        for i in range(0,seg):
            i1s[i] = i1s[i] + vAxi*0.1*leg
            wells += [line(i1s[i][0],i1s[i][1],i1s[i][2],0.9*leg,azn,dip,'perfcluster',pra,prb,prc,rgh,i+n)]
        #place production well lower segments
        for j in range(0,num):
            #vertical
            wells += [line(p1s[j][0],p1s[j][1],oDepth,oDepth-p1s[j][2],0.0,90.0*deg,'pipe',ira,irb,irc,rgh,1+i+j+n)]
            #inclined
            azn, dip, dis = sg.azn_dip(p1s[j],p2s[j])
            vAxi = np.asarray([np.sin(azn)*np.cos(-dip),np.cos(azn)*np.cos(-dip),np.sin(-dip)])
            wells += [line(p1s[j][0],p1s[j][1],p1s[j][2],pLen,azn,dip,'screen',pra,prb,prc,rgh,1+i+j+n)]
            
    # ************************************************************************
    # CGS - conventional wells
    # ************************************************************************
    elif (s.Strategy_Design_type in ['CGS']):
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
        iss = np.asarray([ibs[0],ibs[1],s.Well_TargetDepth_m])
        #production well segments
        p1s = [] #startpoints
        p2s = [] #endpoints
        for i in range(0,num):
            p1s += [p0s[i] - 0.5*vPros[i]*pLen]
            p2s += [p0s[i] + 0.5*vPros[i]*pLen]    
        #production well segments to surface
        pss = []
        for i in range(0,num):
            pss += [np.asarray([p1s[i][0],p1s[i][1],s.Well_TargetDepth_m+0.1])]
        #wellheads
        h = np.asarray([0.0,0.0,s.Well_Spacing_m])
        azn, dip, dis = sg.azn_dip(iss+h,iss)
        wells += [line(iss[0]+h[0],iss[1]+h[1],iss[2]+h[2],dis,azn,dip,'injector',sra,srb,src,rgh,0+n)]
        for i in range(0,num):
            azn, dip, dis = sg.azn_dip(pss[i]+h,pss[i])
            wells += [line(pss[i][0]+h[0],pss[i][1]+h[1],pss[i][2]+h[2],dis,azn,dip,'producer',sra,srb,src,rgh,1+i+n)]
        #place surface segments
        azn, dip, dis = sg.azn_dip(iss,ibs)
        wells += [line(iss[0],iss[1],iss[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,0+n)]
        for i in range(0,num):
            azn, dip, dis = sg.azn_dip(pss[i],p1s[i])
            wells += [line(pss[i][0],pss[i][1],pss[i][2],dis,azn,dip,'pipe',ira,irb,irc,rgh,1+i+n)]
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
        iss = np.asarray([ibs[0],ibs[1],s.Well_TargetDepth_m])
        #production well segments
        p1s = [] #startpoints
        p2s = [] #endpoints
        for i in range(0,num):
            p1s += [p0s[i] - 0.5*vPros[i]*pLen]
            p2s += [p0s[i] + 0.5*vPros[i]*pLen]    
        #production well segments to surface
        pss = []
        for i in range(0,num):
            pss += [np.asarray([p1s[i][0],p1s[i][1],s.Well_TargetDepth_m+0.1])]
        #wellheads
        h = np.asarray([0.0,0.0,s.Well_Spacing_m])
        azn, dip, dis = sg.azn_dip(iss+h,iss)
        wells += [line(iss[0]+h[0],iss[1]+h[1],iss[2]+h[2],dis,azn,dip,'injector',sra,srb,src,rgh,0+n)]
        for i in range(0,num):
            azn, dip, dis = sg.azn_dip(pss[i]+h,pss[i])
            wells += [line(pss[i][0]+h[0],pss[i][1]+h[1],pss[i][2]+h[2],dis,azn,dip,'producer',sra,srb,src,rgh,1+i+n)]
        #place surface segments
        azn, dip, dis = sg.azn_dip(iss,ibs)
        wells += [line(iss[0],iss[1],iss[2],dis,azn,dip,'pipe',ira,irb,irc,rgh,0+n)]
        for i in range(0,num):
            azn, dip, dis = sg.azn_dip(pss[i],p1s[i])
            wells += [line(pss[i][0],pss[i][1],pss[i][2],dis,azn,dip,'pipe',ira,irb,irc,rgh,1+i+n)]
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
    # return
    # ************************************************************************
    return wells