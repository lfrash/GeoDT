# -*- coding: utf-8 -*-
""" 
enumeration of pipe types for cross-referencing
"""
#imports
import numpy as np
import os
if __package__ is None or __package__ == '':
    import vtk
else:
    from . import vtk

def build_box(s,path):
    #setup information
    size = 2*s.Domain_Size_m
    spacing = size/10.0
    depth = s.Well_TargetDepth_m
    grad = s.Rock_ThermalGradient_Kpm
    air = s.Air_Temperature_K
    num = int(2.0*size/spacing)+1
    dep = int((size+depth)/spacing)+1
    zspacing = (size+depth)/(dep-1)
    # pin = s.PIN_PIN_PIN
    # path = path + "/%s_" %(pin)
    label = 'temp_K'
    
    #define slices
    slices = []
    if False: #full 3D volume
        ns = [num,num,num]
        o0 = [-size,-size,-size]
        ss = [spacing,spacing,spacing]
        slices += [['volum',ns,o0,ss]]
    if True: #horizontal slice only
        ns = [num,num,3]
        o0 = [-size,-size,-spacing]
        ss = [spacing,spacing,spacing]
        slices += [['slice',ns,o0,ss]]
    if True: #back wall
        ns = [num,3,dep]
        o0 = [-size,size-spacing,-size]
        ss = [spacing,spacing,zspacing]
        slices += [['yback',ns,o0,ss]]
    if True: #side wall
        ns = [3,num,dep]
        o0 = [-size-spacing,-size,-size]
        ss = [spacing,spacing,zspacing]
        slices += [['xback',ns,o0,ss]]
    if True: #base
        ns = [num,num,3]
        o0 = [-size,-size,-size-spacing]
        ss = [spacing,spacing,spacing]
        slices += [['zbase',ns,o0,ss]]
    if True: #toe
        ns = [3,num,dep]
        o0 = [size-spacing,-size,-size]
        ss = [spacing,spacing,zspacing]
        slices += [['xfron',ns,o0,ss]]
    if True: #front
        ns = [num,3,dep]
        o0 = [-size,size-spacing,-size]
        ss = [spacing,spacing,zspacing]
        slices += [['yfron',ns,o0,ss]]
    if True: #surface
        ns = [num,num,3]
        o0 = [-size,-size,depth-spacing]
        ss = [spacing,spacing,spacing]
        slices += [['zsurf',ns,o0,ss]]
    if True: #soccer pitch
        ns = [2,2,2]
        o0 = [-52.5,-34.0,depth-1.0]
        ss = [105.0,68.0,2.0]
        slices += [['pitch',ns,o0,ss]]
    if True: #centerpitch
        ns = [2,2,2]
        o0 = [-52.5,-34.0,0.5*depth-0.5*size-1.0]
        ss = [105.0,68.0,2.0]
        slices += [['pivot',ns,o0,ss]]
    
    #create vtk files
    for s in slices:
        #fetch
        name = s[0]
        ns = s[1]
        o0 = s[2]
        ss = s[3]
        
        #get temperatures
        data = np.zeros((ns[0],ns[1],ns[2]),dtype=float)
        for i in range(0,ns[2]):
            data[:,:,i] = (depth - (o0[2] + spacing*i)) * grad + air
        
        #zero ground level
        o0[2] += -depth
        
        #build vtk file
        fname = path + name + '.vtk'
        vtk.structVtk(ns,o0,ss,data,label,fname)

"""
def build_pts(self,spacing=25.0,fname='test_gridx',special=''):
    print( '*** constructing temperature grid ***')
    #structured grid of datapoints
    fname = fname + '_therm.vtk'
    size = self.rock.size #TODO modify to focus on area of interest
    num = int(2.0*size/spacing)+1
    dep = int((size+self.rock.ResDepth)/spacing)+1
    label = 'temp_K'
    ns = [num,num,num]
    o0 = [-size,-size,-size]
    ss = [spacing,spacing,spacing]
    if special == 'z': #horizontal slice only
        ns = [num,num,3]
        o0 = [-size,-size,-spacing]
        ss = [spacing,spacing,spacing]
    elif special == 'y': #back wall
        ns = [num,3,dep]
        o0 = [-size,size-2*spacing,-size]
        ss = [spacing,spacing,spacing]
    elif special == 'x': #side wall
        ns = [3,num,dep]
        o0 = [-size,-size,-size]
        ss = [spacing,spacing,spacing]
    elif special == 'b': #base
        ns = [num,num,3]
        o0 = [-size,-size,-size]
        ss = [spacing,spacing,spacing]
    elif special == 't': #toe
        ns = [3,num,dep]
        o0 = [size-2*spacing,-size,-size]
        ss = [spacing,spacing,spacing]
    elif special == 's': #surface
        ns = [num,num,3]
        o0 = [-size,-size,self.rock.ResDepth-2*spacing]
        ss = [spacing,spacing,spacing]
    
    
    #initialize data to initial rock temperature
    #data = np.ones((num,num,num),dtype=float)*self.rock.BH_T
    # data = np.ones((ns[0],ns[1],ns[2]),dtype=float)*self.rock.BH_T
    data = np.zeros((ns[0],ns[1],ns[2]),dtype=float)
    if self.rock.gradient:
        for i in range(0,ns[2]):
            # data[:,:,i] = data[:,:,i] + (size - spacing*i) * self.rock.ResGradient
            data[:,:,i] = (self.rock.ResDepth - (o0[2] + spacing*i)) * self.rock.ResGradient + self.rock.AmbTempK
    
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
        #get points associated with this feature
        try:
            xR = [np.max([int((c0[0]-r0-o0[0])/spacing)-1,0]), np.min([int((c0[0]+r0-o0[0])/spacing)+1,ns[0]])]
            yR = [np.max([int((c0[1]-r0-o0[1])/spacing)-1,0]), np.min([int((c0[1]+r0-o0[1])/spacing)+1,ns[1]])]
            zR = [np.max([int((c0[2]-r0-o0[2])/spacing)-1,0]), np.min([int((c0[2]+r0-o0[2])/spacing)+1,ns[2]])]
        except:
            xR = [0,1]; yR = [0,1]; zR = [0,1]
        print([xR,yR,zR])
        #pipes and wells
        if (int(self.pipes.typ[i]) in [typ('injector'),typ('producer'),typ('pipe'),typ('perfcluster'),typ('screen')]): #pipe, Hazen-Williams
            #well info
            dip = self.wells[self.pipes.fID[i]].dip #dip
            azn = self.wells[self.pipes.fID[i]].azn #azimuth
            vAxi = np.asarray([math.sin(azn)*math.cos(-dip),math.cos(azn)*math.cos(-dip),math.sin(-dip)]) #axial vector
            rc = self.wells[self.pipes.fID[i]].rc #azimuth
            #cycle thruogh all points
            # for x in range(0,len(data)):
            #     for y in range(0,len(data[0])):
            #         for z in range(0,len(data[0,0])):
            #cycle only through local points
            for x in range(xR[0],xR[1]):
                for y in range(yR[0],yR[1]):
                    for z in range(zR[0],zR[1]):
                        #point coordinates
                        xPt = np.asarray([o0[0]+x*spacing, o0[1]+y*spacing, o0[2]+z*spacing])
                        #spherical radius vector
                        pi = xPt-c0
                        #lengthwise distance along well
                        li = np.linalg.norm(np.dot(pi,vAxi))
                        #lengthwise distance from origin x0 as weighting factor
                        lw = np.linalg.norm(np.dot((xPt-x0),vAxi))/(2.0*r0)
                        #radial distance from well
                        ri = ((np.linalg.norm(pi))**2.0 - li**2.0)**0.5
                        #if within length and thermal radius
                        if (li < r0) and (ri < R0) and (ri > rc):
                            #subtract delta T at point based on distance of point versus thermal radius
                            # data[x,y,z] = data[x,y,z] + (0.5*(T1+T0-Tr1-Tr0))*(1.0-(np.log(rc/ri)/np.log(rc/R0)))
                            data[x,y,z] = data[x,y,z] + ((1.0-lw)*(T0-Tr0)+(lw)*(T1-Tr1))*(1.0-(np.log(rc/ri)/np.log(rc/R0)))
                        #if point is inside the well
                        elif (li < r0) and (ri <= rc):
                            # data[x,y,z] = data[x,y,z] + (0.5*(T1+T0-Tr1-Tr0))
                            data[x,y,z] = data[x,y,z] + ((1.0-lw)*(T0-Tr0)+(lw)*(T1-Tr1))
            
        #fractures and planes
        elif (int(self.pipes.typ[i]) in [typ('fracture'),typ('propped')]): #,typ('choke')]): #(int())Y[i][2] == 1: #fracture, effective cubic law
            #fracture info
            dip = self.faces[self.pipes.fID[i]].dip
            azn = self.faces[self.pipes.fID[i]].str
            vNor = np.asarray([math.sin(azn+90.0*deg)*math.sin(dip),math.cos(azn+90.0*deg)*math.sin(dip),math.cos(dip)])
            vLeg = (x1-c0)/np.linalg.norm(x1-c0)
            vWid = np.cross(vNor,vLeg) 
            #cycle thruogh all points
            # for x in range(0,len(data)):
            #     for y in range(0,len(data[0])):
            #         for z in range(0,len(data[0,0])):
            #cycle only through local points
            for x in range(xR[0],xR[1]):
                for y in range(yR[0],yR[1]):
                    for z in range(zR[0],zR[1]):
                        #point coordinates
                        xPt = np.asarray([o0[0]+x*spacing, o0[1]+y*spacing, o0[2]+z*spacing])
                        #spherical radius vector
                        pi = xPt-c0
                        #normal distance from fracture
                        ni = np.linalg.norm(np.dot(pi,vNor))
                        #lengthwise distance from fracture
                        li = np.linalg.norm(np.dot(pi,vLeg))
                        #lengthwise distance from origin x0 as weighting factor
                        lw = np.linalg.norm(np.dot((xPt-x0),vLeg))/(2.0*r0)
                        #widthwise distance from fracture
                        wi = np.linalg.norm(np.dot(pi,vWid))
                        #if within length, width, normal
                        if (ni <= R0) and (li <= r0) and (wi <= (0.5*self.pipes.W[i])):
                            #subtract delta T at point based on distance of point versus thermal radius
                            # data[x,y,z] = data[x,y,z] + (0.5*(T1+T0)-self.rock.BH_T)*(1.0-ni/R0)
                            # data[x,y,z] = data[x,y,z] + ((1.0-lw)*(T0-Tr0)+(lw)*(T1-Tr1))*(1.0-ni/R0)
                            Tri = self.rock.BH_T - xPt[2]*self.rock.ResGradient
                            data[x,y,z] = data[x,y,z] + ((1.0-lw)*T0+(lw)*T1-Tri)*(1.0-ni/R0)

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
        """
