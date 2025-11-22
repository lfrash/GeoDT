""" 
Useful unit conversions for GeoDT
"""
#initializations
import numpy as np
import math
import pylab
if __package__ is None or __package__ == '':
    from units import *
else:
    from .units import *

#initializations
deg=deg
MPa=MPa

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
            nor_dir_radians = dip_dir_radians + np.pi
            pylab.pcolormesh(nor_dir_radians,dip_angle_deg,np.ma.masked_where(np.isnan(criticalDelPpG),criticalDelPpG), 
                        vmin=self.Sh*0.5, vmax=np.max([self.SV,self.SH]),cmap='rainbow_r')
            ax.plot(str_dip[0,0]+np.pi/2,str_dip[0,1]/deg,'.k')
            ax.plot(str_dip[1,0]+np.pi/2,str_dip[1,1]/deg,'ob')
            ax.plot(str_dip[2,0]+np.pi/2,str_dip[2,1]/deg,'xc')
            ax.grid(True)
            ax.set_rgrids([0,30,60,90],labels=[])
            ax.set_thetagrids([0,90,180,270])
            pylab.colorbar()
            ax.set_title("Critical pressure (pole plot)", va='bottom')
            pylab.savefig('check.png', format='png',dpi=128)
            pylab.close()
        return str_dip

#verification
if False:
    #shmin 10 deg down and SSE and normal faulting
    s123 = np.zeros(3,dtype=float)
    s123[0] = 2700.0*9.81*5000.0
    s123[2] = 0.5*s123[0]
    s123[1] = 1.5*s123[2]
    tensor = cauchy()
    tensor.set_sigG_from_Principal(s123[2], s123[1], s123[0], 157.5*deg, 10.0*deg)
    str_dip = tensor.get_conjugates(plots=True)
    s123.sort(axis=0)