from __future__ import division

import numpy as np
import math
import pdb

from geometry import Geometry, Panel, Strip

class WingedConeVehicle(object):
    '''Class that represents a WCV of multiple parts. Each part is an
    instance of the Geometry class, and the vehicle is simply a list of these
    parts. The vehicle is analyzed each part at a time using strip theory.

    We can do this because the analysis assumes the interaction between body
    parts is negligible when using strip theory with the tangent cone methd.
    '''

    def __init__(self, cylinder_R=1, cylinder_L=5, cone_L=7, boattail_L=1,
            boattail_R=0.5, wing_span=3.5, wing_chord=1, wing_origin_L=2.5,
            wing_finish_L=13, wing_thickness=0.3):

        ## Input Parameters
        self.cylinder_R = cylinder_R
        self.cylinder_L = cylinder_L
        self.cone_L = cone_L
        self.boattail_L = boattail_L
        self.boattail_R = boattail_R

        self.wing_span = wing_span
        self.wing_chord = wing_chord
        self.wing_origin_L = wing_origin_L
        self.wing_finish_L = wing_finish_L
        self.wing_thickness = wing_thickness

        # Derive geometry and mass quantities
        self._derive_geometry()
        self._derive_mass()

        # Create geometries
        self.body = self.makeBody(Nlines=11, Npoints=4)
        self.wing1 = self.makeWing(self.wing_span, self.wing_chord,
                                    self.wing_origin_L, self.wing_finish_L,
                                    self.wing_thickness, theta=np.pi/2,
                                    Nlines=3, Npoints=3)
        self.wing2 = self.makeWing(self.wing_span, self.wing_chord,
                                    self.wing_origin_L, self.wing_finish_L,
                                    self.wing_thickness, theta=3*np.pi/2,
                                    Nlines=3, Npoints=3)
        self.tail = self.makeWing(self.tail_span, self.tail_chord,
                                    self.tail_origin_L, self.tail_finish_L,
                                    self.tail_thickness, theta=np.pi)

        self.geometries = [self.body, self.wing1, self.wing2, self.tail]
#        self.geometries = [self.body, self.wing1, self.wing2]

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        self.body.plot(ax, rescale=True)
        self.wing1.plot(ax, rescale=False)
        self.wing2.plot(ax, rescale=False)
        self.tail.plot(ax, rescale=False)

    def _derive_geometry(self):

        # Derived measurements
        self.forebody_L = self.cylinder_L + self.cone_L
        self.total_L = self.forebody_L + self.boattail_L
        self.wing_chord_inner = self.wing_finish_L - self.wing_origin_L
        self.boattail_finish_L = self.forebody_L + self.boattail_L*(
                    (self.cylinder_R)/(self.cylinder_R - self.boattail_R))

        self.tail_span = 3
        self.tail_origin_L = self.total_L-2
        self.tail_chord = (self.total_L - self.tail_origin_L)*0.25
        self.tail_finish_L = self.total_L
        self.tail_thickness = self.wing_thickness*(
                    self.tail_chord/self.wing_chord)

    def _derive_mass(self):

        self.exposed_planform_area = 0
        self.internal_volume = 0
        self.inlet_capture_width = 0
        self.fuel_volume = 0

    # Scale the vehicle mass based on total surface area of vehicle

        mass_model_type = 1;

    #         wing_area_base = 58.3872;
    #         wing_area_ratio = self.wing_area/wing_area_base;
    #         surface_area_base = 93.887;
    #         surface_area_ratio = self.surface_area/surface_area_base;

        exposed_area_base = 20.1272;
        exposed_area_ratio = self.exposed_planform_area/exposed_area_base;

        internal_volume_base = 33.5431;
        internal_volume_ratio = self.internal_volume/internal_volume_base;

        scramjet_width_base = 0.64;
        scramjet_ratio = self.inlet_capture_width / scramjet_width_base;

        fuel_volume_base = 12.5952;
        fuel_volume_ratio = self.fuel_volume / fuel_volume_base;

        # Empty mass in kg
        m_fuselage_base = 2019;
        m_wings_base = 209;
        m_scramjets_base = 473;

        m_tank = 135;
        m_systems = 707;
        m_landing_gear = 142;

        m_fuselage = m_fuselage_base * internal_volume_ratio;
        m_wings = m_wings_base * exposed_area_ratio;
        m_scramjets = m_scramjets_base * scramjet_ratio;
        m_other = (m_tank*fuel_volume_ratio + m_systems + m_landing_gear);

    #         m_body_base = m_fuselage + m_tank + m_systems + m_landing_gear + m_scramjets; # reusable mass
    #         m_fuselage = m_body_base*surface_area_ratio;
    #         m_wings = m_wings_base*surface_area_ratio;

        m_vehicle = m_scramjets + m_fuselage + m_wings + m_other ; # Reusable vehicle mass

            # Payload and fuel mass
        m_payload = 1000;

        fuel_density = 71; # 71 kg/m^3 for hydrogen
        m_fuel = self.fuel_volume*fuel_density;

    #         fuselage_specific_mass = 61.5; # kg/m^3
    #         wing_specific_mass = 148.3; # kg/m^3 
    #         fuselage_mass = fuselage_specific_mass*self.internal_volume;
    #         wing_mass = wing_specific_mass*self.planform_area;
    #         scramjet_mass = 0.25*structural_mass;

        return m_vehicle, m_fuel, m_payload


    def makeCone(self, phi=0.1, chord=1, Nlines=20, Npoints=5):

#        # Strips are stored along 2nd axis (horizontal axis)
        x = np.zeros([Nlines, Npoints])
        y = np.zeros([Nlines, Npoints])
        z = np.zeros([Nlines, Npoints])

        for il in np.arange(Nlines):
            for ip in np.arange(Npoints):

                xc = chord*float(ip)/float(Npoints-1)
                r = xc*math.tan(phi)
                theta = 2*math.pi*float(il)/float(Nlines-1)*0.6

                x[il, ip] = xc
                y[il, ip] = -r*math.sin(theta)
                z[il, ip] = r*math.cos(theta)

        return Geometry(self.matricesToStrips(x, y, z))

    def makeBody(self, Nlines=21, Npoints=5):
        x = np.zeros([Nlines, 5])
        y = np.zeros([Nlines, 5])
        z = np.zeros([Nlines, 5])

        for il in np.arange(Nlines):
            theta = 2*math.pi*float(il)/float(Nlines-1)

            x[il, 0] = 0
            y[il, 0] = 0
            z[il, 0] = 0

            x[il, 1] = self.cone_L
            y[il, 1] = -self.cylinder_R*math.sin(theta)
            z[il, 1] = self.cylinder_R*math.cos(theta)

            x[il, 2] = self.cone_L + self.cylinder_L
            y[il, 2] = -self.cylinder_R*math.sin(theta)
            z[il, 2] = self.cylinder_R*math.cos(theta)

            x[il, 3] = self.cone_L + self.cylinder_L + self.boattail_L
            y[il, 3] = -self.boattail_R*math.sin(theta)
            z[il, 3] = self.boattail_R*math.cos(theta)

            x[il, 4] = self.cone_L + self.cylinder_L + self.boattail_L
            y[il, 4] = 0
            z[il, 4] = 0

        xmat = self._interpBody(x, Nlines, Npoints)
        ymat = self._interpBody(y, Nlines, Npoints)
        zmat = self._interpBody(z, Nlines, Npoints)

        return Geometry(self.matricesToStrips(xmat, ymat, zmat))

    def _interpBody(self, x, Nlines, Npoints):
        Np = Npoints
        NpB = int(np.ceil(Npoints/2))
        NpC = int(np.ceil(Npoints*1.5))
        xout = np.zeros([Nlines, NpC + Np + NpB + 1])
        for il in np.arange(Nlines):
            for ip in np.arange(NpC):
                xout[il, 0+ip] = x[il, 0] + (x[il,1]-x[il,0])*(ip/NpC)
            for ip in np.arange(Np):
                xout[il, NpC+ip] = x[il, 1] + (x[il, 2]-x[il,1])*(ip/Np)
            for ip in np.arange(NpB):
                xout[il, NpC+Np+ip] = x[il, 2]+(x[il, 3]-x[il, 2])*(ip/NpB)
            xout[il, -1] = x[il, 4]
        return xout

    def makeWing(self, wing_span, wing_chord, wing_origin_L, wing_finish_L,
                 wing_thickness, theta=0, Nlines=5, Npoints=4):

        wing_chord_inner = wing_finish_L - wing_origin_L

        x_inner = np.zeros(5)
        z_inner = np.zeros(5)
        x_outer = np.zeros(5)
        z_outer = np.zeros(5)

        tancone = self.cylinder_R / self.cone_L
        tanboat = (self.cylinder_R - self.boattail_R) / self.boattail_L

        w_o_L = wing_origin_L
        w_slope_L = wing_finish_L - w_o_L - wing_chord
        tanwing = wing_span / w_slope_L

        w_o_L2 = wing_origin_L + wing_chord_inner/2.
        w_slope_L2 = wing_finish_L - w_o_L2 - wing_chord/2.
        tanwing2 = wing_span / w_slope_L2

        b_o_L = self.boattail_finish_L - w_o_L
        b_o_L2 = self.boattail_finish_L - w_o_L2

        w0_z = min(w_o_L / ( (1/tancone) - (1/tanwing)), self.cylinder_R)
        w0_z = min(b_o_L / ( (1/tanwing) + (1/tanboat)), w0_z)
        w0_x = w0_z / tanwing
        w2_z = min(w_o_L2 / ( (1/tancone) - (1/tanwing2)), self.cylinder_R)
        w2_z = min(b_o_L2 / ( (1/tanwing2) + (1/tanboat)), w2_z)
        w2_x = w2_z / tanwing2

        ## Manually define the wing edges
        x_inner[0] = w_o_L + w0_x
        z_inner[0] = w0_z
        x_outer[0] = min(w_o_L + wing_span/tanwing, wing_finish_L)
        z_outer[0] = wing_span

        x_inner[1] = max(self.cone_L, x_inner[0])
        z_inner[1] = self.cylinder_R
        x_outer[1] = x_outer[0] + wing_chord/2.
        z_outer[1] = wing_span

        x_inner[2] = w_o_L2 + w2_x
        z_inner[2] = w2_z
        x_outer[2] = x_outer[0] + wing_chord/2.
        z_outer[2] = wing_span

        x_inner[3] = min(self.forebody_L, wing_finish_L)
        z_inner[3] = self.cylinder_R
        x_outer[3] = x_outer[0] + wing_chord/2.
        z_outer[3] = wing_span

        x_inner[4] = wing_finish_L
        z_inner[4] = self.cylinder_R - (wing_finish_L - self.forebody_L)/\
            self.boattail_L*(self.cylinder_R - self.boattail_R)
        x_outer[4] = wing_finish_L
        z_outer[4] = wing_span

        inds = np.argsort(x_inner)
        x_inner = x_inner[inds]
        z_inner = z_inner[inds]
        x_outer = x_outer[inds]
        z_outer = z_outer[inds]

        x_upper = np.zeros([Nlines, 5])
        y_upper = np.zeros([Nlines, 5])
        z_upper = np.zeros([Nlines, 5])
        x_lower = np.zeros([Nlines, 5])
        y_lower = np.zeros([Nlines, 5])
        z_lower = np.zeros([Nlines, 5])
        for il in np.arange(Nlines):
            for ip in np.arange(5):
            ## Interpolate between wing edges
                x_upper[il, ip] = x_inner[ip] + \
                    (x_outer[ip]-x_inner[ip])*(il/(Nlines-1))
                z_upper[il, ip] = z_inner[ip] + \
                    (z_outer[ip]-z_inner[ip])*(il/(Nlines-1))
                x_lower[-1-il, ip] = x_inner[ip] + \
                    (x_outer[ip]-x_inner[ip])*(il/(Nlines-1))
                z_lower[-1-il, ip] = z_inner[ip] + \
                    (z_outer[ip]-z_inner[ip])*(il/(Nlines-1))

            ## Add third dimension with getWingSection
            y_upper[il, :] = self.getWingSection(x_upper[il,:],
                    z_upper[il,:], wing_chord_inner, wing_chord,
                    wing_span, wing_origin_L, wing_finish_L, wing_thickness)
            y_lower[-1-il, :] = -1*self.getWingSection(x_lower[-1-il, :],
                    z_lower[-1-il, :], wing_chord_inner, wing_chord,
                    wing_span, wing_origin_L, wing_finish_L, wing_thickness)

        xl = self._interpWing(x_lower, Nlines, Npoints)
        xu = self._interpWing(x_upper, Nlines, Npoints)
        yl = self._interpWing(y_lower, Nlines, Npoints)
        yu = self._interpWing(y_upper, Nlines, Npoints)
        zl = self._interpWing(z_lower, Nlines, Npoints)
        zu = self._interpWing(z_upper, Nlines, Npoints)

        xmat = np.concatenate([xu, xl], axis=0)
        ymat = np.concatenate([yu, yl], axis=0)
        zmat = np.concatenate([zu, zl], axis=0)

        xrot = np.zeros(xmat.shape)
        yrot = np.zeros(ymat.shape)
        zrot = np.zeros(zmat.shape)
        for ii in np.arange(xmat.shape[0]):
            for jj in np.arange(xmat.shape[1]):
                vec = np.array([xmat[ii, jj], ymat[ii, jj], zmat[ii, jj]])
                RotMat = np.array([ [1, 0, 0],
                                    [0, np.cos(theta), -1*np.sin(theta)],
                                    [0, np.sin(theta), np.cos(theta)]])
                outvec = np.dot(RotMat, vec)
                xrot[ii, jj] = outvec[0]
                yrot[ii, jj] = outvec[1]
                zrot[ii, jj] = outvec[2]

        return Geometry(self.matricesToStrips(xrot, yrot, zrot))

    def _interpWing(self, x, Nlines, Npoints):
        Np1 = Npoints
        Np2 = int(np.ceil(Npoints/2))
        Np3 = int(np.ceil(Npoints))
        Np4 = int(np.ceil(Npoints))
        xout = np.zeros([Nlines, Np1 + Np2 + Np3 + Np4])
        for il in np.arange(Nlines):
            for ip in np.arange(Np1):
                N = 0+ip
                xout[il, N] = x[il, 0]+(x[il, 1]-x[il, 0])*(ip/Np1)
            for ip in np.arange(Np2):
                N = Np1+ip
                xout[il, N] = x[il, 1]+(x[il, 2]-x[il, 1])*(ip/Np2)
            for ip in np.arange(Np3):
                N = Np1+Np2+ip
                xout[il, N] = x[il, 2]+(x[il, 3]-x[il, 2])*(ip/Np3)
            for ip in np.arange(Np4):
                N = Np1+Np2+Np3+ip
                xout[il, N] = x[il, 3]+(x[il, 4]-x[il, 3])*(ip/Np4)
        return xout

    def getWingSection(self, xvec, zvec, wing_chord_inner, wing_chord,
                wing_span, wing_origin_L, wing_finish_L, wing_thickness):
        '''Input is a vector of xcoordinates, this returns the y-coordinates'''

        yvec = np.zeros(xvec.shape)

        for ix, (xi, zi) in enumerate(zip(xvec, zvec)):

            xc = wing_chord_inner - (wing_chord_inner-wing_chord)*zi/wing_span

            x0 = wing_origin_L + ((wing_finish_L-wing_chord-wing_origin_L)*
                                  zi/wing_span)

            assert(abs(x0 + xc - wing_finish_L) < 1e-6)

#            h = xc*0.3
            h = wing_thickness
            x = xi - x0

            if x/xc <= 0.5:
                yvec[ix] = h*x/xc
            else:
                yvec[ix] = h*(1 - x/xc)

        return yvec


    def matricesToStrips(self, x, y, z):
        Nlines, Npoints = x.shape

        strips = []
        for il in np.arange(Nlines-1):
            strip = []
            for ip in np.arange(Npoints-1):
                p1 = [x[il, ip], y[il, ip], z[il, ip]]
                p2 = [x[il, ip+1], y[il, ip+1], z[il, ip+1]]
                p3 = [x[il+1, ip+1], y[il+1, ip+1], z[il+1, ip+1]]
                p4 = [x[il+1, ip], y[il+1, ip], z[il+1, ip]]
                strip.append(Panel(p1, p2, p3, p4))

            strips.append(Strip(strip))

        return strips

if __name__ == "__main__":
    pass
