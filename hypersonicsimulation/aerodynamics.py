#!/usr/bin/env python
import numpy as np
import geometry
import math

class analysis_results:
    """Class to use a structure for analyses results"""
    def __init__(self):
        out = []

class aero_model:
    """Class for the forces on an aerodynamic geometry)"""

    def __init__(self, aero_geom = geometry.geom()):
        lift = 0
        drag = 0
        moment = 0
        self.geom = aero_geom # The geometry to be analyzed in here
        self.atmosphere = atmosphere_model() # The free stream properties
        self.atmosphere.evaluate_props(0)


    def aero_analysis(self):
        nlines = int(self.geom.x.shape[0])
        nstrips = nlines - 1 # Number of strips
        npoints = int(self.geom.x.shape[1]) 
        npanels = npoints -1 # Number of panels per strip

        for ii in range(nlines):
            xs, ys, zs = self.geom.x[ii,:], self.geom.y[ii,:], self.geom.z[ii,:]
            strip_results = self.strip_analysis(xs,ys,zs)


    def strip_analysis(self, xstrip, ystrip, zstrip):
        pass

    def tangent_cone(self, MachIn = 2.0, theta = 0.0):
        M1 = MachIn
        p1 = self.atmosphere.p
        gam = self.atmosphere.gamma
        temp1 = math.sin(theta)*math.sqrt(((gam+1)/2) +
                (M1*math.sin(theta))**(-2))
        beta = math.asin( temp1 )

        M1norm = M1*math.sin(beta)

        temp1 = (math.cos(beta)**2)*(1 + 2*(beta-theta))
        temp2 = 1 + ((gam-1)/2) * M1**2 * (math.sin(beta)**2 -
                2*(beta-theta)**2 *math.cos(theta))**2 

        M2 = math.sqrt( M1**2 * temp1/temp2 )

        print beta*180/math.pi
        print M2

        p2 = p1*(1 + (2*gam/(gam+1))*(M1norm**2-1) )

    def tangent_wedge(self, MachIn = 2.0, theta = 0.0):
        M1 = MachIn
        p1 = self.atmosphere.p
        gam = self.atmosphere.gamma








class atmosphere_model:
    """ Class for evaluating free stream properties in the atmosphere"""

    # Class variables
    # Properties at sea level
    T0 = 288.15 # Temperature in
    p0 = 101325 # pressure in Pascals
    rho0 = 1.225 # density in kg/m^2
    a0 = 340.3 # speed of sound in m/s
    g0 = 9.80665 # acceleration due to gravity m/s^2
    rE = 6378137 # Earth's mean radius at the equator

    def __init__(self):

        self.R = 287.058
        self.gamma = 1.4
        self.evaluate_props(0) # Initialize properties at sea level

    def evaluate_props(self, altitude):
        hkm = altitude # in kilometres
        hm = hkm*1000 # in metres
        if hkm < 0 or hkm > 80 :
            print('h outside of range')
            return

        if hkm <= 11:
            self.T = self.T0 - 0.0065*hm
            self.p = self.p0*(self.T/self.T0)**(5.2559)
        else:
            self.T = 216
            self.p = 22630*math.exp(-0.00015769*(hm-11000))

        self.rho = self.p/(self.R*self.T)
        self.a = math.sqrt(self.gamma*self.R*self.T)
        self.mu = self.sutherland_viscosity(self.T)

    def sutherland_viscosity(self,T = 298):
        mu0 = 1.827*10**-5
        T0 = 291.15
        C = 120 # Sutherland's constant

        mu = mu0*((T0 + C)/(T + C))*(T/T0)**(3/2)

        return mu




if __name__ == "__main__":
    a1 = aero_model()
    a1.tangent_cone(3,0.175)
