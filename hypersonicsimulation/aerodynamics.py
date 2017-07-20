#!/usr/bin/env python
from __future__ import division
import pdb
import numpy as np
import geometry

from math import sin, cos, asin, exp, pi, atan, tan, log
from numpy import sqrt


class AeroModel(object):
    '''Class for the forces on an aerodynamic geometry'''

    def __init__(self):

        self.geometry = None

        self.T0 = 288.15 # Temperature in
        self.p0 = 101325 # pressure in Pascals
        self.rho0 = 1.225 # density in kg/m^2
        self.a0 = 340.3 # speed of sound in m/s
#        self.g0 = 9.80665 # acceleration due to gravity m/s^2
#        self.rE = 6378137 # Earth's mean radius at the equator
        self.R = 287.058
        self.gamma = 1.4

    def addGeometry(self, geometry):
        self.geometry = geometry

    def run(self):
        pass

    def strip_analysis(self, xstrip, ystrip, zstrip):
        pass

    def getAtmosphericProperties(self, altitude):

        fs_props = {}

        hkm = altitude # in kilometres
        hm = hkm*1000 # in metres
        if hkm < 0 or hkm > 80 :
            print('h outside of range')
            return

        if hkm <= 11:
            fs_props['T'] = self.self.T0 - 0.0065*hm
            fs_props['p'] = self.self.p0*(fs_props['T']/self.T0)**(5.2559)
        else:
            fs_props['T'] = 216
            fs_props['p'] = 22630*exp(-0.00015769*(hm-11000))

        fs_props['rho'] = fs_props['p']/(self.R*fs_props['T'])
#        fs_props['a'] = sqrt(self.gamma*self.R*fs_props['T'])
#        fs_props['mu'] = self.sutherland_viscosity(fs_props['T'])

        return fs_props['p'], fs_props['T'], fs_props['rho']

    def getSutherlandViscosity(self, T=298):
        mu0 = 1.827*10**-5
        T0 = 291.15
        C = 120. # Sutherland's constant
        mu = mu0*((T0 + C)/(T + C))*(T/T0)**(3./2)
        return mu

    def tangentCone(self, M_fs, p_fs, T_fs, rho_fs, theta):
        '''Flow across an conical shock wave giving surface properties.

        Using Rasmussen (1967) - this is an approximation for high Mach nubmers
        and small cone angles.
        For exact solutions, should solve the Taylor Maccoll (1933) equations
        directly, which requires a numerical solution or a lookup table.
        '''
        g = self.gamma
        a_fs = sqrt(g*self.R*T_fs)
        th = theta
        sinth = sin(th)
        costh = cos(th)

        beta = asin( sinth*sqrt(((g+1)/2) + (M_fs*sinth)**(-2)) )
        b = beta
        sinb = sin(b)
        Msq = M_fs**2
        cosb = cos(b)

        M_norm = M_fs*sin(b)

        __ = ( cosb**2*(1 + 2*(b-th)) ) / \
                (1 + ((g-1)/2)*Msq*(sin(b)**2 - 2*(b-th)**2*costh**2) )
        M = sqrt( Msq * __ )

        p_other = p_fs*(1 + (2*g/(g+1))*(M_norm**2-1) )

        __ = 1 + (((g+1)*Msq*sinth**2 + 2)/((g-1)*Msq*sinth**2 + 2))* \
                log(((g+1)/2) + 1/(M_fs*sinth)**2)
        Cp = __*sinth**2
        p = Cp*0.5*rho_fs*(M_fs*a_fs)**2

        T = T_fs*(( 1 + (2*g/g+1))*M_norm**2-1) * \
                ( (2+(g-1)*M_norm**2)/((g+1)*M_norm**2) )

        rho = rho_fs*((g+1)*M_norm**2)/(2+(g-1)*M_norm**2)
        pdb.set_trace()

        return M, p, T, rho, beta

    def tangentWedge(self, M_fs, p_fs, T_fs, rho_fs, theta):
        '''Flow solution across an obique shock wave.

        Anderson. Fundamentals of Aerodynamics, 1991'''

        ## Checked against shock tables - working

        g = self.gamma
        th = theta

        beta = self.betaFromTM(th, M_fs)
        M_norm = M_fs*sin(beta);

        M = sqrt((1 + ((g-1)/2)*M_norm**2) / (g*M_norm**2 - (g-1)/2)) / \
             sin(beta-th);

        p = p_fs*(1 + (2*g/(g+1))*(M_norm**2 - 1))
        rho = rho_fs*((g+1)*M_norm**2)/(2 + (g-1)*M_norm**2)
        T = T_fs * (p/p_fs) * (rho_fs/rho);

        return M, p, T, rho, beta


    def betaFromTM(self, theta, M):
        '''Obtains beta (radians) from theta (radians) and M in closed form.

        Rudd and Lewis. Journal of Aircraft Vol. 35, No. 4, 1998'''

        gam = self.gamma
        n = 0  # weak shock
        mu = asin(1/M)  # Mach wave angle
        c = tan(mu)**2
        a = (( gam-1)/2+(gam+1)*c/2)*tan(theta)
        b = (( gam+1)/2+(gam+3)*c/2)*tan(theta)
        d = sqrt(4*(1-3*a*b)**3/((27*a**2*c+9*a*b-2)**2)-1)

        return atan((b+9*a*c)/(2*(1-3*a*b)) - \
            (d*(27*a**2*c+9*a*b-2))/(6*a*(1-3*a*b)) * \
            tan(n*pi/3+1/3*atan(1/d)))

    def thetaFromBM(self, beta, M):
        '''Obtains theta (radians) from beta (radians) and M in closed form.

        Anderson. Fundamentals of Aerodynamics, 1991'''
        g = self.gamma
        __ = (M**2*sin(beta)**2 - 1) / (M**2*(g + cos(2*beta)) + 2)

        return atan(2*(1./tan(beta))*__)


if __name__ == "__main__":
    pass
