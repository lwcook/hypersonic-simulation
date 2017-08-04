#!/usr/bin/env python
from __future__ import division
import pdb
import numpy as np
from scipy.optimize import fsolve

from geometry import Geometry, Strip, Panel

import plotting
import matplotlib.pyplot as plt

from math import sin, cos, asin, exp, pi, atan, tan, log
from numpy import sqrt


class AeroModel(object):
    '''Class for the forces on an aerodynamic geometry'''

    def __init__(self, altitude=20, M=8, alpha=0.0, dynamic_pressure=None):

        ## Sea level properties
        self.T0 = 288.15 # Temperature in
        self.p0 = 101325 # pressure in Pascals
        self.rho0 = 1.225 # density in kg/m^2
        self.a0 = 340.3 # speed of sound in m/s
        self.g0 = 9.80665 # acceleration due to gravity m/s^2
        self.rE = 6378137 # Earth's mean radius at the equator
        self.R = 287.058
        self.gamma = 1.4

        ## Orientation of lift and drag vectors
        alpha = alpha*np.pi/180.
        Rotmat = np.array([ [np.cos(-1*alpha), 0, np.sin(-1*alpha)],
                            [0, 1, 0],
                            [-1*np.sin(-1*alpha), 0, np.cos(-1*alpha)]])
        self.drag_nvec = Rotmat.dot(np.array([1, 0, 0]))
        self.lift_nvec = Rotmat.dot(np.array([0, 0, 1]))
        self.fs_nvec = -1*self.drag_nvec

        assert(abs(np.dot(self.drag_nvec, self.lift_nvec)) < 1e-6)

        if dynamic_pressure is not None:

            def fun(h):
                p, T, rho = self.getAtmosphericProperties(h)
                a = np.sqrt(self.gamma*self.R*T)
                q = 0.5*rho*(M*a)**2
                return q - dynamic_pressure

            altitude = fsolve(fun, x0=altitude)

        ## Free stream properties
        self.p_fs, self.T_fs, self.rho_fs = \
            self.getAtmosphericProperties(altitude)
        self.M_fs = M
        self.a_fs = sqrt(self.gamma*self.R*self.T_fs)
        self.dynamic_pressure = 0.5*self.rho_fs*(self.M_fs*self.a_fs)**2


    def analyze_geometries(self, geoms):
        try:
            iter(geoms)
        except TypeError:
            geoms = [geoms]

        Lift = 0
        Drag = 0
        for geom in geoms:
            L, D = self.analyze_geometry(geom)
            Lift += L
            Drag += D
            print 'Lift: ', L, '    Drag: ', D

        return Lift, Drag

    def analyze_geometry(self, geom, coeffs=False):
        if not type(geom) == Geometry:
            raise TypeError('Analyze requires instance of Geometry class')

        Lift = 0
        Drag = 0
        for strip in geom:
            L, D, strip_info = self.strip_analysis(strip, geom.tangent_method)
            Lift += L
            Drag += D

        if coeffs:
            Cl = Lift/(self.dynamic_pressure*geom.ref_area)
            Cd = Drag/(self.dynamic_pressure*geom.ref_area)
            return Cl, Cd
        else:
            return Lift, Drag

    def strip_analysis(self, strip, tangent_method='cone'):
        if not type(strip) == Strip:
            raise TypeError('Strip analysis requries instance of Strip class')

        plifts, flifts = np.zeros(len(strip)), np.zeros(len(strip))
        pdrags, fdrags = np.zeros(len(strip)), np.zeros(len(strip))
        M1, p1, T1, rho1, i1 = self.M_fs, self.p_fs, self.T_fs, self.rho_fs, 0.
        pstrip, Mstrip = np.zeros(len(strip)), np.zeros(len(strip))
        istrip, betastrip = np.zeros(len(strip)), np.zeros(len(strip))
        sstrip, zstrip = np.zeros(len(strip)), np.zeros(len(strip))

        in_shadow = None
        turbulent_blayer = False
        for ip, panel in enumerate(strip):

            ## Only work with non-zero area panels
            if panel.getArea() < 1e-8:
                continue

            panel_location = np.dot(panel.getCentroid(), self.drag_nvec)
            norm = panel.getNormalVector()
            sign = np.sign(np.dot(norm, self.drag_nvec))

            ## Determine whether panel is in shadow
            if in_shadow is None:
                if np.dot(norm, self.drag_nvec) > 0:
                    in_shadow = True
                else:
                    in_shadow = False
            else:
                if sign == 0:
                    sign = 1
                if sign != prev_sign:
                    in_shadow = not in_shadow

            ## Uset tangent method or PM expansion to find local properties
            incidence = np.pi/2 - np.arccos(
                        np.dot(panel.getNormalVector(), self.fs_nvec))

            if False:
                print 'In Shadow: ', in_shadow
                print 'Incidence: ', incidence
                print 'Area: ', panel.getArea()
                print 'Norm: ', panel.getNormalVector()
                print 'M, p, T, rho: ', M1, p1, T1, rho1
                print ' '

            if not in_shadow:
                if abs(incidence) > 1e-6:
                    if tangent_method.lower() == 'wedge':
                        M, p, T, rho, beta = self.tangentWedge(incidence)
                    else:
                        M, p, T, rho, beta = self.tangentCone(incidence)
                        Mw, pw, Tw, rhow, betaw = self.tangentWedge(incidence)
#                        pdb.set_trace()
                else:
                    M, p, T, rho, beta = M1, p1, T1, rho1, np.pi/2
            else:
                del_inc = abs(incidence - i1)
                M, p, T, rho = self.PMExpansion(del_inc, M1, p1, T1, rho1)

            ## Determine whether boundary layer transitions
            if not turbulent_blayer:
                transition_point = self.blayerTransition(M, p, T, rho)
                if panel_location > transition_point:
                    turbulent_blayer = True

            ## Evaluate skin friction
            tau = self.referenceTempFriction(panel_location, M, p, T, rho,
                    turbulent_blayer, tangent_method)

            ## Update previous panel's properties
            M1, p1, T1, rho1, i1 = M, p, T, rho, incidence
            prev_norm, prev_sign = norm, sign

            pstrip[ip], Mstrip[ip], istrip[ip] = p, M, incidence
            sstrip[ip] = in_shadow

            ## Evaluate panel forces in format (lift, drag)
            plifts[ip] = p*(-1*panel.getAreaVector()).dot(self.lift_nvec)
            pdrags[ip] = p*(-1*panel.getAreaVector()).dot(self.drag_nvec)
            flifts[ip] = tau*(panel.getArea())*sin(incidence)
            fdrags[ip] = tau*(panel.getArea())*cos(incidence)

#        plotting.plot_geometry(Geometry([strip]))
#        plt.show()
#        pdb.set_trace()

        pLift, fLift = sum(plifts), sum(flifts)
        pDrag, fDrag = sum(pdrags), sum(fdrags)
        info = {'Pressures': pstrip, 'Machs': Mstrip}

        return pLift + fLift, pDrag + fDrag, info

    def blayerTransition(self, M, p, T, rho):
        mu = sutherlandViscosity(T)
        logRe = 6.421*np.exp(1.209e-4 * M**2.641)
        ReT = 10**logRe
        return ReT*mu/(rho*M*np.sqrt(self.gamma*self.R*T))

    def referenceTempFriction(self, x, M, p, T, rho, turb, method):
        g = self.gamma
        mu = sutherlandViscosity(T)
        Re = M*np.sqrt(g*self.R*T)*rho*x / mu

        if not turb:
            Pr = 0.72
            r = Pr**(1./2.)
            Tw = T*(1 + r*((g-1)/2)*M**2)
            Tref = T*(0.45 + 0.55*(Tw/T) + 0.16*r*((g-1)/2.)*M**2)
        else:
            Pr = 0.9
            r = Pr**(1./3.)
            Tw = T*(1 + r*((g-1)/2)*M**2)
            Tref = T*(0.5 + 0.5*(Tw/T) + 0.16*r*((g-1)/2.)*M**2)

        rhoref = p/(self.R*Tref)
        muref = 1.458*1e-6 * Tref**(3./2)/(Tref + 110.4)
        Reref = M*np.sqrt(g*self.R*Tref)*rhoref*x/muref

        if not turb:
            Cf = 0.664/np.sqrt(Reref)
        else:
            Cf = 0.02296/(Re**0.139) * (rhoref/rho)**0.861 * (muref/mu)**0.139

        tau = Cf*0.5*self.rho_fs*(self.M_fs*self.a_fs)**2

        return tau

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

    def PMExpansion(self, del_theta, M1, p1, T1, rho1):

        if del_theta == 0:
            return M1, p1, T1, rho1
        else:

            ## Prandtl Meyer function solving
            nu1 = prandtlMeyer(M1, self.gamma)
            nu2 = del_theta + nu1
            fun = lambda M: prandtlMeyer(M, self.gamma) - nu2
            fprime = lambda M: prandtlMeyerGradient(M, self.gamma)
#            print fprime(M1)
#            print (fun(M1+1e-6)-fun(M1))/1e-6
            M2 = float(fsolve(fun, x0=M1, fprime=fprime))
            g = self.gamma

            ## Isentropic expansion, from Anderson Section 9.6
            T2 = T1*(1+((g-1)/2.)*M1**2)/(1+((g-1)/2.)*M2**2)
            p2 = p1*(T2/T1)**( g / (g-1) )
            rho2 = rho1*(T2/T1)**( 1/(g-1) )

            rho3 = p2 / (self.R * T2)
            rho0 = p1 / (self.R * T1)

            return M2, p2, T2, rho2



    def tangentCone(self, theta):
        '''Flow across an conical shock wave giving surface properties.

        Using Rasmussen (1967) - this is an approximation for high Mach nubmers
        and small cone angles.
        For exact solutions, should solve the Taylor Maccoll (1933) equations
        directly, which requires a numerical solution or a lookup table.
        '''
        g = self.gamma
        th = theta
        sinth = sin(th)
        costh = cos(th)

        beta = asin( sinth*sqrt(((g+1)/2) + (self.M_fs*sinth)**(-2)) )
        b = beta
        sinb = sin(b)
        Msq = self.M_fs**2
        cosb = cos(b)

        M_norm = self.M_fs*sin(b)

        __ = ( cosb**2*(1 + 2*(b-th)) ) / \
                (1 + ((g-1)/2)*Msq*(sin(b)**2 - 2*(b-th)**2*costh**2) )
        M = sqrt( Msq * __ )

        p_other = self.p_fs*(1 + (2*g/(g+1))*(M_norm**2-1) )

        __ = 1 + (((g+1)*Msq*sinth**2 + 2)/((g-1)*Msq*sinth**2 + 2))* \
                log(((g+1)/2) + 1/(self.M_fs*sinth)**2)
        Cp = __*sinth**2
        p = Cp*0.5*self.rho_fs*(self.M_fs*self.a_fs)**2

        T = self.T_fs*(( 1 + (2*g/g+1))*M_norm**2-1) * \
                ( (2+(g-1)*M_norm**2)/((g+1)*M_norm**2) )

        rho = self.rho_fs*((g+1)*M_norm**2)/(2+(g-1)*M_norm**2)

        return M, p, T, rho, beta

    def tangentWedge(self, theta):
        '''Flow solution across an obique shock wave.

        Anderson. Fundamentals of Aerodynamics, 1991'''

        ## Checked against shock tables - working

        g = self.gamma
        th = theta

        beta = betaFromTM(th, self.M_fs, self.gamma)
        M_norm = self.M_fs*sin(beta);

        M = sqrt((1 + ((g-1)/2)*M_norm**2) / (g*M_norm**2 - (g-1)/2)) / \
             sin(beta-th);

        p = self.p_fs*(1 + (2*g/(g+1))*(M_norm**2 - 1))
        rho = self.rho_fs*((g+1)*M_norm**2)/(2 + (g-1)*M_norm**2)
        T = self.T_fs * (p/self.p_fs) * (self.rho_fs/rho);

        return M, p, T, rho, beta


def betaFromTM(theta, M, gamma):
    '''Obtains beta (radians) from theta (radians) and M in closed form.

    Rudd and Lewis. Journal of Aircraft Vol. 35, No. 4, 1998'''

    gam = gamma
    n = 0  # weak shock
    mu = asin(1/M)  # Mach wave angle
    c = tan(mu)**2
    a = (( gam-1)/2+(gam+1)*c/2)*tan(theta)
    b = (( gam+1)/2+(gam+3)*c/2)*tan(theta)
    d = sqrt(4*(1-3*a*b)**3/((27*a**2*c+9*a*b-2)**2)-1)

    return atan((b+9*a*c)/(2*(1-3*a*b)) - \
        (d*(27*a**2*c+9*a*b-2))/(6*a*(1-3*a*b)) * \
        tan(n*pi/3+1/3*atan(1/d)))

def thetaFromBM(beta, M, gamma):
    '''Obtains theta (radians) from beta (radians) and M in closed form.

    Anderson. Fundamentals of Aerodynamics, 1991'''
    g = gamma
    __ = (M**2*sin(beta)**2 - 1) / (M**2*(g + cos(2*beta)) + 2)

    return atan(2*(1./tan(beta))*__)

def sutherlandViscosity(T=298):
    mu0 = 1.827*10**-5
    T0 = 291.15
    C = 120. # Sutherland's constant
    mu = mu0*((T0 + C)/(T + C))*(T/T0)**(3./2)
    return mu

def prandtlMeyer(M, gamma=1.4):
    g = gamma
    nu = np.sqrt((g+1)/(g-1)) * np.arctan(np.sqrt(((g-1)/(g+1))*(M**2-1))) - \
            np.arctan(np.sqrt(M**2 -1))
    return nu

def prandtlMeyerGradient(M, gamma=1.4):
    A = np.sqrt((gamma+1)/(gamma-1))
    B = (gamma-1)/(gamma+1)

    grad = A*( 1/(1 + B*(M**2-1)) )*( B**0.5*M / (M**2-1)**0.5 ) - \
            (1/M**2)*(M/(M**2 - 1)**0.5)

    return grad

if __name__ == "__main__":
    pass
