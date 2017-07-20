#!/usr/bin/env python
import numpy as np
import math
import pdb
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Panel(object):
    '''Class for individual panels for use with strip theory. 
    A panel consists of four 3D points stored as numpy arrays.

    The convention for points is from 1 to 4 clockwise around the panel looking
    inwards, with points 1 and 4 facing the leading edge of the vehicle.

    The area/normal vectors point outwards from the vehicle
    '''

    def __init__(self, p1, p2, p3, p4):
        self.p1 = np.array(p1)
        self.p2 = np.array(p2)
        self.p3 = np.array(p3)
        self.p4 = np.array(p4)
        for p in [self.p1, self.p2, self.p3, self.p4]:
            assert p.size == 3

    def __str__(self):
        return '['+str(self.p1)+', '+str(self.p2)+', '+\
               str(self.p3)+', '+str(self.p4)+']'''

    def __repr__(self):
        return '['+str(self.p1)+', '+str(self.p2)+', '+\
               str(self.p3)+', '+str(self.p4)+']'''

    def __add__(self, other):
        if type(other) is Panel:
            return Panel(self.p1 + other.p1, self.p2 + other.p2,
                        self.p3 + other.p3, self.p4 + other.p4)
        else:
            return Panel(self.p1 + other, self.p2 + other,
                        self.p3 + other, self.p4 + other)

    def getAreaVector(self):
        return 0.5*np.cross((self.p3-self.p1), (self.p2-self.p4))

    def getArea(self):
        return np.linalg.norm(self.getAreaVector())

    def getNormalVector(self):
        vec = self.getAreaVector()
        return vec/np.linalg.norm(vec)

    def getAngleWithVector(self, vector):
        norm_vec = self.getNormalVector()
        return np.pi/2 - math.acos(norm_vec.dot(vector))


class Strip(list):
    '''Class for a strip made out of individual panels for use with strip
    theory. This class is just a list of panels that ensure the panels link
    together to form a strip'''
    def __init__(self, listval):
        list.__init__(self, listval)
        for ip, panel in enumerate(listval):
            if type(panel) != Panel:
                raise ValueError('Strip must consist of a list of panels')
            if ip != 0:
                errmsg = '''Strip panels do not join up in longitudinal vehicle
                    direction to form a strip.'''
                if not np.allclose(panel.p1, listval[ip-1].p2):
                    raise ValueError(errmsg)
                if not np.allclose(panel.p4, listval[ip-1].p3):
                    raise ValueError(errmsg)

class Geometry(object):

    def __init__(self):

        self.strips = []

#        # Strips are stored along 2nd axis (horizontal axis)
#        self.x = np.zeros([nlns, npts])
#        self.y = np.zeros([nlns, npts])
#        self.z = np.zeros([nlns, npts])
#        self.makecone()

    # makecone creates a cone geometry by default
    def makecone(self, phi = 0.1, chord = 1, fidelity=13):

        for ip in np.arange(fidelity):
            for il in np.arange(fidelity):

                xc = chord*float(ip)/float(fidelity-1)
                r = xc*math.tan(phi)
                theta = 2*math.pi*float(il)/float(fidelity-1)

                self.x[il, ip] = xc
                self.y[il, ip] = r*math.cos(theta)
                self.z[il, ip] = r*math.sin(theta)

    def matricesToStrips(self, xmat, ymat, zmat):
        ## Matrix should be 2D, with strips stored along 2nd axis (horizontal)
        Nlines, Npoints = xmat.shape
        for ix in np.arange(Nlines):
            for iy in np.arange(Npoints)

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(g1.x,g1.y,g1.z)
        plt.show()


    def makeWCV(self):
        pass



if __name__ == "__main__":
    g1 = Geometry(13)
    g1.plot()

