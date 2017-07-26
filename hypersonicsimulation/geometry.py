#!/usr/bin/env python
from __future__ import division

import numpy as np
import math
import pdb

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

    def __sub__(self, other):
        if type(other) is Panel:
            return Panel(self.p1 - other.p1, self.p2 - other.p2,
                        self.p3 - other.p3, self.p4 - other.p4)
        else:
            return Panel(self.p1 - other, self.p2 - other,
                        self.p3 - other, self.p4 - other)

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
    together to form a strip.

    Convention is that the list starts at the leading edge of the vehicle and
    ends at the trailing edge. '''
    def __init__(self, list_of_panels):
        list.__init__(self, list_of_panels)
        for ip, panel in enumerate(list_of_panels):
            if type(panel) != Panel:
                raise ValueError('Strip must consist of a list of panels')

            errmsg = '''Strip panels do not join up in longitudinal vehicle
                direction to form a strip.'''
            if ip != 0:
                if not np.allclose(panel.p1, list_of_panels[ip-1].p2):
                    raise ValueError(errmsg)
                if not np.allclose(panel.p4, list_of_panels[ip-1].p3):
                    raise ValueError(errmsg)

class Geometry(list):
    '''Class for a geometry made out of strips for use with strip
    theory. This class is a list of strips that ensures the strips link
    together to form an enclosed 3D geometry.

    The convention is that strips are listed clockwise looking from in front
    of the vehicle towards the leading edge. '''

    def __init__(self, list_of_strips):
        list.__init__(self, list_of_strips)
        for ii, strip in enumerate(list_of_strips[0:-1]):
            if type(strip) != Strip:
                raise ValueError('Geometry must consist of a list of strips')

            errmsg = '''Strips do not join up side by side to form an enclosed
                3D geometry'''
            if ii == 0:
                pass
#                if not self._do_strips_join(strip, list_of_strips[-1]):
#                    raise ValueError(errmsg)
            else:
                if not self._do_strips_join(list_of_strips[ii-1], strip):
                    raise ValueError(errmsg)
                if not len(list_of_strips[ii-1]) == len(strip):
                    raise ValueError('Length of strips must all be equal')

    def _do_strips_join(self, strip1, strip2):
        for ii, (panel1, panel2) in enumerate(zip(strip1, strip2)):
            if not (np.allclose(panel1.p4, panel2.p1) and 
                    np.allclose(panel1.p3, panel2.p2)):
                return False

        return True

    def toMatrices(self):

        xmat = np.zeros([len(self)+1, len(self[0])+1])
        ymat = np.zeros([len(self)+1, len(self[0])+1])
        zmat = np.zeros([len(self)+1, len(self[0])+1])

        for il, strip in enumerate(self):
            for ip, panel in enumerate(strip):
                xmat[il, ip] = panel.p1[0]
                ymat[il, ip] = panel.p1[1]
                zmat[il, ip] = panel.p1[2]

                xmat[il, ip+1] = panel.p2[0]
                ymat[il, ip+1] = panel.p2[1]
                zmat[il, ip+1] = panel.p2[2]

        for ip, panel in enumerate(self[-1]):
            xmat[-1, ip] = panel.p4[0]
            ymat[-1, ip] = panel.p4[1]
            zmat[-1, ip] = panel.p4[2]

            xmat[-1, ip+1] = panel.p3[0]
            ymat[-1, ip+1] = panel.p3[1]
            zmat[-1, ip+1] = panel.p3[2]

        return xmat, ymat, zmat


if __name__ == "__main__":
    pass
