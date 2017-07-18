#!/usr/bin/env python
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class geom:
    """A class for geometries"""

    # makecone creates a cone geometry by default
    def makecone(self, phi = 0.1, chord = 1):
        npts = self.npoints
        nlns = self.nlines
        for ip in range(npts):
            for il in range(nlns):
                xc = chord*float(ip)/float(npts-1)
                r = xc*math.tan(phi)
                theta = 2*math.pi*float(il)/float(nlns-1)
                self.x[il,ip] = xc
                self.y[il,ip] = r*math.cos(theta)
                self.z[il,ip] = r*math.sin(theta)

    def __init__(self, fidelity = 13):
        npts = int(fidelity) # Number of points along a strip
        nlns = int(fidelity) # Number of lines bounding the strips
        self.npoints = npts
        self.nlines = nlns
        # Strips are stored along 2nd axis (horizontal axis)
        self.x,self.y = np.zeros((nlns,npts)), np.zeros((nlns,npts))
        self.z = np.zeros((nlns,npts))
        self.makecone()




if __name__ == "__main__":
    g1 = geom(13)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(g1.x,g1.y,g1.z)
    plt.show()
