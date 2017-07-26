from geometry import Geometry

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_geometry(geometry, ax=None, rescale=True):
    '''Plots an instance of the Geometry class'''

    if not type(geometry) == Geometry:
        raise ValueError('Input must be of the geometry class')

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    X, Y, Z = geometry.toMatrices()
    ax.plot_wireframe(X, Y, Z)

    if rescale:
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(),
                              Z.max()-Z.min()]).max() / 2.0
        mid_x = (X.max()+X.min()) * 0.5
        mid_y = (Y.max()+Y.min()) * 0.5
        mid_z = (Z.max()+Z.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
