from geometry import Geometry

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

blue = [50./255., 100./255., 180./255.]
red = [180./255., 34./255., 34./255.]
green = [34./255., 139./255., 34./255.]
grey = [150./255., 150./255., 150./255.]

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

def mpl2tex(figsize=[2.8, 2.2]):
        '''Makes matplotlib's parameters give a plot formatted
        for latex documents
        -figsize: the figures dimension's in inches.
        Default is small enough for single column'''

        # Note: 1.0/72.27 inches per pt
        # [2.8,2.2] fits on the smallest of these (CMAME) and is a good ratio
        # CMAME template has a 390.0 pt wide textwidth - 5.396 inches
        # Thesis: 437.46 - 6.05 inches

        mpl.rcParams.update({"figure.figsize": figsize,
                             "font.family": "serif",
                             "text.usetex": True,
                             "text.latex.preamble": r"\usepackage{amsmath}",
                             "font.size": 8,
                             "font.weight": "light",
                             'axes.labelsize': 9,
                             'axes.titlesize': 8,
                             'legend.fontsize': 8,
                             'xtick.labelsize': 8,
                             'ytick.labelsize': 8,
                             'lines.linewidth': 0.6,
                             'axes.linewidth': 0.75,
                             'patch.linewidth': 0.75,
                             'legend.fontsize': 'medium',
                             'legend.scatterpoints': 1
                             })

def savefig(name='saved_fig', bSaveBase=False,
            base='/phd-thesis/Figs/', bSaveData=False, formatstr='pdf'):
    '''Function that saves the plot as well as the
    underlying data of the currently open figure:
    -name: string that the figure is saved as'''

    date = datetime.datetime.now()
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    subprocess.call(["mkdir", "-p", "./figs/"])
#    plt.savefig('./output/' + str(name) + '.pdf', format='pdf')
    plt.savefig('./figs/' +  str(name) + '_' + str(date.day) +
                months[date.month-1] + '.' + formatstr, format=formatstr)

if __name__ == "__main__":
    pass
