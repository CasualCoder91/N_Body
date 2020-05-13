import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

def potentialSurfaceDensity(dataPath='',showPlot=True,arguments=[]):
    dataBulge = np.loadtxt(dataPath+'potentialSurfaceDensityBulge.dat',delimiter=',')
    dataDisk = np.loadtxt(dataPath+'potentialSurfaceDensityDisk.dat',delimiter=',')

    plt.plot(dataDisk[:,0]*10**(-3), dataDisk[:,1], label='Disk')
    plt.plot(dataBulge[:,0]*10**(-3), dataBulge[:,1], label='Bulge')
    plt.yscale('log')
    plt.ylabel('SMD ['+r'$M_{\odot}*pc^{-2}$]')
    plt.xlabel('R [kpc]')
    plt.title('Surface Mass Density')
    plt.legend()
    if showPlot:
        plt.show()

potentialSurfaceDensity()
