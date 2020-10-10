import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

def clusteringVelocity(dataPath='',showPlot=True,arguments=[]):
    data = np.loadtxt(dataPath+'clusteringVelocity.dat',delimiter=',')
    colors = np.where(data[:,2]==1,'y','k')
    plt.scatter(x=data[:,5], y=data[:,6],c=colors)
    #plt.yscale('log')
    #plt.ylabel('SMD ['+r'$M_{\odot}*pc^{-2}$]')
    #plt.xlabel('R [kpc]')
    plt.title('Proper Motion')

    if showPlot:
        plt.show()

if __name__ == "__main__":
    clusteringVelocity()
