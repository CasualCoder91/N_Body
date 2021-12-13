import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import os
import sys

def potentialPositions(dataPath='',showPlot=True,arguments=[]):
    sampleSize = 5000
    if len(arguments)>0:
        sampleSize = int(arguments[0])
    dataDisk = np.loadtxt(dataPath + 'potentialDiskPositionsSample'+str(sampleSize)+'.dat',delimiter=',')/1000
    dataBulge = np.loadtxt(dataPath + 'potentialBulgePositionsSample'+str(sampleSize)+'.dat',delimiter=',')/1000
    fig = plt.figure()
    #fig.set_size_inches(10, 5)
    ax = fig.add_subplot(1,2,1)
    plt.hist(dataDisk[:,0]+dataBulge[:,0], bins=200)
    #ax.title.set_text('Distribution along x Axis')
    plt.xlabel('$x_{GCA}$ [kpc]')
    #plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))

    ax = fig.add_subplot(1,2,2)
    #ax.title.set_text('Distribution along z Axis')
    plt.hist(dataDisk[:,2]+dataBulge[:,2], bins=200)
    plt.xlabel('$z_{GCA}$ [kpc]')
    ax.set_xlim([-5,5])
    #plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))

    #ax = fig.add_subplot(2,2,(3,4), projection='3d')
    #ax.scatter(dataDisk[:,0]*0.001, dataDisk[:,1]*0.001, dataDisk[:,2]*0.001, s=5,c='red')
    #ax.scatter(dataBulge[:,0]*0.001, dataBulge[:,1]*0.001, dataBulge[:,2]*0.001, s=5,c='blue')
    #ax.view_init(elev=5., azim=180)
    #ax.dist = 8
    #ax.set_xlim([-40,40])
    #ax.set_ylim([-40,40])
    #ax.set_zlim([-7,7])

    if showPlot:
        plt.show()

if __name__ == "__main__":
    potentialPositions()
