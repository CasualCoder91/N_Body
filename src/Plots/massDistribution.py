import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

def massDistribution(dataPath='',showPlot=True,arguments=[]):
    zParam = 1
    if len(arguments)>0:
        zParam = int(arguments[0])
    diskData = np.loadtxt(dataPath + 'massDistributionDisk_z'+str(zParam)+'.dat',delimiter=',')
    bulgeData = np.loadtxt(dataPath + 'massDistributionBulge_z'+str(zParam)+'.dat',delimiter=',')

    x = diskData[:,0]
    y = diskData[:,1]
    z = diskData[:,2]
    N = int(len(z)**.5)
    z = z.reshape(N, N)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(9.5, 5)
    fig.text(0.5, 0.04, 'x [kpc]', ha='center')
    fig.text(0.04, 0.5, 'y [kpc]', va='center', rotation='vertical')
    fig.suptitle('Mass distribution at z = 1.0 [kpc]', fontsize=14)

    ax1.title.set_text('Disk')
    pcm = ax1.imshow(z+10, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
            cmap=cm.hot, norm=LogNorm())
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(pcm,cax=cax,ax=ax1)
    cbar.ax.set_title('[$M_{\odot}$]')

    x = bulgeData[:,0]
    y = bulgeData[:,1]
    z = bulgeData[:,2]
    N = int(len(z)**.5)
    z = z.reshape(N, N)
    ax2.title.set_text('Bulge')
    pcm = ax2.imshow(z+10, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
            cmap=cm.hot, norm=LogNorm())
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(pcm,cax=cax,ax=ax2)
    cbar.ax.set_title('[$M_{\odot}$]')

    if showPlot:
        plt.show()

massDistribution()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(data[:,0], data[:,1], data[:,2], s=10,c='red')
