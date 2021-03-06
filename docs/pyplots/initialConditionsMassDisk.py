import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def initialConditionsMassDisk(dataPath='', showPlot=True, arguments=[]):
    count = 10000
    if len(arguments)>0:
        count = int(arguments[0])

    data = np.loadtxt('initialConditionsMassDisk'+str(count)+'.dat',delimiter=',')
    x = data[:,1]
    fig = plt.figure()

    # histogram on linear scale
    plt.subplots_adjust(hspace = 0.4)
    plt.subplot(211)
    plt.title('Stars with total mass of 10e3 $M_{\odot}$ sampled from disk PDMF')
    plt.xlabel('mass [$M_{\odot}$]')
    plt.ylabel('Count')
    hist, bins, _ = plt.hist(x, bins=200)

    # histogram on log scale.
    # Use non-equal bin sizes, such that they look equal on log scale.
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.subplot(212)
    plt.hist(x, bins=logbins)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('mass [$M_{\odot}$]')
    plt.ylabel('Count')
    if showPlot:
        plt.show()

if __name__ == "__main__":
    initialConditionsMassDisk()
