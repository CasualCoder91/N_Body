import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def initialConditionsMassSalpeter(dataPath='',showPlot=True,arguments=[]):
    n = 10000
    if len(arguments)>0:
        n = int(arguments[0])
    data = np.loadtxt(dataPath+'initialConditionsMassSalpeter'+str(n)+'.dat',delimiter=',')
    x = data[:,1]
    fig = plt.figure(figsize=(8, 3))
    ax = plt.gca()
    # histogram on linear scale
    #plt.subplots_adjust(hspace = 0.4)
    #plt.subplot(211)
    #plt.title('10e3 Stars sampled from Salpeter IMF')
    plt.xlabel('mass [$M_{\odot}$]')
    plt.ylabel('Count')
    #hist, bins, _ = plt.hist(x, bins=400)

    # histogram on log scale.
    # Use non-equal bin sizes, such that they look equal on log scale.
    logbins = np.logspace(np.log10(x[0]),np.log10(np.amax(x)),400)
    #plt.subplot(212)
    plt.hist(x, bins=logbins)
    plt.xscale('log')
    #ax.set_xticks([25, 50, 75])
    plt.yscale('log')
    plt.xlabel('mass [$M_{\odot}$]')
    plt.ylabel('Count')
    plt.tight_layout()
    if showPlot:
        plt.show()

if __name__ == "__main__":
    initialConditionsMassSalpeter()
