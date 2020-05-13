import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import importlib

import potentialCircularVelocity
import initialConditionsMassBulge
import potentialPositions
import massDistribution
import bulgeDispersion
import initialConditionsMassSalpeter
import initialConditionsMassDisk

import sys
import os

#global default values
exePath = os.path.abspath('')
functionName = ''
showPlot=True

def main():
    global dataPath
    global showPlot
    global functionName
    if len(sys.argv) != 1:
        dataPath=exePath+sys.argv[1]
        showPlot = sys.argv[2].lower() == 'true'
        functionName = sys.argv[3]
        plotPath = exePath+sys.argv[4]
        module = importlib.import_module(functionName)
        getattr(module, functionName)(dataPath,showPlot,sys.argv[5:])
    if not showPlot:
        plt.savefig(plotPath+functionName+'.jpg')
        plt.close()

main()
