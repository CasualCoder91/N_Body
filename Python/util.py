import numpy as np
import matplotlib.pyplot as plt

#calculate needed degrees to show whole cluster
clusterDist = 6000 #pc
clusterSize = 1.1 #pc diameter
angle = 2*np.degrees(np.arctan2(0.5*1.1,clusterDist))
print("Minimum angle required: ", angle)
