import numpy as np
import matplotlib.pyplot as plt

#calculate needed degrees to show whole cluster
#param: clusterDist in pc, clusterSize: diameter of cluster in pc
def degNeeded(clusterDist,clusterSize):
    angle = 2*np.degrees(np.arctan2(0.5*1.1,clusterDist))
    print("Minimum angle required: ", angle)
    return

#calculate focus based on angle along xy plane with 0° at galatic center
#param: angle in deg
def calcFocus(angle):
    viewDist = 6000 #pc, r from sun
    viewPoint = np.array([[8300],[0],[27]]) #column vector, location of sun in pc
    x = -viewDist*np.cos(np.radians(angle))+viewPoint[0]
    y = viewDist*np.sin(np.radians(angle))+viewPoint[1]
    z = viewPoint[2]
    print("focus: ",x,y,z)
    return

degNeeded(6000,1.1)
calcFocus(15)
