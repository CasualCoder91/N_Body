import sqlite3
from sqlite3 import Error
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from PIL.Image import core as Image
from scipy.stats import gaussian_kde

def createConnection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)

    return conn

def selectAllPositions(conn, simulationID):
    cur = conn.cursor()
    #cur.execute("SELECT star.id,mass,position.timestep,position.x,position.y,position.z FROM star INNER JOIN velocity on velocity.id_star = star.id INNER JOIN position on position.id_star = star.id where star.id_simulation = ?1 AND position.timestep = velocity.timestep order by position.timestep", (simulationID,))
    cur.execute("SELECT star.id,mass,position.timestep,position.x,position.y,position.z FROM star INNER JOIN position on position.id_star = star.id where star.id_simulation = ?1 order by position.timestep ", (simulationID,)) #LIMIT 10000000 OFFSET 10000000
    rows = np.array(cur.fetchall())
    return rows

def plot_star_series(output,data):
    #loop timesteps
    for i in np.unique(data[:,2]):
        print(i)
        plotData = data[data[:,2] == i]
        name = output+'\starPositions'+str(int(i))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(plotData[:,3], plotData[:,4], plotData[:,5], s=50,c='red')
        plt.savefig(name+'.jpg')
        plt.close(fig)

#returns center of mass and max distance
def plotDimensions(data):
    com =  np.zeros(3)
    com[0] = np.mean(data[:,3])
    com[1] = np.mean(data[:,4])
    com[2] = np.mean(data[:,5])
    maxDist =  np.zeros(3)
    for x in data[:,3]:
        if abs(x-com[0])>maxDist[0]:
            maxDist[0]=abs(x-com[0])
    for y in data[:,4]:
        if abs(y-com[1])>maxDist[1]:
            maxDist[1]=abs(y-com[1])
    for z in data[:,5]:
        if abs(z-com[2])>maxDist[2]:
            maxDist[2]=abs(z-com[2])
    return com,maxDist

def toCylinder(x,y,z):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi, z

def plot2Dxy(output,data):
    #loop timesteps
    for i in np.unique(data[:,2]):
        timestepData = data[data[:,2] == i]
        com,maxDist = plotDimensions(timestepData)
        fig = plt.figure()
        plt.plot(timestepData[:,3], timestepData[:,4], 'ro')
        plt.xlim(com[0]-maxDist[0])
        plt.ylim(com[1]-maxDist[1])
        name = output+'\starPositions'+str(int(i))
        plt.savefig(name+'.jpg')
        plt.close(fig)

def plot2DCylinder(output,data):
    #loop timesteps
    for i in np.unique(data[:,2]):
        #if i>1000:
        timestepData = data[data[:,2] == i]
        com,maxDist = plotDimensions(timestepData)
        fig = plt.figure(figsize=(1,3),dpi=100)
        fig.set_size_inches(27,9)
        fig.subplots_adjust(hspace=0.4, wspace=0.4)

        ax = fig.add_subplot(1, 3, 3)
        ax.tick_params(axis='both', labelsize=20)
        ax.set_xlabel('x [pc]', fontsize=20)
        ax.set_ylabel('y [pc]', fontsize=20)
        xy = np.vstack([timestepData[:,3],timestepData[:,4]])
        z1 = gaussian_kde(xy)(xy)
        # Sort the points by density, so that the densest points are plotted last
        idx = z1.argsort()
        x1, y1, z1 = timestepData[:,3][idx], timestepData[:,4][idx], z1[idx]
        x_mean = np.mean(x1)
        #x_std = np.std(x1)
        y_mean = np.mean(y1)
        #y_std = np.std(y1)
        ax.scatter(x1, y1, c=z1, s=5, edgecolor='')
        #ax.set_xlim([x_mean-50, x_mean+50])
        #ax.set_ylim([y_mean-50, y_mean+50])
        #ax.set_xlim([x_mean-x_std, x_mean+50])
        #ax.set_ylim([y_mean-50, y_mean+50])

        x, y, z = toCylinder(timestepData[:,3], timestepData[:,4],timestepData[:,5])

        ax = fig.add_subplot(1, 3, 1)
        ax.set_xlabel('r [pc]', fontsize=20)
        ax.set_ylabel('z [pc]', fontsize=20)
        ax.tick_params(axis='both', labelsize=20)
        xz = np.vstack([x,z])
        y1 = gaussian_kde(xz)(xz)
        # Sort the points by density, so that the densest points are plotted last
        idx = y1.argsort()
        x1, y1, z1 = x[idx], y1[idx], z[idx]
        ax.scatter(x1, z1, c=y1, s=5, edgecolor='')
        x_mean = np.mean(x1)
        x_std = np.std(x1)
        z_mean = np.mean(z1)
        #z_std = np.std(z1)
        ax.set_xlim([x_mean-x_std*3, x_mean+x_std*3])
        ax.set_ylim([z_mean-12, z_mean+12])

        ax = fig.add_subplot(1, 3, 2)
        ax.set_xlabel(r'$\phi$ [pc]', fontsize=20)
        ax.set_ylabel('z [pc]', fontsize=20)
        ax.tick_params(axis='both', labelsize=20)
        yz = np.vstack([y,z])
        x1 = gaussian_kde(yz)(yz)
        # Sort the points by density, so that the densest points are plotted last
        idx = x1.argsort()
        x1, y1, z1 = x1[idx], y[idx], z[idx]
        ax.scatter(y1, z1, c=x1, s=5, edgecolor='')
        y_mean = np.mean(y1)
        y_std = np.std(y1)
        z_mean = np.mean(z1)
        #z_std = np.std(z1)
        ax.set_xlim([y_mean-y_std*3, y_mean+y_std*3])
        ax.set_ylim([z_mean-12, z_mean+12])


        name = output+'\starPositions'+str(int(i))
        plt.savefig(name+'.jpg', dpi=100)
        plt.close(fig)

def main():
    simulationID = 1
    database = r"E:\Master_Thesis\VS_Project\N_Body\Output\Database\Default.db"
    output = r"E:\Master_Thesis\VS_Project\N_Body\Output\Simulation2.5k"# + str(simulationID)
    # create a database connection
    conn = createConnection(database)
    with conn:
        data = selectAllPositions(conn, simulationID)
        #plot_star_series(output, data)
        plot2DCylinder(output, data)

if __name__ == '__main__':
    main()
