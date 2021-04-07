import sqlite3
from sqlite3 import Error
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from PIL.Image import core as Image
from scipy.stats import gaussian_kde
import scipy.spatial.distance
from scipy import stats
from mayavi import mlab
import matplotlib.cm as cm

# Functions from @Mateen Ulhaq and @karlo
def set_axes_equal(ax: plt.Axes):
    """Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)

def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

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
    cur.execute("""SELECT star.id,mass,position.timestep,position.x,position.y,position.z
        FROM star
        INNER JOIN position on position.id_star = star.id
        where star.id_simulation = ?1
        order by position.timestep""", (simulationID,)) #and star.isCluster=0 LIMIT 10000000 OFFSET 10000000 and star.isCluster=1
    rows = np.array(cur.fetchall())
    return rows

def select2DPositions(conn, simulationID):
    cur = conn.cursor()
    cur.execute("""SELECT star.id,mass,position.timestep,position.aHEQ,position.dHEQ,star.isCluster FROM star
       INNER JOIN position on position.id_star = star.id
       where star.id_simulation = ?1 And position.timestep<3 order by position.timestep""", (simulationID,)) #and star.isCluster=0 LIMIT 10000000 OFFSET 10000000 and star.isCluster=1
    rows = np.array(cur.fetchall())
    return rows

def plot_star_series(output,data,show=False):
    #loop timesteps
    for i in np.unique(data[:,2]):
        print(i)
        plotData = data[data[:,2] == i]
        name = output+'\starPositions'+str(int(i))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(plotData[:,3], plotData[:,4], plotData[:,5], s=1,c='red')
        if(show):
            plt.show()
        plt.savefig(name+'.jpg')
        plt.close(fig)

def plotDensity(output,data,show=False):
    plotData = data[data[:,2] == 0]
    x = plotData[:,3]
    y = plotData[:,4]
    z = plotData[:,5]

    xyz = np.vstack([x,y,z])
    kde = stats.gaussian_kde(xyz)
    density = kde(xyz)

    idx = density.argsort()
    x, y, z, density = x[idx], y[idx], z[idx], density[idx]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_box_aspect([1,1,1])
    cax = ax.scatter(x, y, z, c=density, s=0.1, cmap=cm.hot)
    fig.colorbar(cax)
    # ax.set_proj_type('ortho') # OPTIONAL - default is perspective (shown in image above)
    set_axes_equal(ax) # IMPORTANT - this is also required
    plt.show()

    # Plot scatter with mayavi
    #figure = mlab.figure('DensityPlot')
    #pts = mlab.points3d(x, y, z, density, scale_mode='none', scale_factor=0.01)
    #print(pts)
    #mlab.axes()
    #mlab.show()
    #imgmap = mlab.screenshot(mode='rgba', antialiased=True)
    #plt.axis('off')
    #plt.imsave(arr=imgmap, fname='rbf_concept_3d.pdf')

#returns center of mass and max distance
def plotDimensions(data):
    com =  np.zeros(3)
    maxDist =  np.zeros(3)
    com[0] = np.mean(data[:,3])
    com[1] = np.mean(data[:,4])
    print(np.size(data, 1))
    if(np.size(data, 1)>5):
        com[2] = np.mean(data[:,5])
        for z in data[:,5]:
            if abs(z-com[2])>maxDist[2]:
                maxDist[2]=abs(z-com[2])
    for x in data[:,3]:
        if abs(x-com[0])>maxDist[0]:
            maxDist[0]=abs(x-com[0])
    for y in data[:,4]:
        if abs(y-com[1])>maxDist[1]:
            maxDist[1]=abs(y-com[1])
    return com,maxDist

def toCylinder(x,y,z):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi, z

def plot2Dxy(output,data):
    #loop timesteps
    #for i in np.unique(data[:,2]):
    timestepData = data[data[:,2] == 0]
    #com,maxDist = plotDimensions(timestepData)
    fig = plt.figure()

    fig.set_size_inches(9,9)
    #print(maxDist[0],maxDist[1])
    #plotDist = np.minimum(maxDist[0],maxDist[1])
    #plt.xlim(com[0]-plotDist, com[0]+plotDist)
    #plt.ylim(com[1]-plotDist, com[1]+plotDist)
    colors = np.where(timestepData[:,5]==1,'y','k')
    plt.scatter(timestepData[:,3], timestepData[:,4], c=colors)
    #print(plotDist)
    #name = output+'\starPositions'+str(int(i))
    plt.show()
    #plt.savefig(name+'.jpg')
    #plt.close(fig)

def plotProjection(output,data):
    #loop timesteps
    origin = np.array([8300, 0, 0])
    focus = np.array([9594,-640,-52])-origin
    focus = focus / np.linalg.norm(focus)
    fov = 10*0.0174533
    for i in np.unique(data[:,2]):
        timestepData = data[data[:,2] == i]
        timestepData = timestepData[:,3:]
        #com,maxDist = plotDimensions(timestepData)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect(1)
        for position in timestepData:
            position = position-origin
            length = np.linalg.norm(position)
            position_u = position / length
            angleXY = np.arctan2(focus[1],focus[0]) - np.arctan2(position_u[1],position_u[0])
            #if angleXY > np.pi:
            #    angleXY = angleXY - 2*np.pi
            #if angleXY < -np.pi:
            #    angleXY = angleXY + 2*np.pi
            #lengthXY = np.linalg.norm(position[:2])
            angleXZ = np.arcsin(focus[2]) - np.arcsin(position_u[2])
            #lengthXZ = np.abs(np.cross(position, focus)*focus)[2]
            #print("lengthXZ",lengthXZ)
            #print("position",position)
            #print("angleXY",angleXY)
            #print("lengthXY",lengthXY)
            #print("angleXZ",angleXZ)
            #print("x",lengthXY*np.sin(angleXY))
            #print("y",position[2]*np.sin(angleXZ))
            plt.scatter(length*np.tan(angleXY), length*np.tan(angleXZ),c='red',s=2)
            #plt.xlim(com[0]-maxDist[0])
            #plt.ylim(com[1]-maxDist[1])
        name = output+'\starPositions'+str(int(i))
        plt.show()
        plt.savefig(name+'.jpg')
        plt.close(fig)


def plot2DCylinder(output,data):
    #loop timesteps
    for i in np.unique(data[:,2]):
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
        ax.scatter(x1, y1, c=z1, s=10, edgecolor='')
        #ax.set_xlim([x_mean-50, x_mean+50])
        #ax.set_ylim([y_mean-50, y_mean+50])
        #ax.set_xlim([x_mean-x_std, x_mean+50])
        #ax.set_ylim([y_mean-50, y_mean+50])

        x, y, z = toCylinder(timestepData[:,3], timestepData[:,4],timestepData[:,5])
        y = np.rad2deg(y)

        ax = fig.add_subplot(1, 3, 1)
        ax.set_xlabel('r [pc]', fontsize=20)
        ax.set_ylabel('z [pc]', fontsize=20)
        ax.tick_params(axis='both', labelsize=20)
        xz = np.vstack([x,z])
        y1 = gaussian_kde(xz)(xz)
        # Sort the points by density, so that the densest points are plotted last
        idx = y1.argsort()
        x1, y1, z1 = x[idx], y1[idx], z[idx]
        ax.scatter(x1, z1, c=y1, s=10, edgecolor='')
        x_mean = np.mean(x1)
        x_std = np.std(x1)
        z_mean = np.mean(z1)
        #z_std = np.std(z1)
        ax.set_xlim([x_mean-x_std*3, x_mean+x_std*3])
        ax.set_ylim([z_mean-12, z_mean+12])

        ax = fig.add_subplot(1, 3, 2)
        ax.set_xlabel(r'$\phi$ [°]', fontsize=20)
        ax.set_ylabel('r [pc]', fontsize=20)
        ax.tick_params(axis='both', labelsize=20)
        yx = np.vstack([y,x])
        z1 = gaussian_kde(yx)(yx)
        # Sort the points by density, so that the densest points are plotted last
        idx = z1.argsort()
        x1, y1, z1 = x[idx], y[idx], z1[idx]
        ax.scatter(y1, x1, c=z1, s=10, edgecolor='')
        y_mean = np.mean(y1)
        y_std = np.std(y1)
        x_mean = np.mean(x1)
        x_std = np.std(x1)
        ax.set_xlim([y_mean-y_std*3, y_mean+y_std*3])
        ax.set_ylim([x_mean-x_std*3, x_mean+x_std*3])


        name = output+'\starPositions'+str(int(i))
        plt.savefig(name+'.jpg', dpi=100)
        plt.close(fig)

def main():
    simulationID = 1
    database = r"E:\Master_Thesis\VS_Project\N_Body\Output\Database\default.db"
    output = r"E:\Master_Thesis\VS_Project\N_Body\Output\NGC2244"# + str(simulationID)
    # create a database connection
    conn = createConnection(database)
    with conn:
        data = select2DPositions(conn, simulationID)
        plot2Dxy(output,data)
        #data = selectAllPositions(conn, simulationID)
        #plotDensity(output, data,True)
        #plotProjection(output, data)

if __name__ == '__main__':
    main()
