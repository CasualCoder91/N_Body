import sqlite3
import pandas as pd
from sqlite3 import Error
from st_dbscan import ST_DBSCAN
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from PIL.Image import core as Image
from scipy.stats import gaussian_kde
import scipy.spatial.distance
from scipy import spatial #cKDTree

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
    cur.execute("""SELECT star.id,mass,position.timestep,position2D.x,position2D.y FROM star
       INNER JOIN position on position.id_star = star.id
       INNER JOIN position2D on position.id=position2D.fk_position
       where star.id_simulation = ?1 order by position.timestep""", (simulationID,)) #and star.isCluster=0 LIMIT 10000000 OFFSET 10000000 and star.isCluster=1
    rows = np.array(cur.fetchall())
    return rows

def select2D(conn, simulationID):
    cur = conn.cursor()
    cur.execute("""SELECT simulation.dt, star.id,mass,position.timestep,position2D.x,position2D.y,velocity2D.x,velocity2D.y
                    FROM simulation
                    INNER JOIN star on simulation.id=star.id_simulation
                    INNER JOIN position on position.id_star = star.id
                    INNER JOIN position2D on position.id=position2D.fk_position
                    INNER JOIN velocity on velocity.id_star = star.id AND position.timestep = velocity.timestep
                    INNER JOIN velocity2D on velocity2D.fk_velocity = velocity.id
                    WHERE star.id_simulation = ?1 AND star.isCluster=1 AND position.timestep<3 order by position.timestep""", (simulationID,)) # LIMIT 10000000 OFFSET 10000000 and star.isCluster=1
    rows = np.array(cur.fetchall())
    return rows

def plot(data, labels):
    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']

    for i in range(-1, len(set(labels))):
        if i == -1:
            col = [0, 0, 0, 1]
        else:
            col = colors[i % len(colors)]

        clust = data[np.where(labels==i)]
        plt.scatter(clust[:,0], clust[:,1], c=[col], s=1)
    plt.show()

    return None

def main():
    simulationID = 1
    database = r"E:\Master_Thesis\VS_Project\N_Body\Output\Database\default.db"
    output = r"E:\Master_Thesis\VS_Project\N_Body\Output\NGC2244"# + str(simulationID)
    # create a database connection
    con = createConnection(database)
    with con:
        df = pd.read_sql_query("""SELECT star.id,mass,position.timestep,position2D.x,position2D.y FROM star
           INNER JOIN position on position.id_star = star.id
           INNER JOIN position2D on position.id=position2D.fk_position
           where star.id_simulation = 1 and position.timestep<10 order by position.timestep""", con) #
        print(df.head())
    con.close()
    #df['x'] = (df['x'] - df['x'].min()) / (df['x'].max() - df['x'].min())
    #df['y'] = (df['y'] - df['y'].min()) / (df['y'].max() - df['y'].min())
    # transform to numpy array
    data = df.loc[:, ['timestep','x','y']].values
    print(len(data))
    st_dbscan = ST_DBSCAN(eps1 = 0.005, eps2 = 5, min_samples = 100)
    st_dbscan.fit(data)
    #st_dbscan.fit_frame_split(data, frame_size = 3)
    #st_dbscan.fit(data)
    print(st_dbscan.labels)
    plot(data[:,1:], st_dbscan.labels)


if __name__ == '__main__':
    main()



#data = select2D(conn, simulationID)
#t0 = data[data[:,3] == 0] # data at timestep 0
#t0 = t0[:, [4,5]] # only x and y
#t1 = data[data[:,3] == 1] # data at timestep 1
#t1 = t1[:, [4,5]] # x, y at timestep 1
#tree = spatial.cKDTree(t1)
#mindist, minid = tree.query(t0)
#print(mindist)
#print(minid)
