import sqlite3
from sqlite3 import Error
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from PIL.Image import core as Image

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
    cur.execute("SELECT star.id,mass,position.timestep,position.x,position.y,position.z,velocity.x,velocity.y,velocity.z FROM star INNER JOIN velocity on velocity.id_star = star.id INNER JOIN position on position.id_star = star.id where star.id_simulation = ?1 AND position.timestep = velocity.timestep order by position.timestep", (simulationID,))

    rows = cur.fetchall()

    return rows

def plot_star_series(output,data):
    for i in np.unique(data[:,2]):
        print(i)
        plotData = data[data[:,2] == i]
        name = output+'\starPositions'+str(int(i))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(plotData[:,3], plotData[:,4], plotData[:,5], s=50,c='red')
        plt.savefig(name+'.jpg')
        plt.close(fig)

def main():
    simulationID = 1
    database = r"E:\Master_Thesis\VS_Project\N_Body\Output\Database\Default.db"
    output = r"E:\Master_Thesis\VS_Project\N_Body\Output\Simulation" + str(simulationID)
    # create a database connection
    conn = createConnection(database)
    with conn:
        data = selectAllPositions(conn, simulationID)
        plot_star_series(output, np.array(data))

if __name__ == '__main__':
    main()
