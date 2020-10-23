import sqlite3
import numpy as np
import matplotlib.pyplot as plt

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

def select2dPositions(conn, simulationID):
    cur = conn.cursor()
    cur.execute("""SELECT position.timestep,position2D.x,position2D.y,star.mass FROM star
       INNER JOIN position on position.id_star = star.id
       INNER JOIN position2D on position.id=position2D.fk_position
       where star.id_simulation = ?1 order by position.timestep""", (simulationID,)) #and star.isCluster=0 LIMIT 10000000 OFFSET 10000000 and star.isCluster=1
    rows = np.array(cur.fetchall())
    return rows

def main():
    simulationID = 1
    database = r"E:\Master_Thesis\VS_Project\N_Body\Output\Database\NGC2244.db"
    output = r"E:\Master_Thesis\VS_Project\N_Body\Output\NGC2244"# + str(simulationID)
    # create a database connection
    conn = createConnection(database)
    with conn:
        data = select2dPositions(conn, simulationID)
        for i in np.unique(data[:,0]):
            timestepData = data[data[:,0] == i]
            #plt.scatter(timestepData[:,1], timestepData[:,2])
            #plt.show()
            if(i<9 or (i%4==0 and i<52) or (i%13==0 and i>52)):
                np.savetxt(output+r"\NGC2244_pos_"+str(int(i))+'.dat', timestepData[:,1:], delimiter=',')   # X is an array

if __name__ == '__main__':
    main()
