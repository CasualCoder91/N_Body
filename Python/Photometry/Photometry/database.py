import sqlite3
import numpy as np

from config import simulation_id, database_path

class Database:
    conn = None #DB connection

    def connect(self,db_path):
        try:
            self.conn = sqlite3.connect(db_path)
        except:
            print("Connection failed!")

    def __init__(self):
        self.connect(database_path)

    def select_last_ID(self,table):
        """ 
        Gets highest id in a table

        :param table: name of the table
        :type table: str
        :returns: highest id in the given table
        :rtype: int
        """
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        cur.execute("select id from "+table+" order by id desc LIMIT 1")
        id, = cur.fetchone()
        return id

    def insert_position_HTP(self,timestep,id_star,ascension,declination):
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        cur.execute("INSERT INTO position (id_star,timestep,aHTP,dHTP) VALUES (?1,?2,?3,?4)",
                   (id_star,timestep,ascension,declination))
        self.conn.commit()

    def insert_star(self,timestep,magnitude,ascension,declination,is_observed=1):
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        cur.execute("INSERT INTO star (id_simulation, magnitude, isObserved) VALUES (?1,?2,?3)",
            (simulation_id,magnitude,is_observed))
        id_star = cur.lastrowid
        self.insert_position_HTP(timestep,id_star,ascension,declination)
        return id_star

    def select_2d_stars(self,timestep):
        cur = self.conn.cursor()
        cur.execute("""SELECT position.rH, position.aHTP, position.dHTP, star.mass 
           FROM star
           INNER JOIN position on position.id_star = star.id
           where star.id_simulation = ?1
           and position.timestep = ?2
           and position.rH NOT NULL""", (simulation_id,timestep))
        rows = np.array(cur.fetchall())
        return rows
