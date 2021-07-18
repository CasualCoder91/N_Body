import sqlite3
import numpy as np
from point import Point

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

    def insert_position_HTP(self,timestep,id_star,ascension,declination,commit=True):
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        cur.execute("INSERT INTO position (id_star,timestep,aHTP,dHTP) VALUES (?1,?2,?3,?4)",
                   (id_star,timestep,ascension,declination))
        if commit:
            self.conn.commit()

    def update_position_HTP(self,timestep,id_star,ascension,declination,commit=True):
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        cur.execute("UPDATE position set aHTP=?1, dHTP=?2 where timestep = ?3 and id_star = ?4",
                   (ascension,declination,timestep,id_star))
        if commit:
            self.conn.commit()

    def insert_velocity_HTP(self,timestep,id_star,ascension,declination,commit=True):
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        cur.execute("INSERT INTO velocity (id_star,timestep,aHTP,dHTP) VALUES (?1,?2,?3,?4)",
                   (id_star,timestep,ascension,declination))
        if commit:
            self.conn.commit()

    def insert_velocities_HTP(self,timestep,points):
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        for point in points:
            cur.execute("INSERT INTO velocity (id_star,timestep,aHTP,dHTP) VALUES (?1,?2,?3,?4)",
                        (point.id, timestep, point.velocity[0], point.velocity[1]))
        self.conn.commit()

    def update_velocity_HTP(self,timestep,id_star,ascension,declination):
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        cur.execute("UPDATE velocity set aHTP=?1, dHTP=?2 where timestep = ?3 and id_star = ?4",
                   (ascension,declination,timestep,id_star))
        self.conn.commit()

    def update_velocities_HTP(self,timestep,points):
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        for point in points:
            cur.execute("UPDATE velocity set aHTP=?1, dHTP=?2 where timestep = ?3 and id_star = ?4",
                       (point.velocity[0],point.velocity[1],timestep,point.id))
        self.conn.commit()


    def insert_star(self,timestep,magnitude,ascension,declination,v_ascension=np.nan,v_declination=np.nan,observed=True,commit=True):
        if self.conn is None:
            connect(database_path)
        cur = self.conn.cursor()
        cur.execute("INSERT INTO star (id_simulation, magnitude, isObserved) VALUES (?1,?2,?3)",
            (simulation_id,magnitude,observed))
        id_star = cur.lastrowid
        if np.isfinite(ascension) and np.isfinite(declination):
            self.insert_position_HTP(timestep,id_star,ascension,declination,commit)
        if np.isfinite(v_ascension) and np.isfinite(v_declination):
            self.insert_velocity_HTP(timestep,id_star,v_ascension,v_declination,commit)

        return id_star

    def update_point(self, point, timestep, commit=True):
        if self.conn is None:
            connect(database_path)
        if point.id < 1:
            print("update_point: id has to be greater than 0!")
            return
        if np.isfinite(point.position[0]) and np.isfinite(point.position[1]):
            self.insert_position_HTP(timestep,point.id,point.position[0],point.position[1],False)
        if np.isfinite(point.velocity[0]) and np.isfinite(point.velocity[1]):
            self.insert_velocity_HTP(timestep,point.id,point.velocity[0],point.velocity[1],False)
        if commit:
            self.conn.commit()


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

    def select_points(self, timestep, observed = True):
        cur = self.conn.cursor()
        cur.execute("""SELECT star.id, position.timestep, position.aHTP, position.dHTP, velocity.aHTP, velocity.dHTP, star.isCluster, star.magnitude 
		    FROM star 
		    LEFT JOIN position on position.id_star = star.id 
		    LEFT JOIN velocity on velocity.id_star = star.id AND position.timestep = velocity.timestep 
		    WHERE star.id_simulation = ?1 AND position.timestep = ?2 
		    AND star.isObserved = ?3""", (simulation_id,timestep,observed))
        result = cur.fetchall()
        array = np.ndarray((len(result),),dtype=object)
        for i, line in enumerate(result):
            array[i] = Point(id=line[0],position=line[2:4],velocity=line[4:6],magnitude=line[7])
        return array

    def insert_points(self,points,timestep):
        for point in points:
            self.insert_star(timestep,point.magnitude,point.position[0],point.position[1],
                             point.velocity[0],point.velocity[1],False)
        self.conn.commit()

    def update_points(self,points,timestep):
        for point in points:
            self.update_point(point,timestep,False)
        self.conn.commit()


