/**
 * SQLite Helper. Link between database and objects like stars, analysis, simulation, ...
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#ifndef DATABASE_H
#define DATABASE_H

#include <iostream>
#include <vector>

#include "sqlite3.h"
#include "Star.h"
#include "Projection.h"
#include "Point.h"
#include "Vec2D.h"
#include "Misc.h"
//#include "Constants.h"

class Database{
	const char* default_database_path = "Output/Database/Default.db";
	const char* database_path = "Output/Database/Default.db";
	bool isOpen;

private:

	/**
	@brief Executes sql statement.
	@param sql The statement which shall be executed.
	*/
	bool exec(char* sql);

	int getLastID();
public:
	sqlite3* db;
	Database(const std::string& path);
	~Database();
	/**
	@brief Opens a database connection. In case the database does not exist yet it creates a new one 
	*/
	bool open();
	bool close();
	/** @brief Creates all tables (and database) if they do not exist yet. */
	void setup();
	/** @brief returns the largest id from the given \p table */
	int selectLastID(std::string table);
	/** @brief inserts the given parameters into the "simulation" table and returns the id of the new simulation*/
	int Database::insertSimulation();
	/** @brief inserts the stars (including positions and velocities)*/
	void insertStars(int simulationID, const std::vector<Star>& stars, int timestep=0, bool clusterStars=true);

	void delete_star(int simulation_id, int star_id);

	void delete_stars(const int simulation_id, const std::vector<int> stars_to_delete);

	/** @brief inserts the analysis parameters for the given simulation*/
	int insertAnalysis(int simulationID);
	/** @brief inserts (or replaces/updates) one record of kinetic, potential and total energy */
	void insertAnalysisdtEnergy(int simulation_id,int dt, double kinE, double potE);
	/** @brief inserts (or replaces/updates) one record of the average velocity */
	void insertAnalysisdtVelocity3D(int simulation_id, int dt, double avgVel3DCluster, double disp3DCluster, double avgVel3DFS, double disp3DFS);
	void insertAnalysisdtVelocity2D(int simulation_id, int dt, double avgVel2DCluster, double disp2DCluster, double avgVel2DFS, double disp2DFS);

	void print_analysis_dt_velocity_2D();

	void insert_analysis_dt_min_dist(int simulation_id, int dt, double minimum_distance);

	double select_analysis_dt_min_dist(int simulation_id, int dt);

	/** @brief inserts positions and velocities of given \p stars and \p timestep */
	void timestep(int timestep, const std::vector<Star>& stars);

	void generateHEQ(int simulationID);

	//HTP = Heliocentric Telescope Polar (LookAt Vector is primary axis)
	void generateHTP(int simulationID, bool observed = false);

	void generateMagnitude(int simulationID, bool observed = false);

	void insertPowerLaw(int simulationID, std::vector<double> massLimits, std::vector<double> exponents);
	/** 
	@brief inserts one star (including positions and velocities)
	@note for multiple stars use insertStars instead!
	*/
	void insertStar(int simulationID, Star* star, int& timestep, bool clusterStar=true);
	/**
	@brief inserts the position of one star
	@note extensive use (loop) is not recommended
	*/
	void insertPosition(int& idStar, Vec3D& position, int& timestep);
	/**
	@brief inserts the velocity of one star
	@note extensive use (loop) is not recommended
	*/
	void insertVelocity(int& idStar, Vec3D& velocity, int& timestep);
	/** @brief returns all saved simulations data by default. specify \p ID to retrieve a specific record */
	//std::vector<SimulationData> selectSimulationData(int ID = -1);

	void selectConstants(int ID);

	bool printSimulations();

	/** @brief returns all velocities at the given \p timestep */
	std::vector<Vec3D> selectVelocities3D(int simulationID, int timestep = -1, bool fieldStars= false, bool clusterStars=true);
	std::vector<Vec2D> selectVelocitiesHTP(int simulationID, int timestep = -1, bool fieldStars = false, bool clusterStars = true, bool observed = false);
	/** @brief returns all timesteps for simulation with ID \p simulationID*/
	std::vector<int> selectTimesteps(int simulationID);


	Star select_star(int star_id, int timestep = 0);

	/** @brief returns all stars for a given simulation with velocity and position at the given timestep (pass 0 to retrieve initial values) 
		@param simulationID ID of the Simulation for which to retrieve the stars
		@param timeStep position and velocity of the star depend on the timestep
		@param observed set true to retrieve stars observed with photutils, false to get simulated stars
		@param type: -1 any type, 0: field, 1: cluster
	*/
	std::vector<Star> select_stars(int simulationID, int timeStep, bool observed = false, int type = -1);
	/** @brief saves all stars at all timesteps into a file. Passed \p filePath must exist and is relative to the executable */
	void outputStars(int simulationID, std::string filePath, bool allStars = true, bool clusterStars = false, bool fieldStars = false);

	//pass timestep = -1 to output all timestep (one file per step)
	void outputStars2D(int simulationID, std::string filePath, int timestep = -1);

	Point select_point(int point_id, int simulationID = 1, int timeStep = 0);

	std::vector<Point>select_points(int simulationID = 1, int timeStep = 0, int observed = 0);

	std::vector<std::vector<Point>>select_time_series_points(int simulationID = 1, int timeStep=0, int nTimeSteps = 2, bool observed = false);

	//updates velocity.aHTP and velocity.dHTP
	void update_points(std::vector<Point>& points, int timestep = 0);

	void set_mass(const std::vector<Point>& points);

	void set_fk_star(std::vector<Point>& points);

	void set_extinction(std::vector<Star>& stars);

	void print_clustering_info(int simulation_id);

	void print_simulated_cluster_info(int simulation_id);

	double confidence_score(int simulation_id);

	void set_mapped_star_mass();

};

#endif