/**
 * SQLite Helper. Link between Database and Objects like Stars, Analysis, Simulation.
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
#include "Analysis.h"
#include "SimulationData.h"

using Record = std::vector<std::string>;
using Records = std::vector<Record>;

class Database{
	sqlite3* db;
	static char* dataBaseDataPath;
	bool isOpen;

private:

public:
	Database();
	/**
	@brief Opens a database connection. In case the database does not exist yet it creates a new one 
	@param name file path of the database. If not passed \ref Database.dataBaseDataPath is used. If name is given, that database will be used/created.
	*/
	bool open(char* name ="");
	bool exec(char* sql);
	void setup();
	int getLastID();
	int selectLastID(std::string table);
	int insert(Parameters* parameters);
	void insertStars(int simulationID, std::vector<Star*>& stars, int timestep=0);
	int insertAnalysis(int simulationID, Analysis analysis);
	void insertAnalysisdtEnergy(int analysisID,int dt, double kinE, double potE);
	void insertAnalysisdtVelocity(int analysisID, int dt, double velocity);
	void timestep(int timestep, std::vector<Star*>& stars);
	void insertStar(int simulationID, Star* star, int& timestep);
	void insertPosition(int& idStar, Vec3D& position, int& timestep);
	void insertVelocity(int& idStar, Vec3D& velocity, int& timestep);
	std::vector<SimulationData> selectSimulationData(int ID = -1);
	std::vector<Vec3D> selectVelocities(int timestep);
	std::vector<int> selectTimesteps();
	std::vector<Star*> selectStars(int simulationID, int timeStep);
};

#endif