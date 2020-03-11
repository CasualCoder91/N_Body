#ifndef DATABASE_H
#define DATABASE_H

#include <iostream>
#include <vector>

#include "sqlite3.h"
#include "Star.h"
#include "Analysis.h"
#include "Simulation.h"

using Record = std::vector<std::string>;
using Records = std::vector<Record>;

class Database{
	sqlite3* db;
	const char* dBName = "Data/test.db";
	bool isOpen;

private:

public:
	Database();
	bool open(const char* name ="");
	bool exec(char* sql);
	void setup();
	int getLastID();
	int selectLastID(std::string table);
	int insert(Parameters& parameters);
	void insertStars(int simulationID, std::vector<Star*>& stars, int timestep=0);
	int insertAnalysis(int simulationID, Analysis analysis);
	void insertAnalysisdtEnergy(int analysisID,int dt, double kinE, double potE);
	void insertAnalysisdtVelocity(int analysisID, int dt, double velocity);
	void timestep(int timestep, std::vector<Star*>& stars);
	void insertStar(int simulationID, Star* star, int& timestep);
	void insertPosition(int& idStar, Vec3D& position, int& timestep);
	void insertVelocity(int& idStar, Vec3D& velocity, int& timestep);
	std::vector<Simulation> selectSimulations();
	std::vector<Vec3D> selectVelocities(int timestep);
	std::vector<int> selectTimesteps();
	std::vector<Star*> selectStars(int simulationID, int timeStep);
};

#endif