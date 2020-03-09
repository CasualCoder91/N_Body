#pragma once

#include <iostream>
#include <vector>

#include "sqlite3.h"
#include "Star.h"
#include "Parameters.h"

class Database{
	sqlite3* db;
	const char* dBName = "Data/test.db";

public:
	Database();
	bool open(const char* name ="");
	bool exec(char* sql);
	void setup();
	int insert(Parameters parameters);
	void insert(int simulationID, std::vector<Star*>& stars);

};

