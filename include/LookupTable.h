#pragma once

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "Parameters.h"

class LookupTable {

	std::map<double, double> map;
	std::string fileName;
	std::string filePath = "src/LookupTables/";
	std::string delimiter = ",";

public:
	double get(double key);
	void setMap(std::vector<double> keys, std::vector<double> values);
	void makeFile(std::string fileName = "", std::string header = "");
	LookupTable(std::string filename, std::string delimiter = ",");
private:
	void init();
	bool checkIsDouble(std::string inputString);
};