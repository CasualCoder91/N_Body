/**
 * Used to create, store and retrieve lookup tables. Lookup tables are very usefull (absolutely necessary) to speedup runtime.
 * Expensive calculations (numerical integrations) are saved into tables. Due to obvious limitations only discrete values can be stored, but with enough granularity this procedure becomes justifiable.
 * Currently the bulge velocity dispersion is stored in such a table since a double integral is needed to calculate it.
 * The lookuptables are stored i the subdirectory src/LookupTables/
 *
 * @author Alarich Herzner
 * @version 1.0 15.05.2020
*/

#pragma once

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

extern bool debug;

class LookupTable {

	std::map<double, double> map;
	std::string fileName;
	std::string filePath = "src/LookupTables/";
	std::string delimiter = ",";

public:
	/**@brief returns whether or not the lookup table is currently empty*/
	bool isEmpty();
	/**@brief returns the value for the passed \p key (example: key=radius value=dispersion)*/
	double get(double key);
	/**@brief sets keys and values of the table. Size of the vectors \p keys and \p values must be equal*/
	void setMap(std::vector<double> keys, std::vector<double> values);
	/**@brief creates a file containing keys, values and (if specified) a \p header. If no \p fileName is specified, the one given at construction is used.*/
	void makeFile(std::string fileName = "", std::string header = "");
	/**@brief creates a lookup table object reading its content from the \p filename*/
	LookupTable(std::string filename, std::string delimiter = ",");
private:
	void init();
	bool checkIsDouble(std::string inputString);
};