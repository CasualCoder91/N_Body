#pragma once

#include <string>
#include <filesystem>
#include <vector>

class Plot {

	std::string dataDirectory;
	std::string plotDirectory;
	std::string pythonFile;
	bool showPlot;

public:
	Plot();
	Plot(std::string dataDirectory, std::string plotDirectory, bool showPlot);

	void plot(std::string function,std::vector<std::string>params);
};