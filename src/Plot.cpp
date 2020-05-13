#include "..\include\Plot.h"

Plot::Plot(){
	dataDirectory = "/docs/pyplots/";
	plotDirectory = "/Output/";
	pythonFile = std::filesystem::current_path().string() + "/src/Plots/plotCentral.py";
	showPlot = false;
}

Plot::Plot(std::string dataDirectory, std::string plotDirectory, bool showPlot){
	this->dataDirectory = dataDirectory;
	this->plotDirectory = plotDirectory;
	this->showPlot = showPlot;
	this->pythonFile = std::filesystem::current_path().string() + "/src/Plots/plotCentral.py";
}

void Plot::plot(std::string function, std::vector<std::string>params){
	std::string command = "python "+ pythonFile + " " + dataDirectory + " "+ (showPlot ? "true" : "false") +" " + function + " " + plotDirectory;
	for (std::string param : params) {
		command += " " + param;
	}
	system(command.c_str());
}
