#include "SimulationData.h"


SimulationData::SimulationData(){
	this->filePath = "./simulation.cfg";
	config.Load(filePath);
	if (!config.Get("n_Stars", n_Stars)) {
		std::cout << "n_Stars missing in " + filePath << std::endl;
	}
	if (!config.Get("boxLength", boxLength)) {
		std::cout << "boxLength missing in " + filePath << std::endl;
	}
	if (!config.Get("dt", dt)) {
		std::cout << "dt missing in " + filePath << std::endl;
	}
	if (!config.Get("nTimesteps", nTimesteps)) {
		std::cout << "nTimesteps missing in " + filePath << std::endl;
	}
	if (!config.Get("title", title)) {
		std::cout << "title missing in " + filePath << std::endl;
	}
	if (!config.Get("precission", precission)) {
		std::cout << "precission missing in " + filePath << std::endl;
	}
}

SimulationData::SimulationData(int id, std::string title, int n_Stars, double boxLength, double dt, int nTimesteps) {
	this->simulationID = id;
	this->title = title;
	this->n_Stars = n_Stars;
	this->boxLength = boxLength;
	this->dt = dt;
	this->nTimesteps = nTimesteps;
}

SimulationData::SimulationData(int id){
	this->simulationID = id;
	SimulationData();
}

std::string SimulationData::print() {
	return "ID: " + std::to_string(this->simulationID) +
		" | Title: " + this->title +
		" #Stars: " + std::to_string(this->n_Stars) + '\n';
}

double SimulationData::getBoxLength() {
	return boxLength;
}

int SimulationData::getNStars() {
	return n_Stars;
}

double SimulationData::getdt() {
	return dt;
}

int SimulationData::getNTimesteps() {
	return nTimesteps;
}

double SimulationData::getG() {
	return G;
}

double SimulationData::getSoftening() {
	return softening;
}

double SimulationData::getPrecission() {
	return this->precission;
}