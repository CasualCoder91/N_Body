#include "SimulationData.h"


SimulationData::SimulationData(){
	this->filePath = "./simulation.cfg";
	config.Load(filePath);
	if (!config.Get("nStars", nStars)) {
		std::cout << "nStars missing in " + filePath << std::endl;
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
	if (!config.Get("outputTimestep", outputTimestep)) {
		std::cout << "outputTimestep missing in " + filePath << std::endl;
	}
	if (!config.Get("title", title)) {
		std::cout << "title missing in " + filePath << std::endl;
	}
	if (!config.Get("softening", softening)) {
		std::cout << "softening missing in " + filePath << std::endl;
	}
	if (!config.Get("precission", precission)) {
		std::cout << "precission missing in " + filePath << std::endl;
	}
	if (!config.Get("minMass", minMass)) {
		std::cout << "minMass missing in " + filePath << std::endl;
	}
	if (!config.Get("maxMass", maxMass)) {
		std::cout << "maxMass missing in " + filePath << std::endl;
	}
	if (!config.Get("alpha", alpha)) {
		std::cout << "precission missing in " + filePath << std::endl;
	}
}

SimulationData::SimulationData(int id, std::string title, int nStars, double boxLength, double dt, int nTimesteps) {
	this->simulationID = id;
	this->title = title;
	this->nStars = nStars;
	this->boxLength = boxLength;
	this->dt = dt;
	this->nTimesteps = nTimesteps;
}

SimulationData::SimulationData(int id){
	this->simulationID = id;
	SimulationData();
}

SimulationData::SimulationData(int id, SimulationData* simulationData){
	this->simulationID = id;
	this->boxLength = simulationData->boxLength;
	this->dt = simulationData->dt;
	this->G = simulationData->G;
	this->alpha = simulationData->alpha;
	this->maxMass = simulationData->maxMass;
	this->minMass = simulationData->minMass;
	this->nTimesteps = simulationData->nTimesteps;
	this->nStars = simulationData->nStars;
	this->outputTimestep = simulationData->outputTimestep;
	this->softening = simulationData->softening;
	this->precission = simulationData->precission;
}

std::string SimulationData::print() {
	return "ID: " + std::to_string(this->simulationID) +
		" | Title: " + this->title +
		" #Stars: " + std::to_string(this->nStars) + '\n';
}

double SimulationData::getBoxLength() {
	return boxLength;
}

int SimulationData::getNStars() {
	return nStars;
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

int SimulationData::getOutputTimestep(){
	return this->outputTimestep;
}

double SimulationData::getMinMass(){
	return this->minMass;
}

double SimulationData::getMaxMass(){
	return this->maxMass;
}

double SimulationData::getAlpha(){
	return this->alpha;
}