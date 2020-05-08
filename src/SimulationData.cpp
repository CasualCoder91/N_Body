#include "SimulationData.h"

double SimulationData::G = 4.483e-3;

void SimulationData::getParametersFromConfig(){
	std::cout << "\nReading Parameters from " << filePath << std::endl;
	config.Load(filePath);
	initParameterFromCfg("boxLength", boxLength);
	initParameterFromCfg("nTimesteps", nTimesteps);
	initParameterFromCfg("outputTimestep", outputTimestep);
	initParameterFromCfg("title", title);
	initParameterFromCfg("precission", precission);
	initParameterFromCfg("minMass", minMass);
	initParameterFromCfg("maxMass", maxMass);
	initParameterFromCfg("softening", softening);
	initParameterFromCfg("alpha", alpha);
	initParameterFromCfg("angle", angle);
	initParameterFromCfg("dx", dx);
	initParameterFromCfg("distance", distance);
	initParameterFromCfg("viewPoint", viewPoint);
	initParameterFromCfg("focus", focus);
	initParameterFromCfg("nStars", nStars);
	initParameterFromCfg("G", G);
}

SimulationData::SimulationData(){
	getParametersFromConfig();
}

SimulationData::SimulationData(int id):SimulationData(){
	this->simulationID = id;
}

SimulationData::SimulationData(int id, std::string title,int nStars, double boxLength, double dt, int nTimesteps, int outputTimestep) {
	this->simulationID = id;
	this->title = title;
	this->nStars = nStars;
	this->boxLength = boxLength;
	this->dt = dt;
	this->nTimesteps = nTimesteps;
	this->outputTimestep = outputTimestep;
}

SimulationData::SimulationData(int id, SimulationData* simulationData):SimulationData(simulationData) {
	this->simulationID = id;
}

SimulationData::SimulationData(SimulationData* simulationData){
	this->boxLength = simulationData->boxLength;
	this->dt = simulationData->dt;
	this->alpha = simulationData->alpha;
	this->maxMass = simulationData->maxMass;
	this->minMass = simulationData->minMass;
	this->nTimesteps = simulationData->nTimesteps;
	this->outputTimestep = simulationData->outputTimestep;
	this->softening = simulationData->softening;
	this->precission = simulationData->precission;
	this->title = simulationData->title;
	this->focus = simulationData->focus;
	this->angle = simulationData->angle;
	this->dx = simulationData->dx;
	this->distance = simulationData->distance;
	this->viewPoint = simulationData->viewPoint;
	this->G = simulationData->G;
	this->nStars = simulationData->nStars;
}

std::string SimulationData::print() {
	return "ID: " + std::to_string(this->simulationID) +
		" | Title: " + this->title
		//+" #Stars: " + std::to_string(this->nStars) 
		+ '\n';
}

double SimulationData::getBoxLength() {
	return boxLength;
}

double SimulationData::getdt() {
	return dt;
}

int SimulationData::getNTimesteps() {
	return nTimesteps;
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

double SimulationData::getG() {
	return this->G;
}

int SimulationData::getNStars() {
	return nStars;
}


void SimulationData::initParameterFromCfg(std::string name, double& value){
	if (!config.Get(name, value)) {
		std::cout << name << " missing in " + filePath << std::endl;
		std::cout << "Using default value: " << value << std::endl;
	}
}

void SimulationData::initParameterFromCfg(std::string name, int& value){
	if (!config.Get(name, value)) {
		std::cout << name << " missing in " + filePath << std::endl;
		std::cout << "Using default value: " << value << std::endl;
	}
}

void SimulationData::initParameterFromCfg(std::string name, std::string& value){
	if (!config.Get(name, value)) {
		std::cout << name << " missing in " + filePath << std::endl;
		std::cout << "Using default value: " << value << std::endl;
	}
}

void SimulationData::initParameterFromCfg(std::string name, Vec3D& value){
	if (!config.Get(name, value)) {
		std::cout << name << " missing in " + filePath << std::endl;
		std::cout << "Using default value: " << value.print() << std::endl;
	}
}
