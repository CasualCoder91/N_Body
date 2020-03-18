#include "SimulationData.h"


SimulationData::SimulationData(int id, std::string title, int n_Stars, double boxLength, double dt, int nTimesteps) {
	this->simulationID = id;
	this->title = title;
	this->n_Stars = n_Stars;
	this->boxLength = boxLength;
	this->dt = dt;
	this->nTimesteps = nTimesteps;
}

SimulationData::SimulationData(int id, Parameters* parameters){
	this->simulationID = id;
	this->title = parameters->getTitle();
	this->n_Stars = parameters->getNStars();
	this->boxLength = parameters->getBoxLength();
	this->dt = parameters->getdt();
	this->nTimesteps = parameters->getNTimesteps();
}

std::string SimulationData::print() {
	return "ID: " + std::to_string(this->simulationID) +
		" | Title: " + this->title +
		" #Stars: " + std::to_string(this->n_Stars) + '\n';
}