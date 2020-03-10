#include "Simulation.h"

Simulation::Simulation(int id){
	this->simulationID = id;
}

Simulation::Simulation(int id, std::string title, int n_Stars, double boxLength, double dt, int nTimesteps){
	this->simulationID = id;
	this->title = title;
	this->n_Stars = n_Stars;
	this->boxLength = boxLength;
	this->dt = dt;
	this->nTimesteps = nTimesteps;
}

void Simulation::setID(int id){
	this->simulationID = id;
}

int Simulation::getID(){
	return this->simulationID;
}

std::string Simulation::print(){
	return "ID: " + std::to_string(this->simulationID) +
		" | Title: " + this->title + 
		" #Stars: " + std::to_string(this->n_Stars) +'\n';
}
