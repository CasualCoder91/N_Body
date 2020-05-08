#include "Parameters.h"

void Parameters::getParametersFromConfig(){
    //std::cout << "\nReading Parameters from " << filePath << std::endl;
    //config.Load(filePath);
}

Parameters::Parameters():SimulationData(){
	config.Load(filePath);
    if (!config.Get("bEnergy", bEnergy)) {
        bEnergy = false;
        std::cout << "bEnergy missing in " + filePath + " add bEnergy = true to activate energy analysis." << std::endl;
    }
    if (!config.Get("bAverageVelocity", bAverageVelocity)) {
        bAverageVelocity = false;
        std::cout << "bAverageVelocity missing in " + filePath + " add bAverageVelocity = true to activate 3d velocity analysis." << std::endl;
    }
    if (!config.Get("bAverage2DVelocity", bAverage2DVelocity)) {
        bAverage2DVelocity = false;
        std::cout << "bAverage2DVelocity missing in " + filePath + " add bAverage2DVelocity = true to activate 2d velocity analysis." << std::endl;
    }
    getParametersFromConfig();
}

Parameters::Parameters(Parameters* parameters):SimulationData((SimulationData*)parameters){
	this->bEnergy = parameters->bEnergy;
	this->bAverageVelocity = parameters->bAverageVelocity;
	this->bAverage2DVelocity = parameters->bAverage2DVelocity;
}

std::string Parameters::getTitle() {
    return title;
}

std::string Parameters::getFilePath(){
    return this->filePath;
}

bool Parameters::doEnergyAnalysis(){
    return this->bEnergy;
}

bool Parameters::doAverageVelocity(){
    return this->bAverageVelocity;
}

bool Parameters::doAverage2DVelocity(){
    return this->bAverage2DVelocity;
}