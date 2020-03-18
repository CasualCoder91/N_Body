#include "Parameters.h"


Parameters::Parameters():SimulationData(){

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
    if (!config.Get("softening", softening)) {
        this->softening = 0.16;
        std::cout << "softening missing in " + filePath + " setting default value: " + std::to_string(softening) << std::endl;
    }
    if (!config.Get("G", G)) {
        this->G = 4.483e-3;
        std::cout << "G missing in " + filePath + " setting default value: " + std::to_string(G) << std::endl;
    }
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


