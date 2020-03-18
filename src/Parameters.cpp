#include "Parameters.h"


Parameters::Parameters(){
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

double Parameters::getG() {
    return G;
}
double Parameters::getSoftening() {
    return softening;
}

double Parameters::getPrecission(){
    return this->precission;
}

double Parameters::getBoxLength() {
    return boxLength;
}
double Parameters::getdt() {
    return dt;
}
int Parameters::getNTimesteps() {
    return nTimesteps;
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

int Parameters::getNStars() {
    return n_Stars;
}
