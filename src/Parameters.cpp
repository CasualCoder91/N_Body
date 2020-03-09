#include "Parameters.h"

double Parameters::G = 4.483e-3;
double Parameters::softening = 0.16;
int Parameters::n_Stars;
double Parameters::boxLength; //[pc]
double Parameters::dt;
int Parameters::nTimesteps;
bool Parameters::energyAnalysis;
std::string Parameters::title;


Configuration   config;

Parameters::Parameters(){
    std::string file = "./simulation.cfg";
    config.Load(file);
    if (!config.Get("n_Stars", n_Stars)) {
        std::cout << "n_Stars missing in " + file<< std::endl;
    }
    if (!config.Get("boxLength", boxLength)) {
        std::cout << "boxLength missing in " + file << std::endl;
    }
    if (!config.Get("dt", dt)) {
        std::cout << "dt missing in " + file << std::endl;
    }
    if (!config.Get("nTimesteps", nTimesteps)) {
        std::cout << "nTimesteps missing in " + file << std::endl;
    }
    if (!config.Get("nTimesteps", nTimesteps)) {
        std::cout << "nTimesteps missing in " + file << std::endl;
    }
    if (!config.Get("title", title)) {
        std::cout << "title missing in " + file << std::endl;
    }
}

double Parameters::getG(){
    return G;
}

double Parameters::getSoftening(){
    return softening;
}

int Parameters::getN_Stars(){
    return n_Stars;
}

double Parameters::getBoxLength(){
    return boxLength;
}

double Parameters::getdt(){
    return dt;
}

int Parameters::getNTimesteps(){
    return nTimesteps;
}

bool Parameters::getEnergyAnalysis(){
    return energyAnalysis;
}

std::string Parameters::getTitle(){
    return this->title;
}
