#pragma once

#include <string>

#include "Configuration.h"

class Parameters{
protected:
    std::string title;
    /** @brief Gravitational constant in astronomical units: parsec*solar mass*km^2/s^2*/
    double G;
    double softening;
    double precission;
    int n_Stars;
    double boxLength; //[pc]
    double dt;
    int nTimesteps;
    std::string filePath;
    Configuration config;

    //Analysis Parameters
    bool bEnergy;
    bool bAverageVelocity;
    bool bAverage2DVelocity;
public:
    Parameters();

    double getG();
    double getSoftening();
    double getPrecission();
    double getBoxLength(); //[pc]
    double getdt();
    int getNTimesteps();
    int getNStars();
    std::string getTitle();
    std::string getFilePath();
    bool doEnergyAnalysis();
    bool doAverageVelocity();
    bool doAverage2DVelocity();
};

