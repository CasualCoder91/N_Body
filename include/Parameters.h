#pragma once

#include <string>

#include "Configuration.h"

class Parameters{
protected:
    std::string title;
    /** @brief Gravitational constant in astronomical units: parsec*solar mass*km^2/s^2*/
    double G;
    double softening;
    int n_Stars;
    double boxLength; //[pc]
    double dt;
    int nTimesteps;
    std::string filePath;
    Configuration   config;
public:
    Parameters();

    double getG();
    double getSoftening();
    int getN_Stars();
    double getBoxLength(); //[pc]
    double getdt();
    int getNTimesteps();
    std::string getTitle();
};

