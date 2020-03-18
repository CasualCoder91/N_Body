#pragma once

#include <string>
#include <iostream>
//#include "Parameters.h"
#include "Configuration.h"

class SimulationData
{
protected:
    std::string filePath;
    Configuration config;

    double softening;
    double precission;

    std::string title;
    int n_Stars;
    double boxLength;
    double dt;
    int nTimesteps;
    int simulationID;
    /** @brief Gravitational constant in astronomical units: parsec*solar mass*km^2/s^2*/
    double G;
public:
    SimulationData();
    SimulationData(int id, std::string title, int n_Stars, double boxLength, double dt, int nTimesteps);
    SimulationData(int id);
    std::string print();
    double getBoxLength(); //[pc]
    int getNStars();
    double getdt();
    int getNTimesteps();
    double getG();
    double getSoftening();
    double getPrecission();
};