#pragma once

#include <string>
#include <iostream>
#include "Configuration.h"

class SimulationData
{
protected:
    std::string filePath;
    Configuration config;

    double softening;
    double precission;
    std::string title;
    int nStars;
    double boxLength;
    double dt;
    int nTimesteps;
    int outputTimestep;
    int simulationID;
    /** @brief Gravitational constant in astronomical units: parsec*solar mass*km^2/s^2*/
    double G;
    double minMass;
    double maxMass;
    double alpha;
public:
    SimulationData();
    SimulationData(int id, std::string title, int nStars, double boxLength, double dt, int nTimesteps);
    SimulationData(int id);
    SimulationData(int id, SimulationData* simulationData);
    std::string print();
    double getBoxLength(); //[pc]
    int getNStars();
    double getdt();
    int getNTimesteps();
    double getG();
    double getSoftening();
    double getPrecission();
    int getOutputTimestep();
    double getMinMass();
    double getMaxMass();
    double getAlpha();
};