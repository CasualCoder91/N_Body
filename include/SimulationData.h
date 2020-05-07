#pragma once

#include <string>
#include <iostream>
#include "Configuration.h"

class SimulationData
{
protected:
    static std::string filePath;
    static Configuration config;

    double softening;
    double precission;
    std::string title;
    static int nStars;
    double boxLength;
    double dt; //[day]
    int nTimesteps;
    int outputTimestep;
    int simulationID;
    /** @brief Gravitational constant in astronomical units: parsec/solar mass*km^2/s^2*/
    static double G;
    double minMass;
    double maxMass;
    double alpha;
private:
    void getParametersFromConfig();
public:
    SimulationData();
    //used when loading from Database
    SimulationData(int id, std::string title, int nStars, double boxLength, double dt, int nTimesteps, int outputTimestep);
    SimulationData(int id);
    SimulationData(int id, SimulationData* simulationData);
    std::string print();
    double getBoxLength(); //[pc]
    int getNStars();
    static int getNStarsFromCfg();
    double getdt();
    int getNTimesteps();
    double getG();
    static double getGFromCfg();
    double getSoftening();
    double getPrecission();
    int getOutputTimestep();
    double getMinMass();
    double getMaxMass();
    double getAlpha();
};