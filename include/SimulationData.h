#pragma once

#include <string>

#include "Parameters.h"

class SimulationData
{
protected:
    std::string title;
    int n_Stars;
    double boxLength;
    double dt;
    int nTimesteps;
    int simulationID;
public:
    SimulationData(int id, std::string title, int n_Stars, double boxLength, double dt, int nTimesteps);
    SimulationData(int id, Parameters* parameters);
    std::string print();
};