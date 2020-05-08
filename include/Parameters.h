#pragma once

#include <string>

#include "Configuration.h"
#include "SimulationData.h"

class Parameters : public SimulationData{
protected:
    std::string filePath = "./simulation.cfg";
    Configuration config = Configuration();

    //Analysis Parameters
    bool bEnergy = false;
    bool bAverageVelocity = false;
    bool bAverage2DVelocity = false;

protected:
    void getParametersFromConfig();
public:
    Parameters();
    Parameters(Parameters* parameters);

    std::string getTitle();
    std::string getFilePath();

    bool doEnergyAnalysis();
    bool doAverageVelocity();
    bool doAverage2DVelocity();
};

