#pragma once

#include <string>

#include "Configuration.h"
#include "SimulationData.h"

class Parameters : public SimulationData{
protected:

    //Analysis Parameters
    bool bEnergy;
    bool bAverageVelocity;
    bool bAverage2DVelocity;
public:
    Parameters();

    std::string getTitle();
    std::string getFilePath();
    bool doEnergyAnalysis();
    bool doAverageVelocity();
    bool doAverage2DVelocity();
};

