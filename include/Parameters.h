/**
 * Contains all parameters read from the "simulation.cfg" file. If parameters are not specified in the file they are set to default values.
 * Parameters which are stored into the simulation database table are read into the parent class SimulationData.
 *
 * @author Alarich Herzner
 * @version 0.2 15.05.2020
*/

#pragma once

#include <string>

#include "Configuration.h"
#include "SimulationData.h"

class Parameters : public SimulationData{
protected:
    std::string filePath = "./simulation.cfg";
    Configuration config = Configuration(filePath);

    //Analysis Parameters
    bool bEnergy = false;
    bool bAverageVelocity = false;
    bool bAverage2DVelocity = false;

protected:
    void getParametersFromConfig();
public:
    /**@brief default constructor. Parameter are read from the file*/
    Parameters();
    /**@brief Copy constructor*/
    Parameters(Parameters* parameters);

    std::string getTitle();
    std::string getFilePath();

    bool doEnergyAnalysis();
    bool doAverageVelocity();
    bool doAverage2DVelocity();
};

