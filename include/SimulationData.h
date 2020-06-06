/**
 * Stores all parameters needed for the simulation. 
 * Parameters are loaded from the simulation.cfg or from the database.
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#pragma once

#include <string>
#include <iostream>
#include "Configuration.h"

class SimulationData
{
protected:

    std::string filePath = "./simulation.cfg";
    Configuration config = Configuration();

    double softening = 0.16;
    double precission = 0.5;
    std::string title = "Default Title";
    double boxLength = 1;
    double dt = 1; //[day]
    int nTimesteps = 1000;
    int outputTimestep = 1;
    int simulationID = 1;

    //cluster mass
    double minMass = 0.1;
    double maxMass = 125;
    double alpha = -2.35;
    std::vector<double> massLimits; //broken powerlaw
    std::vector<double> exponents; //broken powerlaw


    double angle = 1; //degrees
    double dx = 10; //pc
    double distance = 8000; //pc
    Vec3D viewPoint = Vec3D(8000, 0, 0);
    Vec3D focus = Vec3D(0, 0, 0);


    //General Parameters, todo: put into separate class
    /** @brief Gravitational constant in astronomical units: parsec/solar mass*km^2/s^2*/
    static double G;// = 4.483e-3;
    int nStars = 0;
    Vec3D clusterLocation = Vec3D(40, 40, 0);

private:
    void getParametersFromConfig();
protected:
    SimulationData();
    void initParameterFromCfg(std::string name, double& value);
    void initParameterFromCfg(std::string name, int& value);
    void initParameterFromCfg(std::string name, std::string& value);
    void initParameterFromCfg(std::string name, Vec3D& value);
public:
    /**@brief: All values except \p id are read from simulation.cfg */
    SimulationData(int id);
    /**@brief: used when loading from database */
    SimulationData(int id, std::string title,int nStars, double boxLength, double dt, int nTimesteps, int outputTimestep);
    SimulationData(int id, SimulationData* simulationData);
    /**@brief: Simulation id is not set! */
    SimulationData(SimulationData* simulationData);
    /**@brief: returns a formated string containing ID, title and amount of cluster stars*/
    std::string print();
    double getBoxLength(); //[pc]
    double getdt();
    int getNTimesteps();
    double getSoftening();
    double getPrecission();
    int getOutputTimestep();
    double getMinMass();
    double getMaxMass();
    double getAlpha();
    double getG();
    int getNStars();
    /** @brief Sets amount of cluster stars. If possible use the simulation.cfg file instead*/
    void setNStars(int N);
};