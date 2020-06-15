/**
 * The heart of the program. Used to setup and run a specific simulation.
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>

#include "Database.h"
#include "Parameters.h"
#include "ProgressBar.h"

class Simulation : public Parameters{

public:
    Database* database;
    MWPotential* potential;
    
public:
    Simulation(int id, Database* database, SimulationData* simulationData);
    Simulation(int id, Database* database, Parameters* parameters);
    ~Simulation();

    void setID(int id);
    int getID();
    /**
    Cluster and field stars are initialized and stored in the database. In each integration step 1. the acceleration
    of the stars are updated based on the potential and in case of the cluster stars based on each other and 2. the timestep is 
    integrated by one of the available integration methods. At defined timestep intervals the current values for position and velocity are stored into the database.
    */
    void run();
private:
    Simulation(int id); //todo: do not need this?
};

#endif