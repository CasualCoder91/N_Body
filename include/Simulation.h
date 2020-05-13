/**
 * Simulation shared between different classes. 
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
    Potential* potential;
    
public:
    Simulation(int id, Database* database);
    Simulation(int id);
    Simulation(int id, SimulationData* simulationData);
    Simulation(int id, Database* database, Parameters* parameters);
    ~Simulation();

    void setID(int id);
    int getID();
    std::string print();
    void run();
private:
};

#endif