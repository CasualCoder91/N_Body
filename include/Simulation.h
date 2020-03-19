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

class Simulation : SimulationData
{
private:
    Database* database;
    //Parameters* parameters;
    
public:
    Simulation(int id, Database* database);
    Simulation(int id);
    Simulation(int id, SimulationData* simulationData);
    void setID(int id);
    int getID();
    /*int getNStars();
    int getNTimesteps();
    double getdt();*/
    std::string print();
    void run();
};

#endif