/**
 * Simulation shared between different classes. 
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "Parameters.h"

class Simulation : public Parameters
{
private:
protected:
    int simulationID;
public:
    Simulation(int id);
    Simulation(int id, std::string title, int n_Stars, double boxLength, double dt, int nTimesteps);
    void setID(int id);
    int getID();
    std::string print();
};

#endif