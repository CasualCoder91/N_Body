/**
 * Parameters shared between different classes. 
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#include "Configuration.h"

#pragma once
class Parameters
{
public:
    Parameters();
protected:
    /** @brief Gravitational constant in astronomical units: parsec*solar mass*km^2/s^2*/
    static std::string title;
	static double G;
    static double softening;
    static int n_Stars;
    static double boxLength; //[pc]
    static double dt;
    static int nTimesteps;
    static bool energyAnalysis;
public:
    double getG();
    double getSoftening();
    int getN_Stars();
    double getBoxLength(); //[pc]
    double getdt();
    int getNTimesteps();
    bool getEnergyAnalysis();
    std::string getTitle();
};

