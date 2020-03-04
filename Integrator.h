#pragma once
#include "Star.h"
#include <vector>

class Integrator
{
public:
	double dt = 0.1; //timestepsize

	Integrator(double dt = 0.01);

	void euler(std::vector<Star*> stars,double dt= 0); // one timestep with euler algorithm
};

