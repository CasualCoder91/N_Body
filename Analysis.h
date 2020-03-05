#pragma once

#include <vector>
#include <chrono> //for timer

#include "Star.h"
#include "Parameters.h"
#include "Integrator.h"
#include "Node.h"
#include "InOut.h"
#include "InitialConditions.h"

class Analysis : Parameters
{
public:
	double static potentialEnergy(std::vector<Star*>& stars);
	double static kineticEnergy(std::vector<Star*>& stars);
	void static scaling(int maxNStars, int nTimesteps, Integrator& integrator);
};

