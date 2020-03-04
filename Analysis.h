#pragma once

#include <vector>
#include <chrono> //for timer

#include "Star.h"
#include "Parameters.h"
#include "Integrator.h"
#include "Node.h"
#include "InOut.h"

class Analysis : Parameters
{
public:
	double static PotentialEnergy(std::vector<Star*>& stars);
	double static KineticEnergy(std::vector<Star*>& stars);
	void static scaling(int maxN, int nTimesteps, std::vector<Star*>& stars, Integrator& integrator);
};

