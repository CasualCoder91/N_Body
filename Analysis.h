#pragma once

#include <vector>

#include "Star.h"
#include "Parameters.h"

class Analysis : Parameters
{
public:
	double static PotentialEnergy(std::vector<Star*>& stars);
	double static KineticEnergy(std::vector<Star*>& stars);
};

