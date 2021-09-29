#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "Vec3D.h"
#include "Star.h"
#include "Projection.h"
#include "Constants.h"

class Extinction {
private:
	std::vector<std::vector<double>> map;

public:
	Extinction();
	void set_extinction(Star& star);
};