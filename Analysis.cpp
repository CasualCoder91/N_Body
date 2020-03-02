#include "Analysis.h"

double Analysis::PotentialEnergy(std::vector<Star*>& stars){
	double potentialEnergy = 0;
	for (unsigned int i = 0; i < stars.size()-1;++i) {
		for (int j = i+1; j < stars.size(); ++j) {
			potentialEnergy -= Parameters::G * stars.at(i)->mass * stars.at(j)->mass / Vec3D::distance(&stars.at(i)->position, &stars.at(j)->position);
		}
	}
	return potentialEnergy;
}

double Analysis::KineticEnergy(std::vector<Star*>& stars){
	double kineticEnergy = 0;
	for (unsigned int i = 0; i < stars.size(); ++i) {
		kineticEnergy += stars.at(i)->mass * (stars.at(i)->velocity * stars.at(i)->velocity);
	}
	kineticEnergy *= 0.5;
	return kineticEnergy;
}
