#include "Integrator.h"

Integrator::Integrator(double dt){
	this->dt = dt;
}

void Integrator::Euler(std::vector<Star*> stars,double dt){
	this->dt = dt;
	#pragma omp parallel for //1:10
	for (int i = 0; i < stars.size(); ++i){
		//star->velocity += star->acceleration * dt;
		//star->position +=  star->velocity * dt;
		stars.at(i)->velocity.x += stars.at(i)->acceleration.x * dt;
		stars.at(i)->velocity.y += stars.at(i)->acceleration.y * dt;
		stars.at(i)->velocity.z += stars.at(i)->acceleration.z * dt;
		stars.at(i)->position.x += stars.at(i)->velocity.x * dt;
		stars.at(i)->position.y += stars.at(i)->velocity.y * dt;
		stars.at(i)->position.z += stars.at(i)->velocity.z * dt;
	}
}
