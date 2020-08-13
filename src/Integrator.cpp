#include "Integrator.h"

double Integrator::secInDay = 86400;
double Integrator::kmInpc = 3.086e-13;

Integrator::Integrator(double dt){
	this->dt = dt;
	this->dt2 = 0.5 * dt;
}

void Integrator::euler(std::vector<Star*> stars, double dt){
	if (dt != 0) {
		this->dt = dt;
	}
	//#pragma omp parallel for //1:10
	for (int i = 0; i < stars.size(); ++i){
		//root->applyForce(stars[i]->position, &stars[i]->acceleration);
		stars[i]->velocity += (stars[i]->acceleration * this->dt * this->secInDay);
		stars[i]->position += (stars[i]->velocity * this->dt*this->kmInpc* this->secInDay);
	}
}


//UNITS!!!
void Integrator::RK4(std::vector<Star*> stars, Node* root, MWPotential* potential, double dt){
	if (dt != 0) {
		this->dt = dt;
		this->dt2 = 0.5 * dt;
	}
	double spaceFactor = this->kmInpc * this->secInDay;
	for (size_t i = 0; i < stars.size(); ++i) {
		Vec3D C[4]; //RK4
		Vec3D K[4]; //RK4
		C[0] = stars[i]->velocity; //carefull not to manipulate C!
		root->applyForce(stars[i]->position, K[0]);
		potential->applyForce(stars[i]->position, K[0]);
		C[1] = (stars[i]->velocity + (this->dt2 * this->secInDay * K[0]));
		root->applyForce(stars[i]->position + (this->dt2 * spaceFactor * C[0]) , K[1]);
		potential->applyForce(stars[i]->position + (this->dt2 * spaceFactor * C[0]), K[1]);
		C[2] = (stars[i]->velocity + (this->dt2 * this->secInDay * K[1]));
		root->applyForce(stars[i]->position + (this->dt2 * spaceFactor * C[1]), K[2]);
		potential->applyForce(stars[i]->position + (this->dt2 * spaceFactor * C[1]), K[2]);
		C[3] = (stars[i]->velocity + (this->dt * this->secInDay * K[2]));
		root->applyForce(stars[i]->position + (this->dt * spaceFactor * C[2]), K[3]);
		potential->applyForce(stars[i]->position + (this->dt * spaceFactor * C[2]), K[3]);
		stars[i]->position += this->dt / 6. * (C[0] + 2. * C[1] + 2. * C[2] + C[3]) * spaceFactor;
		stars[i]->velocity += this->dt / 6. * (K[0] + 2. * K[1] + 2. * K[2] + K[3]) * this->secInDay;
	}
}

void Integrator::RK4(std::vector<Star*> stars, MWPotential* potential, double dt){
	if (dt != 0) {
		this->dt = dt;
		this->dt2 = 0.5 * dt;
	}
	Vec3D C[4]; //RK4
	Vec3D K[4]; //RK4
	for (size_t i = 0; i < stars.size(); ++i) {
		C[0] = stars[i]->velocity; //carefull not to manipulate C!
		potential->applyForce(stars[i]->position, K[0]);
		C[1] = (stars[i]->velocity + (this->dt2 * K[0]));
		potential->applyForce(stars[i]->position + (this->dt2 * C[0]), K[1]);
		C[2] = (stars[i]->velocity + (this->dt2 * K[1]));
		potential->applyForce(stars[i]->position + (this->dt2 * C[1]), K[2]);
		C[3] = (stars[i]->velocity + (this->dt * K[1]));
		potential->applyForce(stars[i]->position + (this->dt * C[2]), K[3]);
		stars[i]->position += this->dt / 6. * (C[0] + 2. * C[1] + 2. * C[2] + C[3]) * this->kmInpc * this->secInDay;
		stars[i]->velocity += this->dt / 6. * (K[0] + 2. * K[1] + 2. * K[2] + K[3]) * this->secInDay;
	}
}
