#include "Integrator.h"

double Integrator::dayInSec = 86400;
double Integrator::kmInpc = 3.086e-13;

Integrator::Integrator(double dt){
	this->dt = dt;
	this->dt2 = 0.5 * dt;
	this->C[0]=new Vec3D();
	this->C[1] = new Vec3D();
	this->C[2] = new Vec3D();
	this->C[3] = new Vec3D();
	this->K[0] = new Vec3D();
	this->K[1] = new Vec3D();
	this->K[2] = new Vec3D();
	this->K[3] = new Vec3D();
}

void Integrator::euler(std::vector<Star*> stars, double dt){
	if (dt != 0) {
		this->dt = dt;
	}
	//#pragma omp parallel for //1:10
	for (int i = 0; i < stars.size(); ++i){
		//root->applyForce(stars.at(i)->position, &stars.at(i)->acceleration);
		stars.at(i)->velocity += (stars.at(i)->acceleration * this->dt * this->dayInSec);
		stars.at(i)->position += (stars.at(i)->velocity * this->dt*this->kmInpc* this->dayInSec);
	}
}


//UNITS!!!
void Integrator::RK4(std::vector<Star*> stars, Node* root, double dt){
	if (dt != 0) {
		this->dt = dt;
		this->dt2 = 0.5 * dt;
	}
	for (int i = 0; i < stars.size(); ++i) {
		C[0] = &stars.at(i)->velocity; //carefull not to manipulate C!
		root->applyForce(stars.at(i)->position, K[0]);
		C[1] = &(stars.at(i)->velocity + (this->dt2 * *K[0]));
		root->applyForce(stars.at(i)->position + (this->dt2 * *C[0]) , K[1]);
		C[2] = &(stars.at(i)->velocity + (this->dt2 * *K[1]));
		root->applyForce(stars.at(i)->position + (this->dt2 * *C[1]), K[2]);
		C[3] = &(stars.at(i)->velocity + (this->dt * *K[1]));
		root->applyForce(stars.at(i)->position + (this->dt * *C[2]), K[3]);
		stars.at(i)->position += this->dt / 6. * (*C[0] + 2. * *C[1] + 2. * *C[2] + *C[3]);
		stars.at(i)->velocity += this->dt / 6. * (*K[0] + 2. * *K[1] + 2. * *K[2] + *K[3]);
	}
}
