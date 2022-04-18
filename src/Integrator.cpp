#include "Integrator.h"

double Integrator::secInDay = 86400;
//double Integrator::kmInpc = 3.086e-13;

Integrator::Integrator(double dt){
	this->dt = dt;
	this->dt2 = 0.5 * dt;
	this->dts = dt * secInDay;
	this->dts2 = 0.5 * dts;
	this->positionDetla = dts * Constants::kmInpc;
	this->positionDetla2 = 0.5 * positionDetla;
}

void Integrator::euler(std::vector<Star>& stars, double dt){
	if (dt != 0) {
		this->dt = dt;
	}
	//#pragma omp parallel for //1:10
	for (int i = 0; i < stars.size(); ++i){
		//root->applyForce(stars[i]->position, &stars[i]->acceleration);
		stars[i].velocity += (stars[i].acceleration * this->dts);
		stars[i].position += (stars[i].velocity * this->positionDetla);
	}
}


void Integrator::RK4(std::vector<Star*> stars, Node* root, MWPotential* potential){
	for (size_t i = 0; i < stars.size(); ++i) {
		Vec3D C[4];
		Vec3D K[4];
		C[0] = stars[i]->velocity; //carefull not to manipulate C!
		root->applyForce(stars[i]->position, K[0]);
		potential->applyForce(stars[i]->position, K[0]);
		C[1] = (stars[i]->velocity + (this->dts2 * K[0]));
		root->applyForce(stars[i]->position + (positionDetla2 * C[0]) , K[1]);
		potential->applyForce(stars[i]->position + (positionDetla2 * C[0]), K[1]);
		C[2] = (stars[i]->velocity + (this->dts2 * K[1]));
		root->applyForce(stars[i]->position + (positionDetla2 * C[1]), K[2]);
		potential->applyForce(stars[i]->position + (positionDetla2 * C[1]), K[2]);
		C[3] = (stars[i]->velocity + (this->dts * K[2]));
		root->applyForce(stars[i]->position + (positionDetla * C[2]), K[3]);
		potential->applyForce(stars[i]->position + (positionDetla * C[2]), K[3]);
		stars[i]->position += this->positionDetla / 6. * (C[0] + 2. * C[1] + 2. * C[2] + C[3]);
		stars[i]->velocity += this->dts / 6. * (K[0] + 2. * K[1] + 2. * K[2] + K[3]);
	}
}

void Integrator::RK4(std::vector<Star*> stars, MWPotential* potential){
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
		stars[i]->position += this->positionDetla / 6. * (C[0] + 2. * C[1] + 2. * C[2] + C[3]);
		stars[i]->velocity += this->dts / 6. * (K[0] + 2. * K[1] + 2. * K[2] + K[3]);
	}
}

void Integrator::Leapfrog(std::vector<Star>& stars, Node* root, MWPotential* potential){
	double spaceFactor = Constants::kmInpc * this->secInDay;
	for (size_t i = 0; i < stars.size(); ++i) {
		Vec3D xHalfStep = stars[i].position + stars[i].velocity * this->positionDetla2;
		Vec3D aHalfStep;
		root->applyForce(xHalfStep, aHalfStep);
		potential->applyForce(xHalfStep, aHalfStep);
		stars[i].velocity = stars[i].velocity + aHalfStep * this->dts;
		stars[i].position = xHalfStep + stars[i].velocity * this->positionDetla2;
	}
}

void Integrator::Leapfrog(std::vector<Star>& stars, MWPotential* potential){
	double spaceFactor = Constants::kmInpc * this->secInDay;
	for (size_t i = 0; i < stars.size(); ++i) {
		Vec3D xHalfStep = stars[i].position + stars[i].velocity * this->dt2 * spaceFactor;
		Vec3D aHalfStep;
		potential->applyForce(xHalfStep, aHalfStep);
		stars[i].velocity = stars[i].velocity + aHalfStep * this->dt * this->secInDay;
		stars[i].position = xHalfStep + stars[i].velocity * this->dt2 * spaceFactor;
	}
}
