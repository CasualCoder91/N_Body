#include "Star.h"

Star::Star(){
	this->mass = 1;
	this->position = Vec3D();
	this->velocity = Vec3D();
	this->acceleration = Vec3D();
}

Star::Star(double mass)
{
	this->mass = mass;
	this->position = Vec3D();
	this->velocity = Vec3D();
	this->acceleration = Vec3D();
}

Star::Star(double mass, Vec3D position){
	this->mass = mass;
	this->position = position;
	this->velocity = Vec3D();
	this->acceleration = Vec3D();
}

Star::Star(double mass, double x, double y, double z){
	this->mass = mass;
	this->position = Vec3D(x, y, z);
	this->velocity = Vec3D();
	this->acceleration = Vec3D();
}

std::string Star::Dump()
{
	return "mass: " + std::to_string(mass) + '\n'
		+ "position: " + position.Print() + '\n'
		+ "velocity: " + velocity.Print() + '\n'
		+ "acceleration: " + acceleration.Print()+'\n';
}

void Star::Reset(){
	this->mass = 0; 
	if(this->position.x)
		this->position.Reset(); 
	if(this->velocity.x)
		this->velocity.Reset(); 
	if(this->acceleration.x)
		this->acceleration.Reset();
}
