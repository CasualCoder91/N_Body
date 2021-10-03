#include "Star.h"

Star::Star(int id){
	this->id = id;
	this->mass = 1;
	this->position = Vec3D();
	this->velocity = Vec3D();
	this->acceleration = Vec3D();
	this->extinction = 0;
}

Star::Star(int id, double mass){
	this->id = id;
	this->mass = mass;
	this->position = Vec3D();
	this->velocity = Vec3D();
	this->acceleration = Vec3D();
	this->extinction = 0;
}

Star::Star(int id, double mass, Vec3D position){
	this->id = id;
	this->mass = mass;
	this->position = position;
	this->velocity = Vec3D();
	this->acceleration = Vec3D();
	this->extinction = 0;
}

Star::Star(int id, double mass, double x, double y, double z){
	this->id = id;
	this->mass = mass;
	this->position = Vec3D(x, y, z);
	this->velocity = Vec3D();
	this->acceleration = Vec3D();
	this->extinction = 0;
}

Star::Star(int id, double mass, double xPos, double yPos, double zPos, double xVel, double yVel, double zVel){
	this->id = id;
	this->mass = mass;
	this->position = Vec3D(xPos, yPos, zPos);
	this->velocity = Vec3D(xVel, yVel, zVel);
	this->acceleration = Vec3D();
	this->extinction = 0;
}

std::string const Star::dump()
{
	return "mass: " + std::to_string(mass) + '\n'
		+ "position: " + position.print() + '\n'
		+ "velocity: " + velocity.print() + '\n'
		+ "acceleration: " + acceleration.print()+'\n';
}

void Star::reset(){
	this->mass = 0; 
	if(this->position.x)
		this->position.reset(); 
	if(this->velocity.x)
		this->velocity.reset(); 
	if(this->acceleration.x)
		this->acceleration.reset();
	this->extinction = 0;
}
