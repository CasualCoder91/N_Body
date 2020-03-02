#include "Vec3D.h"

Vec3D::Vec3D(){
	x = y = z = 0;
}

Vec3D::Vec3D(double x, double y, double z){
	this->x = x;
	this->y = y;
	this->z = z;
}

Vec3D Vec3D::Center(Vec3D a, Vec3D b)
{
	return Vec3D((a.x+b.x)/2, (a.y + b.y)/2, (a.y + b.y)/2);
}

Vec3D Vec3D::RandomVector(double length)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(-1.0, 1.0);
	Vec3D vec = Vec3D(dis(gen), dis(gen), dis(gen));
	double currentLength = vec.Length();
	double factor = length / currentLength;
	vec.x *= factor;
	vec.y *= factor;
	vec.z *= factor;
	return vec;
}

Vec3D Vec3D::RandomAngles(double r){
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(-1.0, 1.0);
	std::uniform_real_distribution<> dis2(0.0, 6.28318530718); //0-2*pi
	double theta = acos(dis(gen));
	double phi = dis2(gen);
	Vec3D vec = Vec3D(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
	return vec;
}

double Vec3D::Length()
{
	return sqrt(this->x*this->x+this->y*this->y+this->z*this->z);
}

Vec3D Vec3D::Normalize(){
	double n = this->Length();
	return Vec3D(x/n,y/n,z/n);
}

Vec3D Vec3D::Cross(Vec3D v1, Vec3D v2){
	return Vec3D(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

std::string Vec3D::Print()
{
	return std::to_string(this->x) +','+ std::to_string(this->y) + ',' + std::to_string(this->z);
}

double Vec3D::Distance(Vec3D* a, Vec3D* b){
	return sqrt(pow(a->x-b->x,2)+ pow(a->y - b->y, 2)+ pow(a->z - b->z, 2));
}

void Vec3D::Reset(){
	this->x = 0; this->y = 0; this->z = 0;
}

Vec3D Vec3D::operator/(double& rhs){
	return Vec3D(this->x/rhs,this->y/rhs,this->z/rhs);
}

Vec3D Vec3D::operator*(double& rhs){
	return Vec3D(this->x * rhs, this->y * rhs, this->z * rhs);
}

double Vec3D::operator*(Vec3D& rhs){
	return this->x * rhs.x + this->y* rhs.y + this->z* rhs.z;
}

Vec3D Vec3D::operator+(Vec3D const& rhs){
	return Vec3D(this->x+rhs.x, this->y + rhs.y, this->z + rhs.z);
}

Vec3D& Vec3D::operator+=(const Vec3D & rhs){
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
	return *this;
}

Vec3D operator*(float f, const Vec3D& v){
	return Vec3D(v.x * f, v.y * f, v.z * f);
}
