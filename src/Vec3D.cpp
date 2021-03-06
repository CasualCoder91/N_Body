#include "Vec3D.h"

Vec3D::Vec3D(){
	x = y = z = 0;
}

Vec3D::Vec3D(double x, double y, double z){
	this->x = x;
	this->y = y;
	this->z = z;
}

Vec3D::~Vec3D(){}

Vec3D Vec3D::center(Vec3D a, Vec3D b)
{
	return Vec3D((a.x+b.x)/2, (a.y + b.y)/2, (a.y + b.y)/2);
}

Vec3D Vec3D::randomVector(double length)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(-1.0, 1.0);
	Vec3D vec = Vec3D(dis(gen), dis(gen), dis(gen));
	double currentLength = vec.length();
	double factor = length / currentLength;
	vec.x *= factor;
	vec.y *= factor;
	vec.z *= factor;
	return vec;
}

Vec3D Vec3D::randomAngles(double r){
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(-1.0, 1.0);
	std::uniform_real_distribution<> dis2(0.0, 6.28318530718); //0-2*pi
	double theta = acos(dis(gen));
	double phi = dis2(gen);
	Vec3D vec = Vec3D(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
	return vec;
}

double Vec3D::length() const {
	return sqrt(this->x*this->x+this->y*this->y+this->z*this->z);
}

double Vec3D::phi(){
	return atan2(y, x);
}

double Vec3D::theta(){
	double r = this->length();
	if (r == 0)
		return 0;
	return acos(z / this->length());
}

Vec3D Vec3D::cartesianToCylindrical(){
	return Vec3D(sqrt(this->x* this->x+ this->y* this->y),this->phi(),this->z);
}

Vec3D Vec3D::cartesianToCylindricalV(Vec3D cartesianPos){
	double phi = atan2(cartesianPos.y,cartesianPos.x);
	return Vec3D(this->x * cos(phi) + this->y * sin(phi), -this->x * sin(phi) + this->y * cos(phi), this->z);
}

Vec3D Vec3D::cartesianToSphericalV(Vec3D cartesianPos){
	double r = cartesianPos.length();
	double R2 = cartesianPos.x * cartesianPos.x + cartesianPos.y * cartesianPos.y;
	double xxd = cartesianPos.x * this->x;
	double yyd = cartesianPos.y * this->y;
	double zzd = cartesianPos.z * this->z;
	return Vec3D((xxd+yyd+zzd)/r,
		(this->x* cartesianPos.y- cartesianPos.x*this->y)/R2,
		(cartesianPos.z*(xxd+yyd)-R2*this->z)/(pow(r,2)*sqrt(R2)));
}

Vec3D Vec3D::normalize(){
	double n = this->length();
	if (n == 0) return Vec3D();
	return Vec3D(x/n,y/n,z/n);
}

Vec3D Vec3D::crossProduct(const Vec3D* v1, const Vec3D* v2){
	return Vec3D(v1->y * v2->z - v1->z * v2->y, v1->z * v2->x - v1->x * v2->z, v1->x * v2->y - v1->y * v2->x);
}

std::string Vec3D::print(){
	return std::to_string(this->x) +','+ std::to_string(this->y) + ',' + std::to_string(this->z);
}

double Vec3D::distance(const Vec3D& a, const Vec3D& b){
	return sqrt(distance2(a,b));
	//return sqrt(pow(a->x-b->x,2)+ pow(a->y - b->y, 2)+ pow(a->z - b->z, 2));
}

double Vec3D::distance2(const Vec3D& a, const Vec3D& b){
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	double dz = a.z - b.z;
	return dx * dx + dy * dy + dz * dz;
}

void Vec3D::reset(){
	this->x = 0; this->y = 0; this->z = 0;
}

Vec3D operator/(const Vec3D& lhs, double& rhs){
	return Vec3D(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs);
}

Vec3D operator*(const Vec3D& lhs, const double& rhs){
	return Vec3D(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
}

Vec3D operator*(const double& lhs, const Vec3D& rhs){
	return rhs*lhs;
}

Vec3D operator+(Vec3D lhs, Vec3D const& rhs){
	lhs += rhs;
	return lhs;
}

Vec3D operator-(Vec3D lhs, Vec3D const& rhs){
	return Vec3D(lhs.x-rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

double operator*(const Vec3D& lhs, const Vec3D& rhs){
	return Vec3D::distance2(lhs, rhs);
}

double Vec3D::operator*(const Vec3D& rhs){
	return this->x * rhs.x + this->y * rhs.y + this->z * rhs.z;
}

Vec3D& Vec3D::operator+=(const Vec3D & rhs){
	this->x += rhs.x;
	this->y += rhs.y;
	this->z += rhs.z;
	return *this;
}

Vec3D operator*(float f, const Vec3D& v){
	return Vec3D(v.x * f, v.y * f, v.z * f);
}
