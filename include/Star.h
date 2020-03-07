/**
 * Star containing information about mass, position, velocity and acceleration.
 *
 * @author Alarich Herzner
 * @version 0.9 02.03.2020
*/

#pragma once
#include "Vec3D.h"
class Star{
public:
	Star();
	Star(double mass);
	Star(double mass, Vec3D position);
	Star(double mass, double x, double y, double z);
	/** @brief Mass of the star in solar mass units */
	double mass;
	/** @brief Position of the star in parsec */
	Vec3D position;
	/** @brief Velocity of the star in km/s */
	Vec3D velocity;
	/** @brief Acceleration of the star in km^2/(parsec*s^2)
	@note Reason for this choice of units: dt is in parsec*s/km*/
	Vec3D acceleration;
	/**
	 @brief creates string containing all member variables.
	 @return created string
	 */
	std::string dump();
	/**
	 @brief set all member variables to 0.
	 @return created string
	 */
	void reset();
};

