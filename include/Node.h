/**
 * Implementation of Barnes-Hut simulation. Ie octree, particles (=Stars) and gravitational force calculations.
 *
 * Each star is assigned a Node/Cube. The star is located within the cube. Starting at the root node each node has up to 8 child nodes.
 * ie the space within the parent cube is split into 8 sub-cubes. If a cube already has a star assigned and a new star gets added within that area
 * the new star aswell as the star already located within the cube get placed within one sub cube respectivly.
 * 
 * Why go through all this hassle? In N-Body simulation the force acting on each Particle depends on the location of all other N-1 particles.
 * But if the distance of the current particle to the center of mass of one of the Cubes/Nodes is large enough relative to the size of the cube.
 * then the center of mass can be used to calculate the force acting upon the particle rather than all the particles contained within that cube.
 * Runtime is N*log(N), which is a big improvmenet over the bruteforce method (order N^2)
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#pragma once
#include <math.h> // pow
#include <iostream>

#include "Structs.h"
#include "Vec3D.h"
#include "Star.h"
#include "Parameters.h"

class Node {
	/** @brief location of the center of the Node. Used to find the octant a star in question lies within. */
	Vec3D center;
	/** @brief The parent of the node. Only for the root node this value is 0 */
	Node *parent;
	/** @brief Amount of stars inside node. For the root node this is equal to the total amount of stars and for a leaf node it is equal to 1 */
	int n_stars;
	/** @brief The star assigned to the node. May be 0. */
	Star *star;
	/** @brief Center of mass. If there is only one star within the node area this equals the stars location */
	Vec3D centerOfMass;
	/** @brief The total mass located within the Nodes area. If there is only one star within the node area this equals the stars mass. */
	double mass;
	/** @brief Defines maximum distance between start and center of mass for force to be approximated, static to save memory. Usualy chosen to be ~1 */
	static double precission; 
	/** @brief Gravitational constant in astronomical units: parsec*solar mass^-1*km^2/s^2. Here needed for force calculation.*/
	static double G;
	/** @brief softening parameter for force calculation*/
	static double softening;
	/** @brief 1 km divided by 1 pc */
	static double kmInpc;
public: //variables
	/** @brief Pointer array to child nodes where indexing according to ::Octant. @note May contain null_ptrs */
	Node* children[8] = {};
	/** @brief Coordinates of the top left front corner of the cube */
	Vec3D top_left_front;
	/** @brief Coordinates of the bottom right back corner of the cube */
	Vec3D bottom_right_back;
	/** @brief If a node has a child node this value is true. Redundant information to save cpu cycles. */
	bool internalNode;
private: //functions
	Node(Vec3D top_left_front, Vec3D bottom_right_back, Node* parent, double G, double softening, double precission);
public:
	Node();
	/** 
	@param top_left_front,  bottom_right_back corners defining the cell.
	@param parent Pointer to the parent of the node. Use null nullptr to create a root node.
	@param parameters Parameters defined in simulation.cfg file. Member variables G, softening and precission are set via parameters.
	*/
	Node(Vec3D top_left_front, Vec3D bottom_right_back, Node* parent, Parameters* parameters);
	~Node();
	/**
	@static
	@brief call this function to get the parameters for the root node.
	@param [in,out] tlf, brb Vectors which store the corners of the root node after the call
	@param stars Stars which are within the cube volume defined by tlf and brb.
	*/
	static void findCorners(Vec3D& tlf, Vec3D& brb, std::vector<Star*>&stars);
	/**
	@static
	@brief Inserts a star into the octree. This always creates at least one new node in the tree.
	@param star the star which is inserted into the tree.
	*/
	void insert(Star* star);
	/**
	@brief Creates the sub node coresponding to the given octant.
	@param octant Corner of the current node where the sub node will be located.
	@return Pointer to the newly created sub node
	*/
	Node* create(Octant octant);
	/**
	@brief The octant of the star is determined, Based on the location of the star relative to the center of the node.
	@param star The star of which the octand will be determined.
	@return Octant (subcell within the current cell) where the star lies within.
	*/
	Octant getOctant(Star* star);
	/**
	@brief determines wether or not the node is the root node of the octree.
	@return true if the caller is root, otherwise false.
	*/
	bool isRoot();
	/**
	@return Outputs the corners (top_left_front,  bottom_right_back) defining the node. Format: Vec3D::print(),Vec3D::print()
	*/
	std::string print();
	/**
	@brief Calculates and sets the center of mass vector and the mass of the calling and all its child nodes
	@note Call by root node to update the whole tree
	*/
	void calculateMassDistribution();
	/**
	@brief Update acceleration of the star based on the force applied to the star from the current node if the approximation is good enough. 
	Otherwise the function is called by all child nodes.
	@param star Acceleration of this star will be updated according to the potential of the stars contained within the node.
	@note Should be called by root node only. Star does not have to be "part" of the tree.
	*/
	void applyForce(Star* star);
	/**
	@brief Sets acceleration based on the force present at position from the current node if the approximation is good enough.
	Otherwise the function is called by all child nodes.
	@param acceleration Acceleration at this position will be updated according to the potential of the stars contained within the node.
	@param position position at which acceleration is calucated according to the potential of the stars contained within the node.
	@note Should be called by root node only. Star does not have to be "part" of the tree.
	*/
	void Node::applyForce(const Vec3D position, Vec3D& acceleration);
};

