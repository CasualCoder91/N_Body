#pragma once
#include <math.h> // pow
#include <iostream>

#include "Structs.h"
#include "Vec3D.h"
#include "Star.h"
#include "Parameters.h"

class Node : Parameters {
	Vec3D center;
	Node *parent;
	int n_stars; // amount of stars inside node
	Star *star;
	Vec3D centerOfMass;
	double mass;
	static double precission; // defines maximum distance between start and center of mass for force to be approximated, static to save memory
public: //variables
	Node* children[8] = {};//pointer array to child nodes where indexing according to enum Octant
	Vec3D top_left_front;
	Vec3D bottom_right_back;
	bool internalNode;
public:
	Node(Vec3D top_left_front, Vec3D bottom_right_back, Node* parent);
	~Node();

	static void findCorners(Vec3D& tlf, Vec3D& brb, std::vector<Star*>&stars);
	void insert(Star* star);
	Node* create(Octant octant);
	Octant getOctant(Star* star);
	bool isRoot();
	std::string print();
	void calculateMassDistribution();
	void applyForce(Star* star); //Update force applied to star from Node
};

