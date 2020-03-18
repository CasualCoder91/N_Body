#include "Node.h"

double Node::precission = 0.5;
double Node::G = 4.483e-3;
double Node::softening = 0.16;

Node::Node(Vec3D top_left_front, Vec3D bottom_right_back, Node* parent, SimulationData* parameters){
	this->top_left_front = top_left_front;
	this->bottom_right_back = bottom_right_back;
	this->center = Vec3D::center(top_left_front, bottom_right_back);
	this->parent = parent;
	this->G = parameters->getG();
	this->softening = parameters->getSoftening();
	this->precission = parameters->getPrecission();
	this->n_stars = 0;
	this->star = nullptr; //needed?
	this->centerOfMass = Vec3D();
	this->mass = 0;
	this->internalNode = false;
}

Node::Node(Vec3D top_left_front, Vec3D bottom_right_back, Node* parent, double G, double softening, double precission){
	this->top_left_front = top_left_front;
	this->bottom_right_back = bottom_right_back;
	this->center = Vec3D::center(top_left_front, bottom_right_back);
	this->parent = parent;
	this->G = G;
	this->softening = softening;
	this->precission = precission;
	this->n_stars = 0;
	this->star = nullptr; //needed?
	this->centerOfMass = Vec3D();
	this->mass = 0;
	this->internalNode = false;
}

Node::~Node(){
	for(Node* child : children){
		delete child;
	}
	//no NOT delete star pointer ^^
}

void Node::findCorners(Vec3D& tlf, Vec3D& brb, std::vector<Star*>& stars){
	tlf = stars.at(0)->position;
	brb = stars.at(0)->position;
	for (Star* star : stars) {
		if (star->position.x > brb.x)
			brb.x = star->position.x;
		else if (star->position.x < tlf.x)
			tlf.x = star->position.x;
		if (star->position.y > tlf.y)
			tlf.y = star->position.y;
		else if (star->position.y < brb.y)
			brb.y = star->position.y;
		if (star->position.z > brb.z)
			brb.z = star->position.z;
		else if (star->position.z < tlf.z)
			tlf.z = star->position.z;
	}
	//Add padding between box and stars
	tlf.x -= abs(tlf.x)*0.1;
	tlf.y += abs(tlf.y) * 0.1;
	tlf.z -= abs(tlf.z) * 0.1;
	brb.x += abs(brb.x) * 0.1;
	brb.y -= abs(brb.y) * 0.1;
	brb.z += abs(brb.z) * 0.1;
}

void Node::insert(Star* star)
{
	if (this->n_stars > 1) {
		Octant octant = this->getOctant(star);
		if (!this->children[static_cast<int>(octant)]) {//if child node does not exist yet create it
			this->children[static_cast<int>(octant)] = this->create(octant);
		}
		this->children[static_cast<int>(octant)]->insert(star);
		this->internalNode = true;
	}
	else if (this->n_stars == 1) {
		Octant octant = this->getOctant(star);
		if (!this->children[static_cast<int>(octant)]) {//if child node does not exist yet create it
			this->children[static_cast<int>(octant)] = this->create(octant);
		}
		this->children[static_cast<int>(octant)]->insert(star);
		octant = this->getOctant(this->star);
		if (!this->children[static_cast<int>(octant)]) {//if child node does not exist yet create it
			this->children[static_cast<int>(octant)] = this->create(octant);
		}
		this->children[static_cast<int>(octant)]->insert(this->star);
		if(this->star)
			this->star=nullptr;
		this->internalNode = true;
	}
	else {
		this->star = star;
		this->mass = star->mass;
		this->centerOfMass = star->position;
	}
	this->n_stars++;
}

Node* Node::create(Octant octant){
	switch(static_cast<int>(octant)) {
		case 0: return new Node(top_left_front, center, this, G, softening, precission);
		case 1: return new Node(Vec3D(center.x, top_left_front.y, top_left_front.z), Vec3D(bottom_right_back.x, center.y, center.z),this, G, softening, precission);
		case 2: return new Node(Vec3D(center.x, center.y, top_left_front.z), Vec3D(bottom_right_back.x,bottom_right_back.y,center.z), this, G, softening, precission);
		case 3: return new Node(Vec3D(top_left_front.x, center.y, top_left_front.z), Vec3D(center.x, bottom_right_back.y, center.z), this, G, softening, precission);
		case 4: return new Node(Vec3D(top_left_front.x, top_left_front.y, center.z), Vec3D(center.x, center.y, bottom_right_back.z), this, G, softening, precission);
		case 5: return new Node(Vec3D(center.x, top_left_front.y, center.z), Vec3D(bottom_right_back.x, center.y, bottom_right_back.z), this, G, softening, precission);
		case 6: return new Node(center, bottom_right_back, this, G, softening, precission);
		case 7: return new Node(Vec3D(top_left_front.x, center.y, center.z), Vec3D(center.x, bottom_right_back.y, bottom_right_back.z), this, G, softening, precission);
	}
	return nullptr;
}

Octant Node::getOctant(Star* star)
{
	//if (star->position.x<top_left_front.x || star->position.x>bottom_right_back.x ||
	//	star->position.y<bottom_right_back.y || star->position.y>top_left_front.y ||
	//	star->position.z<top_left_front.z || star->position.z>bottom_right_back.z) {
	//	return Octant::Invalid;//error: no Octant found
	//}
	bool top = star->position.y > center.y;
	bool front = star->position.z < center.z;
	bool left = star->position.x < center.x;
	if (front) {//front
		if (top) {//top
			if (left) {//left
				return Octant::top_left_front;
			}
			else{//right
				return Octant::top_right_front;
			}
		}
		else{//bottom
			if (left) {//left
				return Octant::bottom_left_front;
			}
			else{//right
				return Octant::bottom_right_front;
			}
		}
	}
	else{ //back
		if (top) {//top
			if (left) {//left
				return Octant::top_left_back;
			}
			else{//right
				return Octant::top_right_back;
			}
		}
		else{//bottom
			if (left) {//left
				return Octant::bottom_left_back;
			}
			else{//right
				return Octant::bottom_right_back;
			}
		}
	}
	return Octant::Invalid;//error: no Octant found
}

bool Node::isRoot(){
	return parent==nullptr;
}

std::string Node::print(){
	return this->top_left_front.print() + ',' + this->bottom_right_back.print() + '\n';
}

void Node::calculateMassDistribution(){
	if (this->internalNode) {
		this->mass = 0;
		this->centerOfMass = Vec3D();
		for (Node* child : this->children) {
			if (child) {
				child->calculateMassDistribution();
				this->mass += child->mass;
				this->centerOfMass = centerOfMass + child->centerOfMass*child->mass;
			}
		}
		this->centerOfMass = this->centerOfMass / this->mass;
	}
	else {
		this->centerOfMass = this->star->position;
	}
	//std::cout << "Node: " << this << std::endl;
	//std::cout << "Star: " << this->star << std::endl;
	//std::cout << star->Dump();
	//std::cout << "Center of Mass: " + this->centerOfMass.print() << std::endl<<std::endl;
}

void Node::applyForce(Star* star){
	Node::applyForce(star->position, &star->acceleration);
}

void Node::applyForce(const Vec3D position,Vec3D* acceleration){
	//if (!this->isRoot())
	//	throw "Only root may call this one.";
	double temp = 0;
	if (!this->internalNode) { //&& this->star != star
		//attractive force -> pointing towars center of mass
		double dx = this->star->position.x - position.x;
		double dy = this->star->position.y - position.y;
		double dz = this->star->position.z - position.z;
		temp = dx * dx + dy * dy + dz * dz;
		if (temp > 0) {
			temp = this->G * this->star->mass / pow(temp + this->softening * this->softening, 3 / 2);
			//star->acceleration += Vec3D(temp * dx, temp * dy, temp * dz);
			acceleration->x += temp * dx;
			acceleration->y += temp * dy;
			acceleration->z += temp * dz;
		}
		return;
	}
	double distStarCOM = Vec3D::distance(&(position), &(this->centerOfMass));
	if ((this->bottom_right_back.x - this->top_left_front.x) / distStarCOM < precission && this->star != star && distStarCOM>0) {
		temp = this->G * this->mass / pow(distStarCOM + softening, 3);
		//star->acceleration += Vec3D(temp * (star->position.x - this->centerOfMass.x), temp * (star->position.y - this->centerOfMass.y), temp * (star->position.z - this->centerOfMass.z));
		acceleration->x += temp * (this->centerOfMass.x - position.x);
		acceleration->y += temp * (this->centerOfMass.y - position.y);
		acceleration->z += temp * (this->centerOfMass.z - position.z);
	}
	else {
		#pragma omp parallel for
		for (int i = 0; i < 8;++i) {
			if (this->children[i]) {
				this->children[i]->applyForce(position, acceleration);
			}
		}
	}
}

