#include "Simulation.h"

Simulation::Simulation(int id, Database* database):SimulationData(id){}

Simulation::Simulation(int id):SimulationData(id) {}

void Simulation::setID(int id){
	this->simulationID = id;
}

int Simulation::getID(){
	return this->simulationID;
}

//int Simulation::getNStars() {
//	return n_Stars;
//}
//
//double Simulation::getdt() {
//	return dt;
//}
//
//int Simulation::getNTimesteps() {
//	return nTimesteps;
//}

std::string Simulation::print(){
	return "ID: " + std::to_string(this->simulationID) +
		" | Title: " + this->title + 
		" #Stars: " + std::to_string(this->n_Stars) +'\n';
}

void Simulation::run(){
	int nextStarIndex = database->selectLastID("star") + 1;
	//Init stars
	InitialConditions initialConditions = InitialConditions(this);
	std::vector<Star*> stars = initialConditions.initStars(nextStarIndex);
	double totalMass = initialConditions.initialMass(stars);
	initialConditions.plummerSphere(stars, 1, totalMass);
	database->insertStars(this->getID(), stars, 0);

	//Integrate
	Integrator rk4 = Integrator(this->getdt());
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	for (int i = 1; i < this->getNTimesteps(); i++) {

		Vec3D tlf = Vec3D(), brb = Vec3D();
		Node::findCorners(tlf, brb, stars);
		Node root = Node(tlf, brb, nullptr, this);
		for (Star* star : stars) {
			root.insert(star);
		}

		root.calculateMassDistribution();

		rk4.euler(stars, &root);

		if (i % 10 == 0) {
			//InOut::writeWithLabel(stars, "./Output/stars" + std::to_string(i) + ".dat");
			//InOut::writeAll(stars, "./Output/stars_all" + std::to_string(i) + ".dat");
			database->timestep(i, stars);
		}
	}
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	//InOut::write(stars,"stars.dat");
	//InOut::write(&root);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
	std::cout << "done" << std::endl;
}
