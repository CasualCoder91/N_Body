#include "Simulation.h"

Simulation::Simulation(int id, Database* database):SimulationData(id){
	this->database = database;
}

Simulation::Simulation(int id):SimulationData(id) {}

Simulation::Simulation(int id, SimulationData* simulationData){
	this->simulationID = id;
}

void Simulation::setID(int id){
	this->simulationID = id;
}

int Simulation::getID(){
	return this->simulationID;
}

std::string Simulation::print(){
	return "ID: " + std::to_string(this->simulationID) +
		" | Title: " + this->title + 
		" #Stars: " + std::to_string(this->nStars) +'\n';
}

void Simulation::run(){
	//Init stars
	InitialConditions initialConditions = InitialConditions(this);

	//Init clusterStars
	int nextStarIndex = database->selectLastID("star") + 1;
	std::vector<Star*> clusterStars = initialConditions.initStars(nextStarIndex);
	double totalMass = initialConditions.initialMassSalpeter(clusterStars,0.08,100);
	initialConditions.plummerSphere(clusterStars, 1, totalMass);
	Vec3D offset = Vec3D(5000, 5000, 0);
	initialConditions.offsetCluster(clusterStars, offset);
	Vec3D clusterVelocity = Vec3D();
	initialConditions.sampleDiskVelocity(clusterVelocity, offset);
	for (Star* star : clusterStars) {
		star->velocity += clusterVelocity;
	}
	database->insertStars(this->getID(), clusterStars, 0);

	//Init otherStars
	std::vector<Star*> otherStars = initialConditions.massDisk(20); //test
	initialConditions.sampleDiskPositions(otherStars, Vec3D(-30000, -30000, -6000), Vec3D(60000, 60000, 12000)); //test
	initialConditions.sampleDiskVelocities(otherStars);

	//Integrate
	Integrator integrator = Integrator(this->getdt());
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	Vec3D tlf = Vec3D(), brb = Vec3D();
	for (int i = 1; i < this->getNTimesteps(); i++) {

		if (clusterStars.size() > 0) {
			Node::findCorners(tlf, brb, clusterStars);
			Node root = Node(tlf, brb, nullptr, this); // Cleanup/Destructor of tree needed in every timestep
			for (Star* star : clusterStars) {
				root.insert(star);
			}
			root.calculateMassDistribution();

			//Force clusterStars
			for (int i = 0; i < clusterStars.size(); ++i) {
				clusterStars.at(i)->acceleration.reset();
				Potential::applyForce(clusterStars.at(i));
				root.applyForce(clusterStars.at(i)->position, &clusterStars.at(i)->acceleration);
			}
		}

		//Force otherStars
		for (int i = 0; i < otherStars.size(); ++i) {
			otherStars.at(i)->acceleration.reset();
			Potential::applyForce(otherStars.at(i));
		}

		if(clusterStars.size()>0)
			integrator.euler(clusterStars);
		integrator.euler(otherStars);

		if (i % this->getOutputTimestep() == 0) {
			InOut::writeWithLabel(otherStars, "Stars/3/otherStars" + std::to_string(i) + ".dat");
			//InOut::writeAll(clusterStars, "./Output/stars_all" + std::to_string(i) + ".dat");
			//database->timestep(i, clusterStars);
			//std::cout << i << std::endl;
		}
	}
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	//InOut::write(clusterStars,"stars.dat");
	//InOut::write(&root);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
	std::cout << "done" << std::endl;
}
