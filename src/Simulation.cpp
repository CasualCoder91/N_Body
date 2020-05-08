#include "Simulation.h"

Simulation::Simulation(int id, Database* database){
	this->database = database;
	this->simulationID = id;
	this->potential = new Potential(this);
}

Simulation::Simulation(int id){
	this->simulationID = id;
	this->database = new Database();
	this->potential = new Potential(this);
}

Simulation::Simulation(int id, SimulationData* simulationData){
	this->simulationID = id;
	this->potential = new Potential(this);
}

Simulation::Simulation(int id, Database* database, Parameters* parameters) :Parameters{ parameters } {
	this->simulationID = id;
	this->database = database;
	this->potential = new Potential(this);
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
	InitialConditions initialConditions = InitialConditions(this,Simulation::potential);

	//Init clusterStars
	int nextStarIndex = database->selectLastID("star") + 1;
	std::vector<Star*> clusterStars = initialConditions.initStars(nextStarIndex);
	double totalMass = initialConditions.initialMassSalpeter(clusterStars,0.08,100);
	initialConditions.plummerSphere(clusterStars, 1, totalMass);
	Vec3D offset = Vec3D(40, 40, 0);
	initialConditions.offsetCluster(clusterStars, offset);
	Vec3D clusterVelocity = Vec3D();
	initialConditions.sampleDiskVelocity(clusterVelocity, offset);
	for (Star* star : clusterStars) {
		star->velocity += clusterVelocity;
	}
	database->insertStars(this->getID(), clusterStars, 0);

	//std::vector<Star*> fieldStars = initialConditions.massDisk(30); //test
	////Init fieldStars
	//if (fieldStars.size() > 0) {
	//	initialConditions.sampleDiskPositions(fieldStars, Vec3D(-30000, -30000, -6000), Vec3D(60000, 60000, 12000)); //test
	//	initialConditions.sampleDiskVelocities(fieldStars);
	//}
	std::vector<Star*> fieldStars = initialConditions.initFieldStars(nextStarIndex); //0.00029088833
	database->insertStars(this->getID(), fieldStars, 0);


	//Integrate
	Integrator integrator = Integrator(this->getdt());
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	Vec3D tlf = Vec3D(), brb = Vec3D();
	for (int i = 0; i < this->getNTimesteps(); i++) {
		ProgressBar progressBar = ProgressBar(0, this->getNTimesteps());
		if (i % this->getOutputTimestep() == 0) {
			//InOut::writeWithLabel(fieldStars, "Stars/5/fieldStars" + std::to_string(i) + ".dat");
			//InOut::writeWithLabel(clusterStars, "Stars/4/clusterStars" + std::to_string(i) + ".dat");
			//InOut::writeAll(clusterStars, "./Output/stars_all" + std::to_string(i) + ".dat");
			database->timestep(i, clusterStars);
			database->timestep(i, fieldStars);
			progressBar.Update(i);
			progressBar.Print();
			//std::cout << i << std::endl;
		}

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
				potential->applyForce(clusterStars.at(i));
				root.applyForce(clusterStars.at(i)->position, &clusterStars.at(i)->acceleration);
			}
		}

		//Force fieldStars
		if (fieldStars.size() > 0) {
			for (int i = 0; i < fieldStars.size(); ++i) {
				fieldStars.at(i)->acceleration.reset();
				potential->applyForce(fieldStars.at(i));
			}
		}

		if(clusterStars.size()>0)
			integrator.euler(clusterStars);
		if (fieldStars.size() > 0) {
			integrator.euler(fieldStars);
		}

	}
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	//InOut::write(clusterStars,"stars.dat");
	//InOut::write(&root);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
	std::cout << "done" << std::endl;
}

Simulation::~Simulation(){
	delete potential;
}
