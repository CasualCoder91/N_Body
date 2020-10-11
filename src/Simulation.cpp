#include "Simulation.h"

Simulation::Simulation(int id){
	this->simulationID = id;
	this->database = new Database();
	this->potential = new MWPotential(this);
}

Simulation::Simulation(int id, Database* database, SimulationData* simulationData){
	this->simulationID = id;
	this->potential = new MWPotential(this);
	this->database = database;
}

Simulation::Simulation(int id, Database* database, Parameters* parameters) :Parameters{ parameters } {
	this->simulationID = id;
	this->database = database;
	this->potential = new MWPotential(this);
}

void Simulation::setID(int id){
	this->simulationID = id;
}

int Simulation::getID(){
	return this->simulationID;
}

void Simulation::run(){

	std::cout << "boxlength: " << Constants::boxLength << std::endl;
	std::cin.get();
	//Init stars
	InitialConditions initialConditions = InitialConditions(Simulation::potential);

	//Init clusterStars
	int nextStarIndex = database->selectLastID("star") + 1;

	//std::vector<Star*> clusterStars = initialConditions.initStars(nextStarIndex,nStars);
	//double totalMass = initialConditions.initialMassSalpeter(clusterStars,minMass,maxMass);
	//initialConditions.plummerSphere(clusterStars, totalMass,boxLength,G);

	std::vector<Star*> clusterStars = InOut::readMcLuster(nextStarIndex, "Data/eintest.txt");
	initialConditions.offsetCluster(clusterStars, clusterLocation);

	double circVel = potential->circularVelocity(&clusterLocation);
	Vec3D clusterVelocity = Vec3D(0,-220,10);
	//initialConditions.sampleDiskVelocity(clusterVelocity, clusterLocation);
	for (Star* star : clusterStars) {
		star->velocity += clusterVelocity;
	}
	database->insertStars(this->getID(), clusterStars, 0,true);
	nextStarIndex = database->selectLastID("star") + 1;

	std::vector<Star*> fieldStars = initialConditions.initFieldStars(nextStarIndex,focus,viewPoint,distance,dx,angleOfView); //0.00029088833
	database->insertStars(this->getID(), fieldStars, 0,false);

	std::cout << "Starting integration" << std::endl;
	//Integrate
	Integrator integrator = Integrator(this->getdt());
	std::cout << "dt = " << integrator.dt << std::endl;
	//std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	Vec3D tlf = Vec3D(), brb = Vec3D();
	ProgressBar progressBar = ProgressBar(0, this->getNTimesteps());

	std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();;
	std::chrono::high_resolution_clock::time_point endTime;

	for (int i = 0; i <= this->getNTimesteps(); i++) {
		if (i>0 && i % this->outputTimestep == 0) {
			database->timestep(i, clusterStars);
			database->timestep(i, fieldStars);
			progressBar.Update(i);
			progressBar.Print();
		}

		if (clusterStars.size() > 0) {
			Node::findCorners(tlf, brb, clusterStars);
			Node root = Node(tlf, brb, nullptr, this); // Cleanup/Destructor of tree needed in every timestep

			for (Star* star : clusterStars) {
				root.insert(star);
			}
			root.calculateMassDistribution();

			//Force clusterStars -- use if Euler
			//for (int j = 0; j < clusterStars.size(); ++j) {
			//	clusterStars[j]->acceleration.reset();
			//	potential->applyForce(clusterStars[j]);
			//	root.applyForce(clusterStars[j]->position, clusterStars[j]->acceleration);
			//}
			//integrator.euler(clusterStars);

			integrator.Leapfrog(clusterStars, &root, this->potential);
		}

		//Force fieldStars -- use if Euler
		//if (fieldStars.size() > 0) {
		//	for (int i = 0; i < fieldStars.size(); ++i) {
		//		fieldStars[i]->acceleration.reset();
		//		potential->applyForce(fieldStars[i]);
		//	}
		//}

		if (fieldStars.size() > 0) {
			//integrator.euler(fieldStars);
			integrator.Leapfrog(fieldStars, this->potential);
		}

	}
	endTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	//InOut::write(clusterStars,"stars.dat");
	//InOut::write(&root);
	std::cout << "Time needed for integration: " << time_span.count() << "seconds" << std::endl;
	std::cout << "done" << std::endl;
}

Simulation::~Simulation(){
	delete potential;
}
