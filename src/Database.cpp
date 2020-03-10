#include "Database.h"

Database::Database(){
	this->isOpen = false;
	this->open();
}

bool Database::open(const char* name){
	if (std::strlen(name) == 0) {
		name = this->dBName;
	}
	else {
		this->dBName = name;
	}
	int rc = sqlite3_open(name, &db);
	if (rc){
		std::cerr << "Error opening SQLite3 database: " << sqlite3_errmsg(db) << std::endl << std::endl;
		sqlite3_close(db);
		this->isOpen = false;
		return true;
	}
	this->isOpen = true;
	return false;
}

bool Database::exec(char* sql){
	if (this->open()) //error during opening connection
		return true;
	char* zErrMsg = 0;
	int rc = sqlite3_exec(db, sql, NULL, 0, &zErrMsg);
	if (rc != SQLITE_OK) {
		fprintf(stderr, "SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
	}
	//sqlite3_close(db);
	return false;
}

void Database::setup(){
	if (this->open()) //error during opening connection
		return;
	char* sql = "CREATE TABLE IF NOT EXISTS simulation("
		"id INTEGER PRIMARY KEY,"
		"title TEXT NOT NULL,"
		"n_stars INTEGER NOT NULL,"
		"boxlength REAL NOT NULL,"
		"dt REAL NOT NULL,"
		"n_timesteps INTEGER NOT NULL);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS star("
		"id INTEGER NOT NULL,"
		"mass REAL NOT NULL,"
		"id_simulation INTEGER NOT NULL,"
		"PRIMARY KEY (id,id_simulation),"
		"FOREIGN KEY (id_simulation) "
			"REFERENCES simulation(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS velocity("
		"id INTEGER PRIMARY KEY,"
		"x REAL NOT NULL,"
		"y REAL NOT NULL,"
		"z REAL NOT NULL,"
		"timestep INTEGER NOT NULL,"
		"id_star INTEGER NOT NULL,"
		"id_simulation INTEGER NOT NULL,"
		"FOREIGN KEY (id_star,id_simulation) "
			"REFERENCES star(id,id_simulation) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS velocity2D("
		"id INTEGER PRIMARY KEY,"
		"x REAL NOT NULL,"
		"y REAL NOT NULL,"
		"z REAL NOT NULL,"
		"fk_velocity INTEGER NOT NULL,"
		"FOREIGN KEY (fk_velocity) "
			"REFERENCES velocity(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS position("
		"id INTEGER PRIMARY KEY,"
		"x REAL NOT NULL,"
		"y REAL NOT NULL,"
		"z REAL NOT NULL,"
		"timestep INTEGER NOT NULL,"
		"id_star INTEGER NOT NULL,"
		"id_simulation INTEGER NOT NULL,"
		"FOREIGN KEY (id_star,id_simulation) "
			"REFERENCES star(id,id_simulation) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS analysis("
		"fk_simulation INTEGER PRIMARY KEY,"
		"doEnergy INTEGER NOT NULL,"
		"doVelocity INTEGER NOT NULL,"
		"doVelocity2D INTEGER NOT NULL,"
		"FOREIGN KEY (fk_simulation) "
			"REFERENCES simulation(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS timeStepAnalysis("
		"id INTEGER PRIMARY KEY,"
		"averageVelocity REAL,"
		"kinE REAL,"
		"potE REAL,"
		"totE REAL,"
		"fk_analysis INTEGER NOT NULL,"
		"FOREIGN KEY (fk_analysis) "
			"REFERENCES analysis(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
}

int Database::insert(Parameters& parameters){
	if (!this->isOpen)
		this->open();
	std::string sql = "INSERT INTO simulation (n_stars,boxLength,dt,n_timesteps,title) VALUES (?1,?2,?3,?4,?5)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, parameters.getN_Stars());
	sqlite3_bind_double(st, 2, parameters.getBoxLength());
	sqlite3_bind_double(st, 3, parameters.getdt());
	sqlite3_bind_int(st, 4, parameters.getNTimesteps());
	char* cstr = new char[parameters.getTitle().length() + 1];
	std::strcpy(cstr, parameters.getTitle().c_str());
	sqlite3_bind_text(st, 5, cstr, -1, SQLITE_TRANSIENT);
	int returnCode = sqlite3_step(st);
	if (returnCode != SQLITE_DONE){
		throw "Could not insert new simulation";
	}
	sql = "SELECT Last_Insert_Rowid()";
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	returnCode = sqlite3_step(st);
	int simulationID = -1;
	if (returnCode == SQLITE_ROW){
		simulationID = sqlite3_column_int(st, 0);
	}
	sqlite3_finalize(st);
	//sqlite3_close(db);
	//parameters.setID(simulationID);
	return simulationID;
}

void Database::insertStars(int simulationID, std::vector<Star*>& stars, int timestep){
	if (!this->isOpen)
		this->open();
	//#pragma omp parallel for
	for (int i = 0; i < stars.size();++i) {
		insertStar(simulationID, stars.at(i), timestep);
	}
}

void Database::insertAnalysis(int simulationID, Analysis analysis){
	std::string sql = "INSERT INTO analysis (fk_simulation,doEnergy,doVelocity) VALUES (?1,?2,?3,?4)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, simulationID);
	sqlite3_bind_int(st, 2, analysis.getbEnergy());
	sqlite3_bind_int(st, 3, analysis.getbAverageVelocity());
	sqlite3_bind_int(st, 4, analysis.getbAverage2DVelocity());
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
}

void Database::timestep(int timestep, std::vector<Star*>& stars){
	if (!this->isOpen)
		this->open();
	#pragma omp parallel for
	for (int i = 0; i < stars.size(); ++i) {
		insertVelocity(stars.at(i)->id,stars.at(i)->velocity, timestep);
		insertPosition(stars.at(i)->id, stars.at(i)->position, timestep);
	}
}

void Database::insertStar(int simulationID, Star* star, int& timestep){
	if (!this->isOpen)
		this->open();
	std::string sql = "INSERT INTO star (id,mass,id_simulation) VALUES (?1,?2,?3)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, star->id);
	sqlite3_bind_double(st, 2, star->mass);
	sqlite3_bind_int(st, 3, simulationID);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
	insertVelocity(star->id, star->velocity, timestep);
	insertPosition(star->id, star->position, timestep);
}

void Database::insertPosition(int& idStar, Vec3D& position, int& timestep){
	if (!this->isOpen)
		this->open();
	std::string sql = "INSERT INTO position (x,y,z,id_star,timestep) VALUES (?1,?2,?3,?4,?5)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_double(st, 1, position.x);
	sqlite3_bind_double(st, 2, position.y);
	sqlite3_bind_double(st, 3, position.z);
	sqlite3_bind_int(st, 4, idStar);
	sqlite3_bind_int(st, 5, timestep);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
}

void Database::insertVelocity(int& idStar, Vec3D& velocity, int& timestep){
	if (!this->isOpen)
		this->open();
	std::string sql = "INSERT INTO velocity (x,y,z,id_star,timestep) VALUES (?1,?2,?3,?4,?5)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_double(st, 1, velocity.x);
	sqlite3_bind_double(st, 2, velocity.y);
	sqlite3_bind_double(st, 3, velocity.z);
	sqlite3_bind_int(st, 4, idStar);
	sqlite3_bind_int(st, 5, timestep);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
}

std::vector<Simulation> Database::selectSimulations(){
	std::vector<Simulation> simulations = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT id,title,n_stars,boxlength,dt,n_timesteps FROM simulation";
	sqlite3_stmt* stmt;
	if (sqlite3_prepare_v2(db, query.c_str(), query.size(), &stmt, nullptr) != SQLITE_OK) {
		// Error reporting and handling
	}
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		simulations.push_back(Simulation(sqlite3_column_int(stmt, 0), reinterpret_cast<const char*>(sqlite3_column_text(stmt,1)), sqlite3_column_int(stmt, 2),
			sqlite3_column_double(stmt, 3), sqlite3_column_double(stmt, 4), sqlite3_column_int(stmt, 5)));
	}
	sqlite3_finalize(stmt);
	return simulations;
}

std::vector<Vec3D> Database::selectVelocities(int timestep){
	std::vector<Vec3D> velocities = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT x,y,z FROM velocity where timestep = ?1";
	sqlite3_stmt* stmt;
	if (sqlite3_prepare_v2(db, query.c_str(), query.size(), &stmt, nullptr) != SQLITE_OK) {
		// Error reporting and handling
	}
	sqlite3_bind_int(stmt, 1, timestep);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		velocities.push_back(Vec3D(sqlite3_column_int(stmt, 0), sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2)));
	}
	sqlite3_finalize(stmt);
	return velocities;
}

std::vector<int> Database::selectTimesteps(){
	std::vector<int> timeSteps = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT DISTINCT timestep FROM velocity";
	sqlite3_stmt* stmt;
	if (sqlite3_prepare_v2(db, query.c_str(), query.size(), &stmt, nullptr) != SQLITE_OK) {
		// Error reporting and handling
	}
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		timeSteps.push_back(sqlite3_column_int(stmt, 0));
	}
	sqlite3_finalize(stmt);
	return timeSteps;
}

std::vector<Star*> Database::selectStars(int simulationID, int timeStep){
	std::vector<Star*> stars = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT star.id,mass,position.x,position.y,position.z,velocity.x,velocity.y,velocity.z " 
		"FROM star INNER JOIN velocity on velocity.fk_star = star.id"
		"INNER JOIN position on position.fk_star = star.id"
		"where position.timestep = ?1 AND velocity.timestep = ?1 AND star.fk_simulation = ?2";
	sqlite3_stmt* stmt;
	if (sqlite3_prepare_v2(db, query.c_str(), query.size(), &stmt, nullptr) != SQLITE_OK) {
		// Error reporting and handling
	}
	sqlite3_bind_int(stmt, 1, timeStep);
	sqlite3_bind_int(stmt, 2, simulationID);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		stars.push_back(new Star(sqlite3_column_int(stmt, 0),sqlite3_column_double(stmt,1), 
			sqlite3_column_double(stmt, 2), sqlite3_column_double(stmt, 3), sqlite3_column_double(stmt, 4),
			sqlite3_column_double(stmt, 5), sqlite3_column_double(stmt, 6), sqlite3_column_double(stmt, 7)));
	}
	sqlite3_finalize(stmt);
	return stars;
}
