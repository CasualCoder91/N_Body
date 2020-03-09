#include "Database.h"

Database::Database(){
	Database::open();
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
		return true;
	}
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
	sqlite3_close(db);
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
		"id INTEGER PRIMARY KEY,"
		"mass REAL NOT NULL,"
		"fk_simulation INTEGER NOT NULL,"
		"FOREIGN KEY (fk_simulation) "
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
		"fk_star INTEGER NOT NULL,"
		"FOREIGN KEY (fk_star) "
			"REFERENCES star(id) "
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
		"fk_star INTEGER NOT NULL,"
		"FOREIGN KEY (fk_star) "
			"REFERENCES star(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS analysis("
		"id INTEGER PRIMARY KEY,"
		"fk_simulation INTEGER NOT NULL,"
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

int Database::insert(Parameters parameters){
	int simulationID = -1;
	std::string sql = "INSERT INTO simulation (n_stars,boxLength,dt,n_timesteps,title) VALUES (?1,?2,?3,?4,?5)";
	if (!this->open()){
		sqlite3_stmt* st;
		sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
		sqlite3_bind_int(st, 1, parameters.getN_Stars());
		sqlite3_bind_double(st, 2, parameters.getBoxLength());
		sqlite3_bind_double(st, 3, parameters.getdt());
		sqlite3_bind_int(st, 4, parameters.getNTimesteps());
		sqlite3_bind_text(st, 5, parameters.getTitle().c_str(), parameters.getTitle().length()+1, NULL);
		int returnCode = sqlite3_step(st);
		if (returnCode != SQLITE_DONE){
			throw "Could not insert new simulation";
		}
		sql = "SELECT Last_Insert_Rowid()";
		sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
		returnCode = sqlite3_step(st);
		if (returnCode == SQLITE_ROW){
			simulationID = sqlite3_column_int(st, 0);
		}
	}
	sqlite3_close(db);
	return simulationID;
}

void Database::insert(int simulationID, std::vector<Star*>& stars){
	sqlite3_stmt* st;
	if (!this->open()) {
		std::string sql; 
		#pragma omp parallel for
		for (int i = 0; i < stars.size();++i) {
			sql = "INSERT INTO star (id,mass,fk_simulation) VALUES (?1,?2,?3)";
			sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
			sqlite3_bind_int(st, 1, stars.at(i)->id);
			sqlite3_bind_double(st, 2, stars.at(i)->mass);
			sqlite3_bind_int(st, 3, simulationID);
			int returnCode = sqlite3_step(st);
			sql = "INSERT INTO velocity (x,y,z,fk_star,timestep) VALUES (?1,?2,?3,?4,?5)";
			sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
			sqlite3_bind_double(st, 1, stars.at(i)->velocity.x);
			sqlite3_bind_double(st, 2, stars.at(i)->velocity.y);
			sqlite3_bind_double(st, 3, stars.at(i)->velocity.z);
			sqlite3_bind_int(st, 4, stars.at(i)->id);
			sqlite3_bind_int(st, 5, 0);
			returnCode = sqlite3_step(st);
			sql = "INSERT INTO position (x,y,z,fk_star,timestep) VALUES (?1,?2,?3,?4,?5)";
			sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
			sqlite3_bind_double(st, 1, stars.at(i)->position.x);
			sqlite3_bind_double(st, 2, stars.at(i)->position.y);
			sqlite3_bind_double(st, 3, stars.at(i)->position.z);
			sqlite3_bind_int(st, 4, stars.at(i)->id);
			sqlite3_bind_int(st, 5, 0);
			returnCode = sqlite3_step(st);
		}
		sqlite3_close(db);
	}
}
