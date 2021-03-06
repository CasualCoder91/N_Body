#include "InOut.h"

std::string InOut::makeDirectory(std::string path){
	setOutputDirectory(path);
	int nError = 0;
#if defined(_WIN32)
	nError = _mkdir(path.c_str()); // can be used on Windows
#else 
	mode_t nMode = 0733; // UNIX style permissions
	nError = mkdir(sPath.c_str(), nMode); // can be used on non-Windows
#endif
	if (nError != 0) {
		int i = 0;
		// handle your error here
	}
	return path;
}

void InOut::setOutputDirectory(std::string& filename){
	size_t found = filename.find('/');
	if (found == std::string::npos)
		filename = outputDirectory + filename;
}

void InOut::write(std::vector<Star*> stars, std::string filename) {
	setOutputDirectory(filename);
	std::ofstream file(filename);
	for (Star* star : stars) {
		file << star->position.print() << '\n';
	}
	file.close();
}

void InOut::writeWithLabel(std::vector<Star*> stars, std::string filename){
	setOutputDirectory(filename);
	std::ofstream file(filename);
	for (int i = 0; i < stars.size();++i) {
		file << stars[i]->position.x << ',' << stars[i]->position.y << ',' << stars[i]->position.z << ',' << i << '\n';
	}
	file.close();
}

void InOut::writeAll(std::vector<Star*> stars, std::string filename){
	setOutputDirectory(filename);
	std::ofstream file(filename);
	#pragma omp parallel for
	for (int i = 0; i < stars.size(); ++i) {
		file << "Star: " << i << '\n';
		file << stars[i]->dump() << '\n';
	}
	file.close();
}

void InOut::write(Node* tree, std::string filename){
	if (!tree->isRoot()) {
		throw "Write function may only be called on root node";
	}
	else {
		setOutputDirectory(filename);
		std::ofstream file(filename);
		writeRecursively(&file, tree);
		file.close();
	}
	return;
}



void InOut::writeRecursively(std::ofstream* file_ptr,Node* node_ptr) {
	for (Node* child : node_ptr->children) {
		if (child) {
			InOut::writeRecursively(file_ptr, child);
		}
	}
	*file_ptr << node_ptr->print();
	return;
}

std::vector<Vec3D> InOut::readVectors(std::string filename){
	std::vector<Vec3D> vectors;
	std::string line;
	std::string delimiter = ",";
	std::ifstream sVectors(filename);
	if (sVectors.is_open()){
		while (std::getline(sVectors, line)){
			std::string firstToken = line.substr(0, line.find(delimiter));
			if (checkIsDouble(firstToken)) {
				size_t pos = line.find(delimiter);
				double x = std::stod(line.substr(0, pos));
				line.erase(0, pos + delimiter.length());
				double y = std::stod(line.substr(0, pos));
				line.erase(0, pos + delimiter.length());
				double z = std::stod(line, nullptr);
				vectors.push_back(Vec3D(x, y, z));
			}
		}
		sVectors.close();
	}
	return vectors;
}

bool InOut::checkIsDouble(std::string inputString){
	char* end;
	double result = strtod(inputString.c_str(), &end);
	if (end == inputString.c_str() || *end != '\0') 
		return false;
	return true;
}

std::vector<Star*> InOut::readMcLuster(int firstID, std::string filename){
	std::vector<Star*> stars;
	std::string line;
	std::string delimiter = "\t";
	std::ifstream sStars(filename);
	int id = firstID;
	if (sStars.is_open()) {
		while (std::getline(sStars, line)) {
			std::string firstToken = line.substr(0, line.find(delimiter));
			if (checkIsDouble(firstToken)) {
				size_t pos = line.find(delimiter);
				double mass = std::stod(line.substr(0, pos));
				line.erase(0, pos + delimiter.length());
				pos = line.find(delimiter);
				double x = std::stod(line.substr(0, pos));
				line.erase(0, pos + delimiter.length());
				pos = line.find(delimiter);
				double y = std::stod(line.substr(0, pos));
				line.erase(0, pos + delimiter.length());
				pos = line.find(delimiter);
				double z = std::stod(line.substr(0, pos));
				line.erase(0, pos + delimiter.length());
				pos = line.find(delimiter);
				double vx = std::stod(line.substr(0, pos));
				line.erase(0, pos + delimiter.length());
				pos = line.find(delimiter);
				double vy = std::stod(line.substr(0, pos));
				line.erase(0, pos + delimiter.length());
				double vz = std::stod(line, nullptr);

				stars.push_back(new Star(id,mass,x, y, z,vx,vy,vz));
				id++;
			}
		}
		sStars.close();
	}
	return stars;
}

void InOut::write(std::vector<double> x, std::vector<double> y, std::string filename, std::string header){
	if (x.size() != y.size()) {
		throw  "Vector size must be equal";
	}
	setOutputDirectory(filename);
	std::ofstream file(filename);
	if (header.size() > 0)
		file << header << '\n';
	//no NOT parallel this one
	for (int i = 0; i < x.size(); ++i) {
		file << x[i] <<','<< y[i] << '\n';
	}
	file.close();
}

void InOut::write(std::vector<Vec3D> line, std::string filename) {
	setOutputDirectory(filename);
	std::ofstream file(filename);
	for (int i = 0; i < line.size(); ++i) {
		file << line[i].print() << '\n';
	}
	file.close();
}

void InOut::write(std::vector<Point> points, std::string filename) {
	setOutputDirectory(filename);
	std::ofstream file(filename);
	for (Point point : points) {
		file << point << '\n';
	}
	file.close();
}

//@Note: replaced by LookupTable
//std::vector<std::vector<double>> InOut::readDoubleMatrix(std::string filname){
//	std::string line;
//	std::ifstream file(filname);
//	std::vector<double> row;
//	std::vector<std::vector<double>> matrix;
//	std::string delimiter = ",";
//
//	while (std::getline(file, line)){
//		std::string firstToken = line.substr(0, line.find(delimiter));
//		if (checkIsDouble(firstToken)) {
//			size_t pos = 0;
//			while ((pos = line.find(delimiter)) != std::string::npos) {
//				row.push_back(std::stod(line.substr(0, pos), nullptr));
//				line.erase(0, pos + delimiter.length());
//			}
//			row.push_back(std::stod(line, nullptr));
//			matrix.push_back(row);
//			row.clear();
//		}
//	}
//
//	return matrix;
//}