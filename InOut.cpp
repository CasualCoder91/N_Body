#include "InOut.h"

void InOut::write(std::vector<Star*> stars, std::string filename) {
	std::ofstream file(filename);
	for (Star* star : stars) {
		file << star->position.x << ',' << star->position.y << ',' << star->position.z << '\n';
	}
	file.close();
}

void InOut::WriteWithLabel(std::vector<Star*> stars, std::string filename){
	std::ofstream file(filename);
	#pragma omp parallel for
	for (int i = 0; i < stars.size();++i) {
		file << stars.at(i)->position.x << ',' << stars.at(i)->position.y << ',' << stars.at(i)->position.z << ',' << i << '\n';
	}
	file.close();
}

void InOut::WriteAll(std::vector<Star*> stars, std::string filename){
	std::ofstream file(filename);
	#pragma omp parallel for
	for (int i = 0; i < stars.size(); ++i) {
		file << "Star: " << i << '\n';
		file << stars.at(i)->dump() << '\n';
		++i;
	}
	file.close();
}

void InOut::write(Node* tree){
	if (!tree->IsRoot()) {
		throw "Write function may only be called on root node";
	}
	else {
		std::ofstream file("tree.dat");
		WriteRecursively(&file, tree);
		file.close();
	}
	return;
}

void InOut::WriteRecursively(std::ofstream* file_ptr,Node* node_ptr) {
	for (Node* child : node_ptr->children) {
		if (child) {
			InOut::WriteRecursively(file_ptr, child);
		}
	}
	*file_ptr << node_ptr->Print();
	return;
}

void InOut::write(std::vector<double> x, std::vector<double> y, std::string filename){
	if (x.size() != y.size()) {
		throw  "Vector size must be equal";
	}
	std::ofstream file(filename);
	#pragma omp parallel for
	for (int i = 0; i < x.size(); ++i) {
		file << x.at(i) <<','<< y.at(i) << '\n';
	}
	file.close();
}
