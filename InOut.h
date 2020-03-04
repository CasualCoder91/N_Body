#pragma once
#include <fstream>
#include <vector>
#include <string>

#include "Star.h"
#include "Node.h"

namespace InOut{
	void write(std::vector<Star*> stars,std::string filename);
	void writeWithLabel(std::vector<Star*> stars, std::string filename);
	void writeAll(std::vector<Star*> stars, std::string filename);
	void write(Node* tree);
	void writeRecursively(std::ofstream* file_ptr, Node* node_ptr);
	void write(std::vector<double> x, std::vector<double> y, std::string filename);
};

