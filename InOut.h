/**
 * Collection of functions handling reading from and writing to files.
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#pragma once
#include <fstream>
#include <vector>
#include <string>

#include "Star.h"
#include "Node.h"

namespace InOut{
    /**
     @brief Writes all star coordinates into a file. Format according to Vec3D::print()
     @param stars Vector of star pointers. All elements are written to file
     @param filename The name and location of the file. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     */
	void write(std::vector<Star*> stars,std::string filename);
    /**
     @brief Writes all star coordinates and star id into a file. Format according to Vec3D::print(). Id is the memory location of the respective star.
     @param stars Vector of star pointers. All elements are written to file
     @param filename The name and location of the file. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     */
	void writeWithLabel(std::vector<Star*> stars, std::string filename);
    /**
     @brief Writes all star member variables into a file using Star::dump()
     @param stars Vector of star pointers. All elements are written to file.
     @param filename The name and location of the file. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     @note meant for debugging purposes.
     */
	void writeAll(std::vector<Star*> stars, std::string filename);
    /**
     @brief Writes octree cell coorindates (top left front and bottom right back) into a file. Format: Vec3D::print(),Vec3D::print().
     @param tree Pointer to the root of the tree.
     @param filename The name and location of the file. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     @note **Only** root node pointer is acceptable parameter.
     */
	void write(Node* tree, std::string filename);
    /**
     @brief Writes std::vector elements to a file. Format: first x element, first y emelent
     @param x,y std::vector meant to be coordinates. vector size must be equal.
     @param filename The name and location of the file. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     */
	void write(std::vector<double> x, std::vector<double> y, std::string filename);
    /**
     @brief Function used by write(Node* tree, std::string filename).
     @attention do **not** call this function
     */
    void writeRecursively(std::ofstream* file_ptr, Node* node_ptr);
};

