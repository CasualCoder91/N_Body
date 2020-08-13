/**
 * Collection of functions handling reading from and writing to files.
 * Writing into the path of the executable is not possible. If the passed path contains only the filename the file will be written to /Output.
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#pragma once
#include <fstream>
#include <vector>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <direct.h>

#include "Star.h"
#include "Node.h"

namespace InOut{

    /**@brief By default all files are writen within this subdirectory */
    static std::string outputDirectory = "Output/";

    /**@brief Creates given directories if they do not exist */
    std::string makeDirectory(std::string path);
    /**
     @brief if parameter is just a filename (contains no /) it is modified to poit to the default Directory 
     @param [in,out] filePath if it just contains the filename the path is set to the default directory.
    */
    static void setOutputDirectory(std::string& filePath);
    /**
     @brief Writes all star coordinates into a filePath. Format according to Vec3D::print()
     @param stars Vector of star pointers. All elements are written to filePath
     @param filename The name and location of the filePath. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     */
	void write(std::vector<Star*> stars,std::string filename);
    /**
     @brief Writes all star coordinates and star id into a filePath. Format according to Vec3D::print(). Id is the memory location of the respective star.
     @param stars Vector of star pointers. All elements are written to filePath
     @param filename The name and location of the filePath. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     */
	void writeWithLabel(std::vector<Star*> stars, std::string filename);
    /**
     @brief Writes all star member variables into a filePath using Star::dump()
     @param stars Vector of star pointers. All elements are written to filePath.
     @param filename The name and location of the filePath. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     @note meant for debugging purposes.
     */
	void writeAll(std::vector<Star*> stars, std::string filename);
    /**
     @brief Writes octree cell coorindates (top left front and bottom right back) into a filePath. Format: Vec3D::print(),Vec3D::print().
     @param tree Pointer to the root of the tree.
     @param filename The name and location of the filePath. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     @note **Only** root node pointer is acceptable parameter.
     */
	void write(Node* tree, std::string filename);
    /**
     @brief Writes std::vector elements to a filePath. Format: first x element, first y emelent
     @param x,y std::vector meant to be coordinates. vector size must be equal.
     @param filename The name and location of the filePath. Can be used as path relative to excelcuateble ie "folder/filename.dat"
     @param header If given this will be written at the start of the generated file.
     */
	void write(std::vector<double> x, std::vector<double> y, std::string filename, std::string header = "");
    /**
     @brief Function used by write(Node* tree, std::string filename).
     @attention do **not** call this function
     */
    void writeRecursively(std::ofstream* file_ptr, Node* node_ptr);
    /**
     @brief Each Vec3D entry is written into one separate line into the passed \p filename.
     */
    void write(std::vector<Vec3D> line, std::string filename);
    /**
     @brief reads Vec3D objects from file (format for each line: x,y,z).
     @todo: test it ;-)
     */
    std::vector<Vec3D> readVectors(std::string filename);
    /** @brief checks if the given \p inputString can be converted into a double. */
    bool checkIsDouble(std::string inputString);

    std::vector<Star*> readMcLuster(int firstID, std::string filename);

};

