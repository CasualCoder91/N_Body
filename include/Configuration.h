//https://wiki.calculquebec.ca/w/C%2B%2B_:_fichier_de_configuration/en

#pragma once

#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Vec3D.h"

class Configuration
{
public:

    Configuration(const std::string& file);
    // clear all values
    void Clear();

    // load a configuration filePath
    bool Load(const std::string& File);

    // check if value associated with given key exists
    bool Contains(const std::string& key) const;

    // get value associated with given key
    bool Get(const std::string& key, std::string& value) const;
    bool Get(const std::string& key, int& value) const;
    bool Get(const std::string& key, long& value) const;
    bool Get(const std::string& key, double& value) const;
    bool Get(const std::string& key, bool& value) const;
    bool Get(const std::string& key, Vec3D& value) const;
    double GetDouble(const std::string& key) const;
    int GetInt(const std::string& key) const;
    Vec3D GetVec3D(const std::string& key) const;
    std::string GetString(const std::string& key) const;
    bool GetBool(const std::string& key) const;
    std::vector<double> GetDoubleVector(const std::string& key) const;

private:
    // the container
    std::map<std::string, std::string> data;

    // remove leading and trailing tabs and spaces
    static std::string Trim(const std::string& str);
};

