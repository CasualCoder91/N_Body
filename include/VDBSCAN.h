#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
//#include <map>

#include "Point.h"

class VDBSCAN {
public:

    static const int NOISE = -2;
    static const int NOT_CLASSIFIED = -1;

    int minPts;
    double epsSpace;
    double epsTime = 2;
    //std::vector<Point> points;
    int size;
    std::vector<std::vector<int> > adjPoints;
    std::vector<bool> visited;
    std::vector<std::vector<int> > cluster;
    int clusterIdx;

    VDBSCAN(double epsSpace, double epsTime, int minPts);

    void run(std::vector<Point>& points);

    void depthFirstSearch(int now, int c, std::vector<Point>& points);

    void checkNearPoints(std::vector<Point>& points);
    
    bool isCoreObject(int pointID, std::vector<Point>& points);

    std::vector<std::vector<int> > getCluster();

};

//class InputReader {
//private:
//    ifstream fin;
//    vector<Point> points;
//public:
//    InputReader(string filename) {
//        fin.open(filename);
//        if (!fin) {
//            cout << filename << " file could not be opened\n";
//            exit(0);
//        }
//        parse();
//    }
//    void parse() {
//        int pointID;
//        double x, y;
//        while (!fin.eof()) {
//            fin >> pointID >> x >> y;
//            points.push_back({ x,y,0, NOT_CLASSIFIED });
//        }
//        points.pop_back();
//    }
//    vector<Point> getPoints() {
//        return points;
//    }
//};
//
//class OutputPrinter {
//private:
//    ofstream fout;
//    vector<vector<int> > cluster;
//    string filename;
//    int n;
//public:
//    OutputPrinter(int n, string filename, vector<vector<int> > cluster) {
//        this->n = n;
//        this->cluster = cluster;
//
//        // remove ".txt" from filename
//        if (filename.size() < 4) {
//            cout << filename << "input file name's format is wrong\n";
//            exit(0);
//        }
//        for (int i = 0; i < 4; i++) filename.pop_back();
//        this->filename = filename;
//
//        // sort by size decending order
//        sort(cluster.begin(), cluster.end(), [&](const vector<int> i, const vector<int> j) {
//            return (int)i.size() > (int)j.size();
//            });
//    }
//    void print() {
//        for (int i = 0; i < n; i++) {
//            fout.open(filename + "_cluster_" + to_string(i) + ".txt");
//
//            for (int j = 0; j < cluster[i].size(); j++) {
//                fout << cluster[i][j] << endl;
//            }
//
//            fout.close();
//        }
//    }
//};

//int main(int argc, const char* argv[]) {
//    if (argc != 5) {
//        std::cout << "Please follow this format. clustering.exe [intput] [n] [eps] [minPts]";
//        return 0;
//    }
//
//    std::string inputFileName(argv[1]);
//    std::string n(argv[2]);
//    std::string epsSpace(argv[3]);
//    std::string minPts(argv[4]);
//
//    InputReader inputReader(inputFileName);
//
//    DBCAN dbScan(stoi(n), stod(epsSpace), stoi(minPts), inputReader.getPoints());
//    dbScan.run();
//
//    OutputPrinter outputPrinter(stoi(n), inputFileName, dbScan.getCluster());
//    outputPrinter.print();
//
//    return 0;
//}
