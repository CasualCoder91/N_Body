#pragma once
#include <iostream>
#include <cstdlib> // For rand()
#include <fstream>
#include <string>
#include <Windows.h>
#include <vector>
#include <ctime> // For using srand(time(0))

#include "Point.h"

class ClosestPair {
	//std::vector<Point> p;

	/* Calculate Distance between two points */
	static __forceinline double calDist(Point p1, Point p2);
	static std::vector<Point> merge(const std::vector<Point>& left, const  std::vector<Point>& right, std::string coord);
	static std::vector<Point> mergeSort(const std::vector<Point>& vect, std::string coord);
	static std::vector<Point> bruteForce(const std::vector<Point>& p);
	static std::vector<Point> closestPair(std::vector<Point>& points);
public:
	ClosestPair();
	~ClosestPair() = default;
	// return a pair with minimum distance
	static std::vector<Point> run(std::vector<Point>& points);
};




/*int main()
{

	string line;
	ifstream pointFile("input.txt");
	if (pointFile.is_open())
	{
		getline(pointFile, line);
		pNum = stoi(line);
		int i = 0;

		while (getline(pointFile, line))
		{
			double x = stod(line.substr(0, line.find(" ")));
			double y = stod(line.substr(line.find(" ")));
			Point pTemp;
			pTemp.x = x;
			pTemp.y = y;
			p.push_back(pTemp);
			i++;
		}
		pointFile.close();
	}
	else
		cout << "There was a problem opening or reading pair file" << endl;


	// Devide & Conquer Strategy
	clock_t beginDevide = clock();
	cout << "Devide & Conquer:" << endl;
	if (p.size() < 2) { // In case we have only two points
						// return two points and their distance
		cout << "p1: " << p[0].x << " " << p[0].y << " , p2: " << p[1].x << " " << p[1].y << calDist(p[0], p[1]) << endl;
	}
	else {
		vector<Point> sortedPx;
		sortedPx = mergeSort(p, "x");
		vector<Point> closestP = closestPair(sortedPx);
		clock_t endDevide = clock();
		double elapsedSecsDevide = double(endDevide - beginDevide) / CLOCKS_PER_SEC;
		cout << "Min Dist: " << calDist(closestP[0], closestP[1]) << endl;
		cout << "Closest Pair: (" << closestP[0].x << "," << closestP[0].y << ") (" << closestP[1].x << "," << closestP[1].y << ")" << endl;
		cout << "Total Running Time of Devide & Conquer: " << elapsedSecsDevide << endl;
		cout << "Devide & Conquer Counter: " << counterDevide << endl;
	}
	// End of Devide & Conquer Strategy

	char ch;
	cin >> ch;
	return 0;
}*/