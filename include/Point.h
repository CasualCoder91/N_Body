#pragma once
#include <cmath>

class Point {
public:
    double x, y, vx, vy;
    int id, nNeighbors, cluster;
    double getDistance(const Point& pt2);
    double getDelta(const Point& pt2);

    Point(int id, double x, double y);
};