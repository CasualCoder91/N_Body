#pragma once
#include <cmath>

class Point {
public:
    double x, y, vx, vy;
    int id, nNeighbors, cluster;
    bool clusterStar;

    double getDistance(const Point& pt2);
    double getDelta(const Point& pt2);

    Point();
    Point(int id, double x, double y, bool clusterStar);
};