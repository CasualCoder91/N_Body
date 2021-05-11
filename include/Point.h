#pragma once
#include <cmath>
#include <string>
#include <iostream>

class Point {
public:
    double x, y, vx, vy;
    int id, nNeighbors, cluster;
    bool clusterStar;
    double magnitude;
    static const std::string header;

    double getDistance(const Point& pt2);
    double getVelDelta(const Point& pt2);

    Point();
    Point(int id, double x, double y, bool clusterStar, double magnitude);

    friend std::ostream& operator<<(std::ostream& o, Point const& point) {
        o << point.id << ',' << point.cluster << ',' << point.clusterStar << ',' << point.x << ',' << point.y << ',' << point.vx << ',' << point.vy;
        return o;
    }

};