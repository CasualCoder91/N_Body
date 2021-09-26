#pragma once
#include <cmath>
#include <string>
#include <iostream>

class Point {
public:
    double x, y;
    double velocity[2] = { 0,0 };
    int id, nNeighbors, cluster, fk_star;
    bool clusterStar;
    double magnitude;
    static const std::string header;

    Point();
    Point(int id, double x, double y, bool clusterStar, double magnitude);
    Point(int id, double x, double y, double vx, double vy,  bool clusterStar, double magnitude);

    double getDistance(const Point& pt2);
    double getVelDelta(const Point& pt2);
    double distance_origin();

    friend std::ostream& operator<<(std::ostream& o, Point const& point) {
        o << point.id << ',' << point.cluster << ',' << point.clusterStar << ',' << point.x << ',' << point.y << ',' << point.velocity[0] << ',' << point.velocity[1];
        return o;
    }

};