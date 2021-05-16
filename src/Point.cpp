#include "Point.h"

const std::string Point::header = "id,cluster,clusterStar,x,y,vx,vy";

double Point::getDistance(const Point& pt2) {
    return sqrt((x - pt2.x) * (x - pt2.x) + (y - pt2.y) * (y - pt2.y));
}
double Point::getVelDelta(const Point& pt2) {
    return sqrt((velocity[0] - pt2.velocity[0]) * (velocity[0] - pt2.velocity[0]) + (velocity[1] - pt2.velocity[1]) * (velocity[1] - pt2.velocity[1]));
}

Point::Point(){
    this->cluster = 0;
}

Point::Point(int id, double x, double y, bool clusterStar, double magnitude){
    this->id = id;
    this->x = x;
    this->y = y;
    this->clusterStar = clusterStar;
    this->magnitude = magnitude;
    this->cluster = -1; //-1 = NOT_CLASSIFIED (Position & Velocity) | 0 = NOT_CLASSIFIED (Velocity)
    this->nNeighbors = 0;
    this->velocity[0] = 0;
    this->velocity[1] = 0;
}

Point::Point(int id, double x, double y, double vx, double vy, bool clusterStar, double magnitude)
{
    this->id = id;
    this->x = x;
    this->y = y;
    this->clusterStar = clusterStar;
    this->magnitude = magnitude;
    this->cluster = -1; //-1 = NOT_CLASSIFIED (Position & Velocity) | 0 = NOT_CLASSIFIED (Velocity)
    this->nNeighbors = 0;
    this->velocity[0] = vx;
    this->velocity[1] = vy;
}
