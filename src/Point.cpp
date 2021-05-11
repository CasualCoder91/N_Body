#include "Point.h"

const std::string Point::header = "id,cluster,clusterStar,x,y,vx,vy";

double Point::getDistance(const Point& pt2) {
    return sqrt((x - pt2.x) * (x - pt2.x) + (y - pt2.y) * (y - pt2.y));
}
double Point::getVelDelta(const Point& pt2) {
    return sqrt((vx - pt2.vx) * (vx - pt2.vx) + (vy - pt2.vy) * (vy - pt2.vy));
}

Point::Point(){}

Point::Point(int id, double x, double y, bool clusterStar, double magnitude){
    this->id = id;
    this->x = x;
    this->y = y;
    this->clusterStar = clusterStar;
    this->magnitude = magnitude;
    this->cluster = -1; //NOT_CLASSIFIED
    this->nNeighbors = 0;
    this->vx = 0;
    this->vy = 0;
}
