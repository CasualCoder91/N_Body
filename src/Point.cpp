#include "Point.h"

double Point::getDistance(const Point& pt2) {
    return sqrt((x - pt2.x) * (x - pt2.x) + (y - pt2.y) * (y - pt2.y));
}
double Point::getDelta(const Point& pt2) {
    return sqrt((vx - pt2.vx) * (vx - pt2.vx) + (vy - pt2.vy) * (vy - pt2.vy));
}

Point::Point(int id, double x, double y){
    this->id = id;
    this->x = x;
    this->y = y;
}
