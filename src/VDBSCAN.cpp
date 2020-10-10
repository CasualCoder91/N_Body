#include "VDBSCAN.h"

VDBSCAN::VDBSCAN(double epsSpace, double epsTime, int minPts) {
    this->epsSpace = epsSpace;
    this->epsTime = epsTime;
    this->minPts = minPts;
    this->clusterIdx = -1;
}

void VDBSCAN::run(std::vector<Point>& points) {
    this->size = (int)points.size();
    adjPoints.resize(size);
    checkNearPoints(points);

    for (int i = 0; i < size; i++) {
        if (points.at(i).cluster != NOT_CLASSIFIED) continue;

        if (isCoreObject(i, points)) {
            depthFirstSearch(i, ++clusterIdx, points);
        }
        else {
            points.at(i).cluster = NOISE;
        }
    }

    cluster.resize(clusterIdx + 1);
    for (int i = 0; i < size; i++) {
        if (points.at(i).cluster != NOISE) {
            cluster.at(points.at(i).cluster).push_back(i);
        }
    }
}

void VDBSCAN::depthFirstSearch(int now, int c, std::vector<Point>& points) {
    points.at(now).cluster = c;
    if (!isCoreObject(now, points)) return;

    for (int next : adjPoints.at(now)) {
        if (points.at(next).cluster != NOT_CLASSIFIED) continue;
        depthFirstSearch(next, c, points);
    }
}

void VDBSCAN::checkNearPoints(std::vector<Point>& points) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j) continue;
            if (points.at(i).getDistance(points.at(j)) <= epsSpace && points.at(i).getDelta(points.at(j)) <= epsTime) {
                points.at(i).nNeighbors++;
                adjPoints.at(i).push_back(j);
            }
        }
    }
}

bool VDBSCAN::isCoreObject(int idx, std::vector<Point>& points) {
    return points.at(idx).nNeighbors >= minPts;
}

std::vector<std::vector<int> > VDBSCAN::getCluster() {
    return cluster;
}