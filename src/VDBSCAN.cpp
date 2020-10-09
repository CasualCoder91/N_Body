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
        if (points[i].cluster != NOT_CLASSIFIED) continue;

        if (isCoreObject(i, points)) {
            depthFirstSearch(i, ++clusterIdx, points);
        }
        else {
            points[i].cluster = NOISE;
        }
    }

    cluster.resize(clusterIdx + 1);
    for (int i = 0; i < size; i++) {
        if (points[i].cluster != NOISE) {
            cluster[points[i].cluster].push_back(i);
        }
    }
}

void VDBSCAN::depthFirstSearch(int now, int c, std::vector<Point>& points) {
    points[now].cluster = c;
    if (!isCoreObject(now, points)) return;

    for (auto& next : adjPoints[now]) {
        if (points[next].cluster != NOT_CLASSIFIED) continue;
        depthFirstSearch(next, c, points);
    }
}

void VDBSCAN::checkNearPoints(std::vector<Point>& points) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j) continue;
            if (points[i].getDistance(points[j]) <= epsSpace && points[i].getDelta(points[j]) <= epsTime) {
                points[i].nNeighbors++;
                adjPoints[i].push_back(j);
            }
        }
    }
}

bool VDBSCAN::isCoreObject(int idx, std::vector<Point>& points) {
    return points[idx].nNeighbors >= minPts;
}

std::vector<std::vector<int> > VDBSCAN::getCluster() {
    return cluster;
}