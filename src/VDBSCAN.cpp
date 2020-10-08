#include "VDBSCAN.h"

VDBSCAN::VDBSCAN(int n, double epsSpace, int minPts, std::vector<Point> points) {
    this->n = n;
    this->epsSpace = epsSpace;
    this->minPts = minPts;
    this->points = points;
    this->size = (int)points.size();
    adjPoints.resize(size);
    this->clusterIdx = -1;
}

void VDBSCAN::run() {
    checkNearPoints();

    for (int i = 0; i < size; i++) {
        if (points[i].cluster != NOT_CLASSIFIED) continue;

        if (isCoreObject(i)) {
            depthFirstSearch(i, ++clusterIdx);
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

void VDBSCAN::depthFirstSearch(int now, int c) {
    points[now].cluster = c;
    if (!isCoreObject(now)) return;

    for (auto& next : adjPoints[now]) {
        if (points[next].cluster != NOT_CLASSIFIED) continue;
        depthFirstSearch(next, c);
    }
}

void VDBSCAN::checkNearPoints() {
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

bool VDBSCAN::isCoreObject(int idx) {
    return points[idx].nNeighbors >= minPts;
}

std::vector<std::vector<int> > VDBSCAN::getCluster() {
    return cluster;
}