#include "ClosestPair.h"

__forceinline double ClosestPair::calDist(Point p1, Point p2) {
	return sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2));
}

std::vector<Point> ClosestPair::merge(const std::vector<Point>& left, const  std::vector<Point>& right, std::string coord) {
	int i = 0, j = 0, k = 0; //i:return vect j:left k:right
	int lsize = left.size();
	int rsize = right.size();
	int size = lsize + rsize;
	std::vector<Point> sortedVect(size);
	for (i = 0; i < size; i++) {
		if (j < lsize && k < rsize) {
			// if coord="x" or coord="y"
			if (coord == "x") { // Sort By X-coordinate
				if (left[j].x <= right[k].x) {
					sortedVect[i] = left[j];
					j++;
				}
				else {
					sortedVect[i] = right[k];
					k++;
				}
			}
			else { // Sort By Y-coordinate
				if (left[j].y <= right[k].y) {
					sortedVect[i] = left[j];
					j++;
				}
				else {
					sortedVect[i] = right[k];
					k++;
				}
			}

		}
		else if (j < lsize) {
			sortedVect[i] = left[j];
			j++;
		}
		else if (k < rsize) {
			sortedVect[i] = right[k];
			k++;
		}
	}
	return sortedVect;
}

std::vector<Point> ClosestPair::mergeSort(const std::vector<Point>& vect, std::string coord) {
	int size = vect.size();
	int mid = size / 2;
	std::vector<Point> left;
	std::vector<Point> right;
	if (size <= 1)//base case
		return vect;
	for (int i = 0; i < mid; i++)
	{
		left.push_back(vect[i]);
		right.push_back(vect[i + mid]);
	}
	if (size % 2 != 0)
		right.push_back(vect.back());//Push last element of vector to the right vector in case size of vect is odd

	left = mergeSort(left, coord);
	right = mergeSort(right, coord);

	return merge(left, right, coord);
}

std::vector<Point> ClosestPair::bruteForce(const std::vector<Point>& p) {
	int psize = p.size();
	if (psize < 2)
		return p;
	else
	{
		double min = calDist(p[0], p[1]);
		Point p1 = p[0], p2 = p[1];
		for (int i = 0; i < psize; i++)
			for (int j = i + 1; j < psize; j++)
			{
				double distpp = calDist(p[i], p[j]);
				if (distpp < min) {
					min = distpp;
					//Two min distance pair:
					p1 = p[i];
					p2 = p[j];
				}
			}
		std::vector<Point> minPoints;
		minPoints.push_back(p1);
		minPoints.push_back(p2);
		return minPoints;
	}
}

ClosestPair::ClosestPair(){}

std::vector<Point> ClosestPair::run(std::vector<Point>& points){
	std::vector<Point> sortedPx = mergeSort(points, "x");
	return closestPair(sortedPx);
}

std::vector<Point> ClosestPair::closestPair(std::vector<Point>& points) {

	if (points.size() <= 3) {
		// return min distance between these using brute-force
		return bruteForce(points);
	}
	int pointsSize = points.size();
	int half = pointsSize / 2;
	std::vector<Point> leftPairs;
	std::vector<Point> rightPairs;
	std::vector<Point> leftResult, rightResult;

	if (pointsSize % 2 == 0) {
		for (int i = 0; i < half; i++)
		{
			leftPairs.push_back(points[i]);
			rightPairs.push_back(points[i + half]);
		}
	}
	else
	{
		for (int i = 0; i <= half; i++)
		{
			leftPairs.push_back(points[i]);
		}
		for (int i = half + 1; i < pointsSize; i++)
		{
			rightPairs.push_back(points[i]);
		}
	}

	leftResult = closestPair(leftPairs);
	rightResult = closestPair(rightPairs);
	double minLeft = calDist(leftResult[0], leftResult[1]);
	double minRight = calDist(rightResult[0], rightResult[1]);
	double minLR = min(minLeft, minRight);
	std::vector<Point> stripPoints;
	// Since we have a sorted points based on their x,
	// It suffices to check only points around px
	for (int i = half - 1; i > 0; i--)
	{
		if (points[half].x - points[i].x < minLR)
			stripPoints.push_back(points[i]);
		else
			break; // In case px>min
	}
	for (int i = half + 1; i < pointsSize; i++)
	{
		if (points[i].x - points[half].x < minLR)
			stripPoints.push_back(points[i]);
		else
			break; // In case px>min
	}
	stripPoints.push_back(points[half]);
	std::vector<Point> sortedStripPointsy;
	std::vector<Point> stripClosestPair;
	if (stripPoints.size() > 1) // At least one 2 points be in the strip area
	{
		// If there is a point in Strip Area
		// sort stripPoint by their y-coordinates
		double min = minLR;
		sortedStripPointsy = mergeSort(stripPoints, "y");
		int sizePointStrip = sortedStripPointsy.size();
		for (int i = 0; i < sizePointStrip - 1; i++)
			for (int j = i + 1; j < sizePointStrip && (sortedStripPointsy[j].y - sortedStripPointsy[i].y < min); j++)
			{
				double pairDistance = calDist(sortedStripPointsy[i], sortedStripPointsy[j]);
				if (pairDistance < min)
				{
					min = pairDistance;
					stripClosestPair.clear();
					stripClosestPair.push_back(sortedStripPointsy[i]);
					stripClosestPair.push_back(sortedStripPointsy[j]);
				}
			}

		if (stripClosestPair.size() > 0) // If closest pair is in the strip line
			return stripClosestPair;
		else
			return (minLeft < minRight) ? leftResult : rightResult;
	}
	else
		return (minLeft < minRight) ? leftResult : rightResult;
}