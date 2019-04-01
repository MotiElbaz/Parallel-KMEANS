#pragma once

#define INPUT "input.txt"
#define TRUE 1
#define FALSE 0

struct Point
{
	int clusterId;
	double x;
	double y;
	double vX;
	double vY;
	double minDis;
}typedef Point;

Point* resetMinDis(Point* allPoints, int size);
int groupPoints(int numOfClusters, Cluster** allClusters, Point** allPoints, int numOfPoints);
Point* calculatePostion(Point* allPoints, int size, double moment);
Point* readFile(int* numOfPoints, int* numOfClusters, int* maxIters, int* interval, double* qualityMeasure, double* deltaInterval);
void printPoints(Point* allPoints, int size);