#pragma once

#define OUTPUT "output.txt"

struct Cluster
{
	int id;
	double x;
	double y;
	int clusterSize;
	double diameter;
}typedef Cluster;

Cluster* buildClusters(Point* allPoints, int numOfPoints, int numOfClusters,Cluster* allClusters);
Cluster* calculateDiameter(Cluster* allClusters, Point* allPoints, int numOfPoints);
void writeToFile(Cluster* allClusters, int size, double moment, double quailty);
double calculateQualityMeasure(Cluster* allClusters, int numOfClusters);