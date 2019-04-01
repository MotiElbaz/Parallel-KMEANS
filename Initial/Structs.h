#pragma once

#include "mpi.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Structs.h"

#define MASTER 0
#define OUTPUT "output.txt"
#define INPUT "input.txt"
#define TRUE 1
#define FALSE 0

struct Cluster
{
	int id, clusterSize;
	double x, y, z;
	double diameter;
}typedef Cluster;

struct Point
{
	int clusterId;
	double x, y, z;
	double vX, vY, vZ;
	double minDis;
}typedef Point;

MPI_Datatype createPointType();
MPI_Datatype createClusterType();
double KMEANS(Point* allPoints, int numOfPoints, Cluster* allClusters, int numOfClusters, double dT, int LIMIT, double QM, int T, double *interval);
Cluster* mergeClusters(Cluster* allClusters, int numOfClusters);
Point* calculatePostion(Point* allPoints, int numOfPoints, double interval);
double calculateQualityMeasure(Cluster* allClusters, int numOfClusters);
Cluster* calculateDiameter(Cluster* allClusters, Point* allPoints, int numOfPoints);
void recalculateClusterCenters(Cluster *allClusters, int numOfClusters, Point *allPoints, int numOfPoints);
int groupPoints(Cluster *allClusters, int numOfClusters, Point *allPoints, int numOfPoints);
Cluster* buildClusters(Point* allPoints, int numOfPoints, int numOfClusters);
Point* readFile(int* numOfPoints, int* numOfClusters, int* T, double* dT, int* LIMIT, double* QM);
void writeToFile(Cluster* allClusters, int size, double interval, double quailty);
void printClusters(Cluster* allClusters, int size);
void printPoints(Point* allPoints, int size);