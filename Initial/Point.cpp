#include "Functions.h"
#include "Point.h"

Point* resetMinDis(Point* allPoints, int size)
{
	int i;
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < size; i++)
	{
		allPoints[i].minDis = -1;
	}
	return allPoints;
}
int groupPoints(int numOfClusters, Cluster** allClusters, Point** allPoints, int numOfPoints)
{
	int ifNoPointMoved = TRUE;
	for (int j = 0; j < numOfClusters; j++)
	{
		double avgX = 0, avgY = 0;
		for (int q = 0; q < numOfPoints; q++)
		{
			int originalId = allPoints[j]->clusterId, nextId;
			double distance = sqrt(pow(allPoints[q]->x - allClusters[j]->x, 2) + pow(allPoints[q]->y - allClusters[j]->y, 2));
			if (allPoints[q]->minDis == -1)
			{
				allPoints[q]->minDis == distance;
				allPoints[q]->clusterId = allClusters[j]->id;
				allClusters[j]->clusterSize++;
				avgX += allPoints[q]->x;
				avgY += allPoints[q]->y;
			}
			else if (allPoints[q]->minDis > distance)
			{
				allPoints[q]->minDis == distance;
				allPoints[q]->clusterId = allClusters[j]->id;
				allClusters[j]->clusterSize++;
				avgX += allPoints[q]->x;
				avgY += allPoints[q]->y;
			}
			nextId = allPoints[j]->clusterId;
			if (originalId != nextId)
			{
				ifNoPointMoved = FALSE;
			}
		}
		allClusters[j]->x = avgX / allClusters[j]->clusterSize;
		allClusters[j]->y = avgY / allClusters[j]->clusterSize;
	}
	return ifNoPointMoved;
}
Point* calculatePostion(Point* allPoints, int size, double moment)
{
	int i;
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < size; i++)
	{
		allPoints[i].x = allPoints[i].x + moment*allPoints[i].vX;
		allPoints[i].y = allPoints[i].y + moment*allPoints[i].vY;
	}
	return allPoints;
}
Point* readFile(int* numOfPoints, int* numOfClusters, int* maxIters, int* interval, double* qualityMeasure, double* deltaInterval)
{
	FILE* f;
	int i;
	f = fopen(INPUT, "r");
	fscanf(f, "%d %d %d %lf %d %lf\n", numOfPoints, numOfClusters, maxIters, qualityMeasure, interval, deltaInterval);
	Point* allPoints = (Point*)malloc(*numOfPoints * sizeof(Point));
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < *numOfPoints; i++)
	{
		fscanf(f, "%lf %lf %lf %lf\n", &allPoints[i].x, &allPoints[i].y, &allPoints[i].vX, &allPoints[i].vY);
		allPoints[i].minDis = -1;
		allPoints[i].clusterId = -1;
	}
	fclose(f);
	return allPoints;
}
void printPoints(Point* allPoints, int size)
{
	int i;
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < size; i++)
	{
		printf("%lf %lf %lf %lf \n", allPoints[i].x, allPoints[i].y, allPoints[i].vX, allPoints[i].vY);
	}
}