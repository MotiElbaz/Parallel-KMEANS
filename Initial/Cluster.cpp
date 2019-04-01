#include "Cluster.h"
#include "Functions.h"

Cluster* buildClusters(Point* allPoints, int numOfPoints, int numOfClusters,Cluster* allClusters)
{
	int i;
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < numOfClusters; i++)
	{
		allClusters[i].x = allPoints[i].x;
		allClusters[i].y = allPoints[i].y;
		allClusters[i].id = i;
		allClusters[i].clusterSize = 0;
		allClusters[i].diameter = 0;
	}
	return allClusters;
}

Cluster* calculateDiameter(Cluster* allClusters, Point* allPoints, int numOfPoints)
{
#pragma omp parallel
	for (int j = 0; j < numOfPoints; j++)
	{
#pragma omp for schedule(dynamic,1)
		for (int q = j + 1; q < numOfPoints; q++)
		{
			if (allPoints[j].clusterId == allPoints[q].clusterId)
			{
				double distance = sqrt(pow(allPoints[q].x - allPoints[j].x, 2) + pow(allPoints[q].y - allPoints[j].y, 2));
				if (allClusters[allPoints[j].clusterId].diameter < distance)
				{
					allClusters[allPoints[j].clusterId].diameter = distance;
				}
			}
		}
	}
}

void writeToFile(Cluster* allClusters, int size, double moment, double quailty)
{
	FILE* f;
	int i;
	f = fopen(OUTPUT, "w");
	fprintf(f, "First occurrence at t = %lf with q = %lf\n", &moment, &quailty);
	fprintf(f, "Centers of the clusters:\n");
	for (i = 0; i < size; i++)
	{
		fprintf(f, "%lf %lf\n", &allClusters[i].x, &allClusters[i].y);
	}
	fclose(f);
}

double calculateQualityMeasure(Cluster* allClusters, int numOfClusters)
{
	double q = 0.0, countOfVar = 0;
	for (int j = 0; j < numOfClusters; j++)
	{
		double temp = allClusters[j].diameter;
		for (int q = j + 1; q < numOfClusters; q++)
		{
			double distance = sqrt(pow(allClusters[q].x - allClusters[j].x, 2) + pow(allClusters[q].y - allClusters[j].y, 2));
			q += temp / distance;
			countOfVar++;
		}
	}
	q /= countOfVar;
	return q;
}