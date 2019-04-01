#include <stdio.h>
#include "mpi.h"
#include <omp.h>
#include "Structs.h"

int main(int argc, char *argv[])
{
	double t1, t2, time, dT, QM, qualityMeasure, interval = 0.0;
	int T, numOfPoints, numOfClusters, LIMIT, numOfProcs, myRank, partSize, restSize;
	Point* allPoints;
	Cluster* allClusters;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
	MPI_Status status;
	MPI_Datatype PointMPIType = createPointType();
	MPI_Datatype ClusterMPIType = createClusterType();

	if (myRank == MASTER)
	{
		int i;
		t1 = MPI_Wtime();
		allPoints = readFile(&numOfPoints, &numOfClusters, &T, &dT, &LIMIT, &QM);
		allClusters = buildClusters(allPoints, numOfPoints, numOfClusters);
		partSize = numOfPoints / numOfProcs;
		restSize = numOfPoints % numOfProcs;
		for (i = 1; i < numOfProcs; i++)
		{
			MPI_Send(&partSize, 1, MPI_INT, i, MASTER, MPI_COMM_WORLD);
			MPI_Send(&numOfClusters, 1, MPI_INT, i, MASTER, MPI_COMM_WORLD);
			MPI_Send(allPoints + partSize*(i - 1), partSize, PointMPIType, i, MASTER, MPI_COMM_WORLD);
			MPI_Send(allClusters, numOfClusters, ClusterMPIType, i, MASTER, MPI_COMM_WORLD);
			MPI_Send(&LIMIT, 1, MPI_INT, i, MASTER, MPI_COMM_WORLD);
			MPI_Send(&T, 1, MPI_INT, i, MASTER, MPI_COMM_WORLD);
			MPI_Send(&dT, 1, MPI_DOUBLE, i, MASTER, MPI_COMM_WORLD);
			MPI_Send(&QM, 1, MPI_DOUBLE, i, MASTER, MPI_COMM_WORLD);
		}
		allPoints = allPoints + partSize*(i - 1);
		numOfPoints = partSize + restSize;
	}
	else {
		MPI_Recv(&numOfPoints, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&numOfClusters, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);

		allPoints = (Point*)malloc(numOfPoints * sizeof(Point));
		allClusters = (Cluster*)malloc(numOfClusters * sizeof(Cluster));

		MPI_Recv(allPoints, numOfPoints, PointMPIType, MASTER, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(allClusters, numOfClusters, ClusterMPIType, MASTER, 0, MPI_COMM_WORLD, &status);

		MPI_Recv(&LIMIT, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&T, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&dT, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&QM, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);
	}

	qualityMeasure = KMEANS(allPoints, numOfPoints, allClusters, numOfClusters, dT, LIMIT, QM, T, &interval);

	if (myRank == MASTER)
	{
		t2 = MPI_Wtime();
		time = t2 - t1;
		printf("Time : %lf\n", time);
		writeToFile(allClusters, numOfClusters, interval, qualityMeasure);
	}

	MPI_Finalize();
	return 0;
}

MPI_Datatype createPointType()
{
	struct Point point;
	int count = 8;
	MPI_Datatype PointMPIType;
	MPI_Datatype type[8] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE , MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE };
	int blocklen[8] = { 1, 1, 1, 1, 1, 1 , 1, 1 };
	MPI_Aint disp[8];

	// Create MPI user data type for partical
	disp[0] = (char *)&point.clusterId - (char *)&point;
	disp[1] = (char *)&point.x - (char *)&point;
	disp[2] = (char *)&point.y - (char *)&point;
	disp[3] = (char *)&point.z - (char *)&point;
	disp[4] = (char *)&point.vX - (char *)&point;
	disp[5] = (char *)&point.vY - (char *)&point;
	disp[6] = (char *)&point.vZ - (char *)&point;
	disp[7] = (char *)&point.minDis - (char *)&point;
	MPI_Type_create_struct(8, blocklen, disp, type, &PointMPIType);
	MPI_Type_commit(&PointMPIType);
	return PointMPIType;
}

MPI_Datatype createClusterType()
{
	struct Cluster cluster;
	int count = 6;
	MPI_Datatype ClusterMPIType;
	MPI_Datatype type[6] = { MPI_INT,MPI_INT, MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE };
	int blocklen[6] = { 1, 1, 1, 1, 1 , 1 };
	MPI_Aint disp[6];

	// Create MPI user data type for partical
	disp[0] = (char *)&cluster.id - (char *)&cluster;
	disp[1] = (char *)&cluster.clusterSize - (char *)&cluster;
	disp[2] = (char *)&cluster.x - (char *)&cluster;
	disp[3] = (char *)&cluster.y - (char *)&cluster;
	disp[4] = (char *)&cluster.z - (char *)&cluster;
	disp[5] = (char *)&cluster.diameter - (char *)&cluster;
	MPI_Type_create_struct(6, blocklen, disp, type, &ClusterMPIType);
	MPI_Type_commit(&ClusterMPIType);
	return ClusterMPIType;
}

double KMEANS(Point* allPoints, int numOfPoints, Cluster* allClusters, int numOfClusters, double dT, int LIMIT, double QM, int T, double *interval)
{
	double qualityMeasure = 0;
	int i = 0, n = 0, ifPointMoved = FALSE, globalPointMoved = FALSE, myRank, numOfProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
	MPI_Status status;
	MPI_Datatype ClusterMPIType = createClusterType();
	do {
		*interval = n*dT;
		allPoints = calculatePostion(allPoints, numOfPoints, *interval);
		do {
			globalPointMoved = FALSE;
			ifPointMoved = groupPoints(allClusters, numOfClusters, allPoints, numOfPoints);
			recalculateClusterCenters(allClusters, numOfClusters, allPoints, numOfPoints);
			MPI_Allreduce(&ifPointMoved, &globalPointMoved, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if (myRank == MASTER)
			{
				allClusters = mergeClusters(allClusters, numOfClusters);
			}
			else {
				MPI_Send(allClusters, numOfClusters, ClusterMPIType, MASTER, 0, MPI_COMM_WORLD);
				MPI_Recv(allClusters, numOfClusters, ClusterMPIType, MASTER, 0, MPI_COMM_WORLD, &status);
			}
			if (globalPointMoved == FALSE)
			{
				break;
			}
		} while (++i < LIMIT);
		allClusters = calculateDiameter(allClusters, allPoints, numOfPoints);		//calculate diameters with kuda
		if (myRank == MASTER) {
			qualityMeasure = calculateQualityMeasure(allClusters, numOfClusters);
			for (int j = 1; j < numOfProcs; j++) {
				MPI_Send(&qualityMeasure, 1, MPI_DOUBLE, j, MASTER, MPI_COMM_WORLD);
				MPI_Send(allClusters, numOfClusters, ClusterMPIType, j, MASTER, MPI_COMM_WORLD);
			}
		}
		else {
			MPI_Recv(&qualityMeasure, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(allClusters, numOfClusters, ClusterMPIType, MASTER, 0, MPI_COMM_WORLD, &status);
		}
		n++;
	} while (qualityMeasure > QM && *interval < T);
	return qualityMeasure;
}

Cluster* mergeClusters(Cluster* allClusters, int numOfClusters)
{
	int j = 0, myRank, numOfProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
	MPI_Status status;
	MPI_Datatype ClusterMPIType = createClusterType();
	Cluster* bufferClusters = (Cluster*)malloc(sizeof(Cluster)*numOfClusters);
	for (j = 1; j < numOfProcs; j++)
	{
		MPI_Recv(bufferClusters, numOfClusters, ClusterMPIType, j, 0, MPI_COMM_WORLD, &status);
		for (int k = 0; k < numOfClusters; k++)
		{
			allClusters[k].clusterSize += bufferClusters[k].clusterSize;
			allClusters[k].x += bufferClusters[k].x;
			allClusters[k].x /= 2;
			allClusters[k].y += bufferClusters[k].y;
			allClusters[k].y /= 2;
			allClusters[k].z += bufferClusters[k].z;
			allClusters[k].z /= 2;
			if (allClusters[k].diameter < bufferClusters[k].diameter)
			{
				allClusters[k].diameter = bufferClusters[k].diameter;
			}
		}
	}
	for (j = 1; j < numOfProcs; j++) {
		MPI_Send(allClusters, numOfClusters, ClusterMPIType, j, MASTER, MPI_COMM_WORLD);
	}
	return allClusters;
}

Point* calculatePostion(Point* allPoints, int numOfPoints, double interval)
{
	int i;
	for (i = 0; i < numOfPoints; i++)
	{
		allPoints[i].x = allPoints[i].x + interval*allPoints[i].vX;
		allPoints[i].y = allPoints[i].y + interval*allPoints[i].vY;
		allPoints[i].z = allPoints[i].z + interval*allPoints[i].vZ;
	}
	return allPoints;
}

double calculateQualityMeasure(Cluster* allClusters, int numOfClusters)
{
	int countOfVar = 0;
	double quality = 0.0;

	for (int j = 0; j < numOfClusters; j++)
	{
		double diameter = allClusters[j].diameter;
		for (int q = 0; q < numOfClusters; q++)
		{
			if (j != q)
			{
				double distance = sqrt(pow(allClusters[q].x - allClusters[j].x, 2) + pow(allClusters[q].y - allClusters[j].y, 2) + pow(allClusters[q].z - allClusters[j].z, 2));
				diameter /= distance;
				quality += diameter;
				countOfVar++;
			}
		}
	}
	quality /= countOfVar;
	return quality;
}

Cluster* calculateDiameter(Cluster* allClusters, Point* allPoints, int numOfPoints)
{
	int j, q;
#pragma omp parallel for
	for (j = 0; j < numOfPoints - 1; j++)
	{
		for (q = j + 1; q < numOfPoints; q++)
		{
			if (allPoints[j].clusterId == allPoints[q].clusterId)
			{
				double distance = sqrt(pow(allPoints[q].x - allPoints[j].x, 2) + pow(allPoints[q].y - allPoints[j].y, 2) + pow(allPoints[q].z - allPoints[j].z, 2));
				if (allClusters[allPoints[j].clusterId].diameter < distance)
				{
					allClusters[allPoints[j].clusterId].diameter = distance;
				}
			}
		}
	}
	return allClusters;
}

void recalculateClusterCenters(Cluster *allClusters, int numOfClusters, Point *allPoints, int numOfPoints)
{
	int i = 0;
	for (int j = 0; j < numOfPoints; j++)
	{
		i = allPoints[j].clusterId;
		allClusters[i].x += allPoints[j].x;
		allClusters[i].y += allPoints[j].y;
		allClusters[i].z += allPoints[j].z;
	}
	for (i = 0; i < numOfClusters; i++)
	{
		if (allClusters[i].clusterSize != 0)
		{
			allClusters[i].x /= allClusters[i].clusterSize;
			allClusters[i].y /= allClusters[i].clusterSize;
			allClusters[i].z /= allClusters[i].clusterSize;
		}
	}
}

int groupPoints(Cluster *allClusters, int numOfClusters, Point *allPoints, int numOfPoints)
{
	int ifPointMoved = FALSE, i, j;
	for (i = 0; i < numOfPoints; i++)
	{
		for (j = 0; j < numOfClusters; j++)
		{
			int originalId = allPoints[i].clusterId;
			double distance = sqrt(pow(allPoints[i].x - allClusters[j].x, 2) + pow(allPoints[i].y - allClusters[j].y, 2) + pow(allPoints[i].z - allClusters[j].z, 2));
			if (allPoints[i].minDis == -1)
			{
				allPoints[i].minDis = distance;
				allPoints[i].clusterId = allClusters[j].id;
				allClusters[j].clusterSize++;
				ifPointMoved = TRUE;
			}
			else if (allPoints[i].minDis > distance)
			{
				allPoints[i].minDis = distance;
				if (originalId != allClusters[j].id)
				{
					allClusters[originalId].clusterSize--;
					allClusters[j].clusterSize++;
					allPoints[i].clusterId = allClusters[j].id;
					ifPointMoved = TRUE;
				}
			}
		}
	}
	return ifPointMoved;
}

Cluster* buildClusters(Point* allPoints, int numOfPoints, int numOfClusters)
{
	int i;
	Cluster* allClusters = (Cluster*)malloc(sizeof(Cluster)*numOfClusters);
	for (i = 0; i < numOfClusters; i++)
	{
		allClusters[i].x = allPoints[i].x;
		allClusters[i].y = allPoints[i].y;
		allClusters[i].z = allPoints[i].z;
		allClusters[i].id = i;
		allClusters[i].clusterSize = 0;
		allClusters[i].diameter = 0;
	}
	return allClusters;
}

Point* readFile(int* numOfPoints, int* numOfClusters, int* T, double* dT, int* LIMIT, double* QM)
{
	FILE* f;
	int i;
	f = fopen(INPUT, "r");
	fscanf(f, "%d %d %d %lf %d %lf\n", numOfPoints, numOfClusters, T, dT, LIMIT, QM);
	Point* allPoints = (Point*)malloc(*numOfPoints * sizeof(Point));
	for (i = 0; i < *numOfPoints; i++)
	{
		fscanf(f, "%lf %lf %lf", &allPoints[i].x, &allPoints[i].y, &allPoints[i].z);
		fscanf(f, "%lf %lf %lf\n", &allPoints[i].vX, &allPoints[i].vY, &allPoints[i].vZ);
		allPoints[i].clusterId = -1;
		allPoints[i].minDis = -1;
	}
	fclose(f);
	return allPoints;
}

void writeToFile(Cluster* allClusters, int size, double interval, double quailty)
{
	FILE* f;
	int i;
	f = fopen(OUTPUT, "w");
	fprintf(f, "First occurrence at t = %lf with q = %lf\n", interval, quailty);
	fprintf(f, "Centers of the clusters:\n");
	for (i = 0; i < size; i++)
	{
		fprintf(f, "%lf %lf %lf\n", allClusters[i].x, allClusters[i].y, allClusters[i].z);
	}
	fclose(f);
}

void printClusters(Cluster* allClusters, int size)
{
	int i;
	for (i = 0; i < size; i++)
	{
		printf("%lf %lf %lf %d %d ", allClusters[i].x, allClusters[i].y, allClusters[i].z, allClusters[i].id, allClusters[i].clusterSize);
		printf("%lf\n", allClusters[i].diameter);
	}
}

void printPoints(Point* allPoints, int size)
{
	int i;
	for (i = 0; i < size; i++)
	{
		printf("%d : %lf %lf %lf\n", i + 1, allPoints[i].x, allPoints[i].y, allPoints[i].z);
	}
}