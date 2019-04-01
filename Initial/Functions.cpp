#include "Functions.h"

double* buildMomentsArray(int interval, double deltaInterval)
{
	int i;
	double* moments = (double*)malloc((int)(1 + interval / deltaInterval) * sizeof(double));
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i <= (int)(interval / deltaInterval); i++)
	{
		moments[i] = i*deltaInterval;
	}
	return moments;
}