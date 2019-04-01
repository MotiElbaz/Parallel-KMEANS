#pragma once

#include "mpi.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Structs.h"

#define INPUT "input.txt"
#define OUTPUT "output.txt"
#define MASTER 0
#define TRUE 1
#define FALSE 0

double* buildMomentsArray(int interval, double deltaInterval);