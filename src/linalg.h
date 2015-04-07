#ifndef LINALG_H
#define LINALG_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double dot (double *v1, double *v2);
double dist (double *v1, double *v2);
double norm (double *vec);
void cross (double *v1, double *v2, double *result);
void rotate (double *axis, double theta, double *vec, double *result);
void test_linalg();

#endif
