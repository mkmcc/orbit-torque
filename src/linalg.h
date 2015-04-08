#ifndef LINALG_H
#define LINALG_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* all of the below assume 3D cartesian geometry */


/* return dot product of v1 and v2. */
double dot(double *v1, double *v2);

/* return the distance between points p1 and p2 */
double dist(double *p1, double *p2);

/* return the norm of vector vec */
double norm(double *vec);

/* compute cross product v1 x v2 and store in result */
void cross(double *v1, double *v2, double *result);


/* rotate vector vec through an angle theta about the axis axis.
   store result in the final argument result
*/
void rotate(double *axis, double theta, double *vec, double *result);

#endif
