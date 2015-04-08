#include "linalg.h"

/* all of the below assume 3D cartesian geometry */

/* return dot product of v1 and v2. */
double dot(double *v1, double *v2)
{
  int i;
  double sum = 0.0;
  for (i=0; i<3; i++)
    sum += v1[i] * v2[i];

  return sum;
}

/* return the norm of vector vec */
double norm(double *vec)
{
  return (sqrt(dot(vec, vec)));
}

/* return the distance between points pt1 and pt2 */
double dist(double *pt1, double *pt2)
{
  double tmp[3];

  tmp[0] = pt1[0] - pt2[0];
  tmp[1] = pt1[1] - pt2[1];
  tmp[2] = pt1[2] - pt2[2];

  return norm(tmp);
}

/* compute cross product v1 x v2 and store in result */
void cross(double *v1, double *v2, double *result)
{
  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return;
}


/* rotate vector vec through an angle theta about the axis axis.
   store result in the final argument result
 */
void rotate(double *axis, double theta, double *vec, double *result)
{
  int i, j;

  double c = cos(theta), s = sin(theta);
  double n = norm(axis);
  double wx = axis[0]/n, wy = axis[1]/n, wz = axis[2]/n;

  double mat[3][3];

  /* build the rotation matrix */
  /* wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle */
  /*  */
  mat[0][0] = c + wx*wx*(1.0-c);
  mat[0][1] = wx*wy*(1.0-c) - wz*s;
  mat[0][2] = wy*s + wx*wz*(1.0-c);

  mat[1][0] = wz*s + wx*wy*(1.0-c);
  mat[1][1] = c + wy*wy*(1.0-c);
  mat[1][2] = -wx*s + wy*wz*(1.0-c);

  mat[2][0] = -wy*s + wx*wz*(1.0-c);
  mat[2][1] = wx*s+wy*wz*(1.0-c);
  mat[2][2] = c + wz*wz*(1.0-c);

  /* multiply vector by the rotation matrix */
  for (j=0; j<3; j++){
    result[j] = 0.0;
    for (i=0; i<3; i++){
      result[j] += mat[j][i] * vec[i];
    }
  }

  return;
}
