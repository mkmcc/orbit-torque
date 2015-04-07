#include "linalg.h"

/* return dot product of v1 and v2 */
double dot (double *v1, double *v2)
{
  int i;
  double sum = 0.0;
  for (i=0; i<3; i++)
    sum += v1[i] * v2[i];

  return sum;
}

double norm (double *vec)
{
  return (sqrt(dot(vec, vec)));
}

double dist (double *vec1, double *vec2)
{
  double tmp[3];

  tmp[0] = vec1[0] - vec2[0];
  tmp[1] = vec1[1] - vec2[1];
  tmp[2] = vec1[2] - vec2[2];

  return norm(tmp);
}

/* compute cross product v1 x v2 and store in result */
void cross (double *v1, double *v2, double *result)
{
  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return;
}

void rotate (double *axis, double theta, double *vec, double *result)
{
  int i, j;

  double c = cos(theta), s = sin(theta);
  double n = norm(axis);
  double wx = axis[0]/n, wy = axis[1]/n, wz = axis[2]/n;

  double mat[3][3];

  mat[0][0] = c + wx*wx*(1.0-c);
  mat[0][1] = wx*wy*(1.0-c) - wz*s;
  mat[0][2] = wy*s + wx*wz*(1.0-c);

  mat[1][0] = wz*s + wx*wy*(1.0-c);
  mat[1][1] = c + wy*wy*(1.0-c);
  mat[1][2] = -wx*s + wy*wz*(1.0-c);

  mat[2][0] = -wy*s + wx*wz*(1.0-c);
  mat[2][1] = wx*s+wy*wz*(1.0-c);
  mat[2][2] = c + wz*wz*(1.0-c);

  for (j=0; j<3; j++){
    result[j] = 0.0;
    for (i=0; i<3; i++){
      result[j] += mat[j][i] * vec[i];
    }
  }

  return;
}

void test_linalg()
{
  double xhat[3] = {1.0, 0.0, 0.0};
  double yhat[3] = {0.0, 1.0, 0.0};
  double zhat[3] = {0.0, 0.0, 1.0};

  double test[3];

  if (dot(xhat, xhat) == 1.0 &&
      dot(yhat, yhat) == 1.0 &&
      dot(zhat, zhat) == 1.0 &&
      dot(xhat, yhat) == 0.0 &&
      dot(xhat, zhat) == 0.0 &&
      dot(yhat, zhat) == 0.0)
    printf("dot product ok\n");

  cross(xhat, yhat, test);
  if (test[0] == 0.0 && test[1] == 0.0 && test[2] == 1.0)
    printf("cross product ok\n");

  double v1[3] = {0.947583, 0.117293, -0.203023};
  double v2[3] = {-1.04086, 1.14846, -0.131156};

  double theta = 2.95772;

  rotate(v2, theta, v1, test);
  printf("%f %f %f\n",
         test[0] - (-0.253408),
         test[1] - (-0.931318),
         test[2] - (0.146014));


  double w1[3] = {0.15087, 1.02944, 2.38964};
  double w2[3] = {0.593419, 0.309299, 0.719491};

  theta = 6.27989;

  rotate(w2, theta, w1, test);
  printf("%f %f %f\n",
         test[0] - (0.150881),
         test[1] - (1.03383),
         test[2] - (2.38774));

  return;
}
