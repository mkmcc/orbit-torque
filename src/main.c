#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "linalg.h"
#include "ath_array.h"
#include "ath_error.h"
#include "kep.h"

#define MPI_PARALLEL

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

double alpha(double a, double e, double theta);
double beta (double a, double e, double theta);
double gamma_sq(double a, double e, double theta);

#ifdef MPI_PARALLEL
static int rank, size;
#endif

int main (int argc, char *argv[])
{
  int i, j, n=32;
  double a=0.5, e, theta;

  double emin = 0.5, emax = 1.0;
  double tmin = 0.0, tmax = 1.0*M_PI;

  /* crude progress bar */
  char progress[64] = "[0%...25%...50%...75%...100%]", buf[64] = "";
  int ind, old=0;

  /* data array */
  double ***data_array;

  FILE *fp;

  int imin=0, imax=n;

#ifdef MPI_PARALLEL
  int err;
  double ***global_data_array;

  MPI_Init (&argc, &argv);                      /* starts MPI              */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id  */
  MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of
                                                   processes */

  global_data_array = (double***) calloc_3d_array(3, n, n, sizeof(double));

  if (n % size != 0)
    ath_error("n must be divisible by num processors,\n");

  imin = rank * (n/size);
  imax = (rank+1) * (n/size);
#endif /* MPI_PARALLEL */

  data_array = (double***) calloc_3d_array(3, n, n, sizeof(double));


  for(i=imin; i<imax; i++){
    for(j=0; j<n; j++){
      e     = emin + (i+0.5)/n * (emax-emin);
      theta = tmin + (j+0.5)/n * (tmax-tmin);

      data_array[0][i][j] = e;
      data_array[1][i][j] = theta;
      data_array[2][i][j] = gamma_sq(a, e, theta);

      /* crude progress bar */
      ind = floor( (1.0*(i-imin)*n+j+1)/((imax-imin)*n) * strlen(progress));
      if (ind > old){
        old = ind;
        strncpy(buf, progress, ind);
        buf[ind] = '\0';
#ifdef MPI_PARALLEL
        fprintf(stderr, "[proc %d]: %s\n", rank, buf);
#else
        fprintf(stderr, "%s\n", buf);
#endif
      }
    }
  }

#ifdef MPI_PARALLEL
  err = MPI_Reduce(data_array, global_data_array, 3*n*n, MPI_DOUBLE,
                   MPI_SUM, 0, MPI_COMM_WORLD);
#endif /* MPI_PARALLEL */


#ifdef MPI_PARALLEL
  if (rank == 0){
#endif
  fp = fopen("gamma-sq.dat", "w");
  for (i=0; i<n; i++){
    for (j=0; j<n; j++){
      fprintf(fp, "%f\t%f\t%f\n",
              data_array[0][i][j],
              data_array[1][i][j],
              data_array[2][i][j]);
    }
  }
  fclose(fp);
#ifdef MPI_PARALLEL
  }
#endif

  free_3d_array((void***) data_array);

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif

  return 0;
}

double gamma_sq(double a, double e, double theta)
{
  double al = alpha(a, e, theta),
    be = beta(a, e, theta);
  return -1.0 * al * be / (a * (1.0-e*e));
}

double alpha(double a, double e, double theta)
{
  double tau[3];
  double lim = 0.1;

  int i, nfit=2;
  double ia[nfit], ta[nfit];

  double axis[3], slope, err;

  /* ahat for orbit 2 (= the rotated one) */
  axis[0] = -cos(theta);
  axis[1] = -sin(theta);
  axis[2] = 0.0;

  for (i=0; i<nfit; i++)
    ia[i] = 0.0 + 1.0 * lim * (i+1) / (nfit);

  for (i=0; i<nfit; i++){
    torque_converged(a, e, theta, ia[i], 0.0, tau);
    ta[i] = dot(tau, axis);
  }

  slope = (ta[nfit-1] - ta[0]) / (ia[nfit-1] - ia[0]);

  err = 0.0;
  for (i=1; i<nfit-1; i++)
    err += fabs(ta[i]-slope*ia[i]);

  if (err > 1.0e-2 * fabs(ta[nfit-1])){
    fprintf(stderr,
            "[alpha]: nonlinear fit at [a,e] = [%f,%f] (%f %%)\n",
            a, e, 100.0 * err / ta[nfit-1]);

    for(i=0; i<nfit; i++)
      fprintf(stderr, "{%f,\t%f},\n", ia[i], ta[i]);
  }

  return slope;
}

double beta(double a, double e, double theta)
{
  double tau[3];
  double lim = 0.1;

  int i, nfit=2;
  double ib[nfit], tb[nfit];

  double axis[3], slope, err;

  /* bhat for orbit 2 (= the rotated one) */
  axis[0] = -sin(theta);
  axis[1] =  cos(theta);
  axis[2] = 0.0;

  for (i=0; i<nfit; i++)
    ib[i] = 0.0 + 1.0 * lim * (i+1) / (nfit);

  for (i=0; i<nfit; i++){
    torque_converged(a, e, theta, 0.0, ib[i], tau);
    tb[i] = dot(tau, axis);
  }

  slope = (tb[nfit-1] - tb[0]) / (ib[nfit-1] - ib[0]);

  err = 0.0;
  for (i=1; i<nfit-1; i++)
    err += fabs(tb[i]-slope*ib[i]);

  if (err > 1.0e-2 * fabs(tb[nfit-1]))
    fprintf(stderr,
            "[ beta]: nonlinear fit at [a,e] = [%f,%f] (%f %%)\n",
            a, e, 100.0 * err / tb[nfit-1]);

  return slope;
}
