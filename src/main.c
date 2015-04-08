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



/* -------------------------------------------------------------------------- */
/* prototypes and global variables */
/*  */
double alpha(double a, double e, double theta);
double beta (double a, double e, double theta);
double gamma_sq(double a, double e, double theta);

#ifdef MPI_PARALLEL
static int rank, size;
#endif



/* -------------------------------------------------------------------------- */
/* main() */
/*  */
int main(int argc, char *argv[])
{
  /* make an n-by-n grid over eccentricity and theta spanning
     emin..emax and tmin..tmax
   */
  int i, j, n=128;
  double a=0.5, e, theta;

  double emin = 0.05,     emax = 1.0;
  double tmin = 0.1*M_PI, tmax = 0.9*M_PI;

  int imin=0, imax=n;           /* loop bounds for this processor */


  /* crude progress bar */
  char progress[64] = "[0%...25%...50%...75%...100%]", buf[64] = "";
  int ind, old=0;


  /* array to store the data */
  double ***data_array;
  FILE *fp;


  /* MPI stuff */
#ifdef MPI_PARALLEL
  int err;
  double *send_buf, *recv_buf;

  MPI_Init(&argc, &argv);                     /* starts MPI              */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);        /* get current process id  */
  MPI_Comm_size(MPI_COMM_WORLD, &size);        /* get number of processes */

  send_buf = (double*) calloc_1d_array(3*n*n, sizeof(double));
  recv_buf = (double*) calloc_1d_array(3*n*n, sizeof(double));

  if (n % size != 0)
    ath_error("n must be divisible by num processors,\n");

  /* distribute work by cutting in i (eccentricity) */
  /* this is better than cutting in j, but definitely not
     great... work is dominated by the low-eccentricity guys.  need to
     think more about load-balancing if this is really going to be
     efficient.
   */
  imin = rank * (n/size);
  imax = (rank+1) * (n/size);
#endif /* MPI_PARALLEL */


  /* allocate memory for data storage and populate it! */
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


  /* if running over MPI, combine the data.  do it as crudely and
     inefficiently as possible
   */
#ifdef MPI_PARALLEL
  /* pack into a 1d array for sending */
  for(i=imin; i<imax; i++){
    for(j=0; j<n; j++){
      ind = (i*n + j) * 3;
      send_buf[ind  ] = data_array[0][i][j];
      send_buf[ind+1] = data_array[1][i][j];
      send_buf[ind+2] = data_array[2][i][j];
    }
  }

  err = MPI_Allreduce(send_buf, recv_buf, 3*n*n, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);

  if (err != 0)
    ath_error("mpi_allreduce failed on proc %d with error %d", rank, err);

  /* unpack */
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      ind = (i*n + j) * 3;
      data_array[0][i][j] = recv_buf[ind  ];
      data_array[1][i][j] = recv_buf[ind+1];
      data_array[2][i][j] = recv_buf[ind+2];
    }
  }

  free_1d_array((void*) send_buf);
  free_1d_array((void*) recv_buf);
#endif /* MPI_PARALLEL */


  /* write out the data.  only on the first processor if running over
     MPI. */
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


  /* clean up and quit. */
  free_3d_array((void***) data_array);

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif

  return 0;
}
/* end main() */
/* -------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------- */
/* science functions */

/* gamma^2 = -alpha * beta.  include the factor of j^2 here for
   simplicity
 */
double gamma_sq(double a, double e, double theta)
{
  double al = alpha(a, e, theta),
         be =  beta(a, e, theta);
  return -1.0 * al * be / (a * (1.0-e*e));
}


/* number of points to evaluate tau(i) for determining the slope.
   only one is needed to get an answer; two needed for error
   checking.  use more for more stringent error checking, or to plot
   tau(i) for debugging
 */
#define nfit 2


/* alpha = d(tau.a)/d(ia) */
/*  */
double alpha(double a, double e, double theta)
{
  double tau[3];
  double lim = 0.1;

  int i;
  double ia[nfit], ta[nfit];

  double axis[3], slope, err;

  /* ahat for orbit 2 (= the rotated one)
   */
  axis[0] = -cos(theta);
  axis[1] = -sin(theta);
  axis[2] = 0.0;

  /* make a list of points {i, tau(i)}
   */
  for (i=0; i<nfit; i++)
    ia[i] = 0.0 + 1.0 * lim * (i+1) / (nfit);

  for (i=0; i<nfit; i++){
    torque_converged(a, e, theta, ia[i], 0.0, tau);
    ta[i] = dot(tau, axis);
  }

  /* estimate the slope and the error due to non-linearity
   */
  slope = (ta[nfit-1]) / (ia[nfit-1]);

  err = 0.0;
  for (i=0; i<nfit-1; i++)
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



/* beta = d(tau.b)/d(ib) */
/*  */
double beta(double a, double e, double theta)
{
  double tau[3];
  double lim = 0.1;

  int i;
  double ib[nfit], tb[nfit];

  double axis[3], slope, err;

  /* bhat for orbit 2 (= the rotated one)
   */
  axis[0] = -sin(theta);
  axis[1] =  cos(theta);
  axis[2] = 0.0;

  /* make a list of points {i, tau(i)}
   */
  for (i=0; i<nfit; i++)
    ib[i] = 0.0 + 1.0 * lim * (i+1) / (nfit);

  for (i=0; i<nfit; i++){
    torque_converged(a, e, theta, 0.0, ib[i], tau);
    tb[i] = dot(tau, axis);
  }

  /* estimate the slope and the error due to non-linearity
   */
  slope = (tb[nfit-1]) / (ib[nfit-1]);

  err = 0.0;
  for (i=0; i<nfit-1; i++)
    err += fabs(tb[i]-slope*ib[i]);

  if (err > 1.0e-2 * fabs(tb[nfit-1]))
    fprintf(stderr,
            "[ beta]: nonlinear fit at [a,e] = [%f,%f] (%f %%)\n",
            a, e, 100.0 * err / tb[nfit-1]);

  return slope;
}
