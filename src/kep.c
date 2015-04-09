#include "kep.h"

/* torque on orbit k1 due to orbit k2 */
/* - output is the (cartesian vector) torque, stored in the last
     argument
 */
void torque(double **k1, double *m1, double **k2, double *m2,
            int npts,
            double *torque)
{
  double dt[3];
  int i;

  torque[0] = torque[1] = torque[2] = 0.0;

  /* loop over points in k1 and sum differential torques */
  for (i=0; i<npts; i++){
    dtorque(k1[i], m1[i], k2, m2, npts, dt);

    torque[0] += dt[0];
    torque[1] += dt[1];
    torque[2] += dt[2];
  }

  return;
}



/* torque on point p1 due to orbit kep */
/* - output is the (cartesian vector) torque, stored in the last
     argument
 */
void dtorque(double *p1, double dm1, double **kep, double *dm,
             int npts,
             double *torque)
{
  double df[3], dr[3], dt[3];
  double fact, rsoft = 1.0e-10;
  int i;

  torque[0] = torque[1] = torque[2] = 0.0;

  /* loop over points in kep and sum differential torques */
  for (i=0; i<npts; i++){
    /* dr vector points *to* p1 */
    dr[0] = p1[0] - kep[i][0];
    dr[1] = p1[1] - kep[i][1];
    dr[2] = p1[2] - kep[i][2];

    fact = dm1*dm[i] / pow(rsoft + norm(dr), 3);

    /* df is the force that p1 has on kep[i] */
    /* (ie, points toward p1) */
    df[0] = dr[0] * fact;
    df[1] = dr[1] * fact;
    df[2] = dr[2] * fact;

    /* dt = r cross f */
    /* maybe should be using -f here?  don't think it matters... */
    cross(p1, df, dt);

    torque[0] += dt[0];
    torque[1] += dt[1];
    torque[2] += dt[2];
  }

  return;
}



/* make a set of points representing a keplerian orbit */
/* - points are evenly spaced in eccentric anomaly */
/* - inputs:
     1. a, e: semi-major axis and eccentricity
     2. alpha: rotation angle about the z-axis.  (alpha=0 corresponds
        to an eccentricity vector pointing along the -x direction)
     3. ia, ib: ann-marie's definitions.  only one can be non-zero
     4. npts: number of points spaced along the orbit
     5. shift: shift points by 1/2 spacing.  helps with convergence.

  - outputs:
    1. kep: list of cartesian coordinates for the points.  this is a
       2D array.  argument is *** because it's a pointer to the
       array, which is admittedly confusing.
    2. dm: list of masses attached to each point.  pointer to a 1D
       array, which is **
 */
void make_ellipse(double a,  double e,  double alpha,
                  double ia, double ib, int npts, int shift,
                  double ***kep, double **dm)
{
  double *eccanomaly, *trueanomaly, *r, *v, *dl;
  double **cartesian, *mass;
  int i;
  double axis[3], old[3], rot;
  double tot;


  /* -------------------------------------------------- */
  /* allocate memory */
  cartesian   = (double**) calloc_2d_array(npts, 3, sizeof(double));

  mass        = (double*)  calloc_1d_array(npts, sizeof(double));

  /* these are temporary arrays */
  eccanomaly  = (double*)  calloc_1d_array(npts, sizeof(double));
  trueanomaly = (double*)  calloc_1d_array(npts, sizeof(double));
  r           = (double*)  calloc_1d_array(npts, sizeof(double));
  v           = (double*)  calloc_1d_array(npts, sizeof(double));
  dl          = (double*)  calloc_1d_array(npts, sizeof(double));


  /* -------------------------------------------------- */
  /* preliminaries: make points in polar coords */

  /* make the eccentric anomaly */
  if (shift > 0){
    for (i=0; i<npts; i++)
      eccanomaly[i] = 2.0 * M_PI * i/npts;
  } else{
    for (i=0; i<npts; i++)
      eccanomaly[i] = 2.0 * M_PI * (i+0.5)/npts;
  }

  /* make the true anomaly */
  for (i=0; i<npts; i++){
    trueanomaly[i] = 2.0 * atan2(sqrt(1.0-e) * cos(eccanomaly[i]/2),
                                 sqrt(1.0+e) * sin(eccanomaly[i]/2));
  }

  /* make r and v */
  for (i=0; i<npts; i++){
    r[i] = a*(1.0-e*e)/(1.0-e*cos(trueanomaly[i]));
    v[i] = sqrt(2.0/r[i] - 1.0/a);
  }


  /* -------------------------------------------------- */
  /* make the cartesian points, assuming eccentricity vector points
     along negative x-axis, and b vector points along negative
     y-axis */
  for(i=0; i<npts; i++){
    cartesian[i][0] = r[i] * cos(trueanomaly[i]);
    cartesian[i][1] = r[i] * sin(trueanomaly[i]);
    cartesian[i][2] = 0.0;
  }


  /* tilt along ia or ib */
  if (fabs(ia) > 0.0 && fabs(ib) > 0.0)
    ath_error("[make_ellipse]: ia and ib cannot both be nonzero.\n");

  if (fabs(ia) > 0.0){
    axis[0] = -1.0;
    axis[1] =  0.0;
    axis[2] =  0.0;

    rot = ia;
  }
  else if (fabs(ib) > 0.0){
    axis[0] =  0.0;
    axis[1] = -1.0;
    axis[2] =  0.0;

    rot = ib;
  }
  if (fabs(ia) > 0.0 || fabs(ib) > 0.0){
    for (i=0; i<npts; i++){
      old[0] = cartesian[i][0];
      old[1] = cartesian[i][1];
      old[2] = cartesian[i][2];

      rotate(axis, rot, old, cartesian[i]);
    }
  }

  /* rotate through alpha about the z-axis */
  if (fabs(alpha) > 0.0){
    axis[0] = 0.0;
    axis[1] = 0.0;
    axis[2] = 1.0;

    for (i=0; i<npts; i++){
      old[0] = cartesian[i][0];
      old[1] = cartesian[i][1];
      old[2] = cartesian[i][2];

      rotate(axis, alpha, old, cartesian[i]);
    }
  }


  /* -------------------------------------------------- */
  /* find the mass of each segment */

  /* start with edge cases */
  dl[0]      = dist(cartesian[1], cartesian[npts-1]);
  dl[npts-1] = dist(cartesian[0], cartesian[npts-2]);

  /* make distances from a symmetric difference */
  for (i=1; i<npts-1; i++)
    dl[i] = dist(cartesian[i+1], cartesian[i-1]);

  /* dm = dl/v; normalize total mass to 1 */
  tot = 0.0;
  for (i=0; i<npts; i++){
    mass[i] = dl[i] / v[i];
    tot += mass[i];
  }
  for (i=0; i<npts; i++)
    mass[i] /= tot;


  /* -------------------------------------------------- */
  /* store outputs in kep and dm */
  *dm  = mass;
  *kep = cartesian;


  /* -------------------------------------------------- */
  /* free memory and return */
  free_1d_array((void*) eccanomaly);
  free_1d_array((void*) trueanomaly);
  free_1d_array((void*) r);
  free_1d_array((void*) v);
  free_1d_array((void*) dl);

  /* do *not* free these! */
  /* free_1d_array((void*) mass); */
  /* free_2d_array((void**) cartesian); */

  return;
}



/* helper function used by torque_converged.  not meant to be called
   directly
 */
void torque_helper(double a, double e, double theta, double ia, double ib,
                   int npts,
                   double *tau)
{
  double **kep1, **kep2, *dm1, *dm2;

  make_ellipse(a, e, 0.0,    ia,  ib, npts, 0, &kep1, &dm1);
  /* rotate kep2 to some angle, and tilt kep1 by some inclination */
  make_ellipse(a, e, theta, 0.0, 0.0, npts, 1, &kep2, &dm2);

  /* torque on kep2 due to kep1 */
  torque(kep2, dm2, kep1, dm1, npts, tau);

  free_2d_array((void**) kep1);
  free_2d_array((void**) kep2);

  free_1d_array((void*) dm1);
  free_1d_array((void*) dm2);
}



/* torque between two orbits with the same semi-major axis a and
   eccentricity e, with eccentricity vectors separated by the angle
   theta.  one of the orbit is inclined by either ia or ib, following
   ann-marie's definitions.  repeatedly doubles the resolution until
   the torque converges.  result is stored as a cartesian vector in
   the final argument tau.
 */
void torque_converged(double a, double e, double theta, double ia, double ib,
                      double *tau)
{
  double oldtau[3], err;
  int j;

  int i=7, npts=128, imax=13;
  torque_helper(a, e, theta, ia, ib, npts, oldtau);

  for (; i<=imax; i++){
    npts *= 2;
    torque_helper(a, e, theta, ia, ib, npts, tau);

    err = 0.0;
    for (j=0; j<3; j++)
      err += fabs(tau[j]-oldtau[j]);
    for (j=0; j<3; j++)
      oldtau[j] = tau[j];

    if (err <= 1.0e-4 * norm(tau))
      break;
  }

  if (i>imax){
    fprintf(stderr,
            "[torque_converged]: did not converge at [a,e,theta,ia,ib] = [%f,%f,%f,%f,%f]\n",
            a, e, theta, ia, ib);
  }

  return;
}
