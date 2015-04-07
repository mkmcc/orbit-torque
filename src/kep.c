#include "kep.h"

void torque(double **k1, double *m1, double **k2, double *m2,
            int npts,
            double *torque)
{
  double dt[3];
  int i;

  torque[0] = torque[1] = torque[2] = 0.0;

  for (i=0; i<npts; i++){
    dtorque(k1[i], m1[i], k2, m2, npts, dt);

    torque[0] += dt[0];
    torque[1] += dt[1];
    torque[2] += dt[2];
  }

  return;
}

void dtorque(double *p1, double dm1, double **kep, double *dm,
             int npts,
             double *torque)
{
  double df[3], dr[3], dt[3];
  double fact, rsoft = 1.0e-10;
  int i;

  torque[0] = torque[1] = torque[2] = 0.0;

  for (i=0; i<npts; i++){
    dr[0] = p1[0] - kep[i][0];
    dr[1] = p1[1] - kep[i][1];
    dr[2] = p1[2] - kep[i][2];

    fact = dm1*dm[i] / pow(rsoft + norm(dr), 3);

    df[0] = dr[0] * fact;
    df[1] = dr[1] * fact;
    df[2] = dr[2] * fact;

    cross(p1, df, dt);

    torque[0] += dt[0];
    torque[1] += dt[1];
    torque[2] += dt[2];
  }

  return;
}


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

  tot = 0.0;
  for (i=0; i<npts; i++){
    mass[i] = dl[i] / v[i];
    tot += mass[i];
  }
  for (i=0; i<npts; i++)
    mass[i] /= tot;

  *dm  = mass;
  *kep = cartesian;


  /* -------------------------------------------------- */
  /* free memory and return */
  free_1d_array((void*) eccanomaly);
  free_1d_array((void*) trueanomaly);
  free_1d_array((void*) r);
  free_1d_array((void*) v);
  free_1d_array((void*) dl);

  /* do not free these! */
  /* free_1d_array((void*) mass); */
  /* free_2d_array((void**) cartesian); */

  return;
}

void torque_helper(double a, double e, double theta, double ia, double ib,
                   int npts,
                   double *tau)
{
  double **kep1, **kep2, *dm1, *dm2;

  make_ellipse(a, e, 0.0,    ia,  ib, npts, 0, &kep1, &dm1);
  make_ellipse(a, e, theta, 0.0, 0.0, npts, 1, &kep2, &dm2);

  torque(kep1, dm1, kep2, dm2, npts, tau);

  free_2d_array((void**) kep1);
  free_2d_array((void**) kep2);

  free_1d_array((void*) dm1);
  free_1d_array((void*) dm2);
}

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

    if (err <= 1.0e-2 * norm(tau))
      break;
  }

  if (i>imax){
    fprintf(stderr,
            "[torque_converged]: did not converge at [a,e,theta,ia,ib] = [%f,%f,%f,%f,%f]\n",
            a, e, theta, ia, ib);
  }

  return;
}
