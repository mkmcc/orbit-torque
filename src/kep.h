#ifndef KEP_H
#define KEP_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "ath_array.h"
#include "ath_error.h"
#include "linalg.h"

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
                  double ***kep, double **dm);



/* torque on orbit k1 due to orbit k2 */
/* - output is the (cartesian vector) torque, stored in the last
     argument
 */
void torque(double **k1, double *m1, double **k2, double *m2,
            int npts,
            double *torque);



/* torque on point p1 due to orbit kep */
/* - output is the (cartesian vector) torque, stored in the last
     argument
 */
void dtorque(double *p1, double dm1, double **kep, double *dm,
             int npts,
             double *torque);



/* helper function.  not meant to be called directly
 */
void torque_helper(double a, double e, double theta, double ia, double ib,
                   int npts,
                   double *tau);



/* torque between two orbits with the same semi-major axis a and
   eccentricity e, with eccentricity vectors separated by the angle
   theta.  one of the orbit is inclined by either ia or ib, following
   ann-marie's definitions.  repeatedly doubles the resolution until
   the torque converges.  result is stored as a cartesian vector in
   the final argument tau.
 */
void torque_converged(double a, double e, double theta, double ia, double ib,
                      double *tau);

#endif
