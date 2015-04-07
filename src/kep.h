#ifndef KEP_H
#define KEP_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "ath_array.h"
#include "ath_error.h"
#include "linalg.h"

void make_ellipse(double a,  double e,  double alpha,
                  double ia, double ib, int npts, int shift,
                  double ***kep, double **dm);

void torque(double **k1, double *m1, double **k2, double *m2,
            int npts,
            double *torque);

void dtorque(double *p1, double dm1, double **kep, double *dm,
             int npts,
             double *torque);

void torque_helper(double a, double e, double theta, double ia, double ib,
                   int npts,
                   double *tau);

void torque_converged(double a, double e, double theta, double ia, double ib,
                      double *tau);

#endif
