
#ifndef SHARPDATA_MONOMAT_H
#define SHARPDATA_MONOMAT_H

#include <R.h>
#include <Rinternals.h>

void F77_SUB(kern)(double *x, double *xi, double *h, double *kernel, 
    double *dkernel, double *d2kernel);
void F77_SUB(afun)(double *xobs, double *xgrid, int *n, int *m, 
    double *amat, double *loclin, double *dloclin, double *h, int *d);

#endif /* SHARPDATA_MONOMAT_H */

