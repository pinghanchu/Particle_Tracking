/* Forward declarations for tracking-specific helper functions from ttools.c */

#ifndef TTOOLS_H
#define TTOOLS_H

#include "typedefs.h"
#include <optv/parameters.h>

int pix_in_next (target  next[], int num,
    double x, double y, double dl, double dr, double du, double dd, int found[POSI]);

int candsearch_in_pix(target  next[], int num, double x, double y,
    double dl, double dr, double du, double dd, int p[4], control_par *cpar);
int candsearch_in_pixrest(target  next[], int num, double x, double y,
    double dl, double dr, double du, double dd, int p[4], control_par *cpar);

void sortwhatfound (foundpix item[16], int *zaehler, int num_cams);
void searchquader(double X, double Y, double Z,
    double xr[4], double xl[4], double yd[4], double yu[4],
    track_par *tpar, control_par *cpar);
void predict(double x1, double y1, double x2, double y2,
    double *x3, double *y3);

#endif
