#ifndef EPI_H
#define EPI_H

#include "typedefs.h"
#include <optv/calibration.h>
#include <optv/tracking_frame_buf.h>
#include <optv/parameters.h>

typedef struct {
  int  	pnr;
  double  tol, corr;
} candidate;

double epi_line(double xl, double yl, Exterior Ex1, Interior I1, Glass G1,
    Exterior Ex2, Interior I2, Glass G2);
int  epi_mm(int cam, double xl, double yl, Exterior Ex1, Interior I1, Glass G1,
    Exterior Ex2, Interior I2, Glass G2, mm_np mmp, volume_par *vpar,
    double *xmin, double *ymin, double *xmax, double *ymax);
int  epi_mm_2D(double xl, double yl, Exterior Ex1, Interior I1, Glass G1,
    mm_np mmp, volume_par *vpar, double *xout, double *yout, double *zout);
void find_candidate_plus_msg(coord_2d crd[], target pix[], int num,
    double xa, double ya, double xb, double yb,
    int n, int nx, int ny, int sumg, candidate cand[], int *count, int i12,
    volume_par *vpar, control_par *cpar);
void find_candidate_plus(coord_2d crd[], target pix[], int num,
    double xa, double ya, double xb, double yb,
    int n, int nx, int ny, int sumg, candidate cand[], int *count, int nr,
    volume_par *vpar, control_par *cpar);

#endif
