
#ifndef IMG_COORD_H
#define IMG_COORD_H

#include <optv/calibration.h>
#include <optv/parameters.h>

void img_coord (int cam, double X, double Y, double Z, 
    Exterior Ex, Interior I, Glass G, ap_52 ap, mm_np mm, double *x, double *y);

#endif

