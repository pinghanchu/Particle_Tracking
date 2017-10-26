/* Forward declarations for functions in tools.c should anyone need them */

#ifndef TOOLS_H
#define TOOLS_H

#include <optv/tracking_frame_buf.h>
#include "typedefs.h"

FILE *fopen_r(char *filename);

void quicksort_con(n_tupel *con, int num);
void sort(int n, float a[], int b[]);
void quicksort_target_y(target *pix, int num);
void quicksort_coord2d_x(coord_2d *crd, int num);
int kill_in_list (int nr, int num, int ms_x, int ms_y);

int nearest_neighbour_geo(coord_2d crd[], int num, double x, double y,
    double eps);

#endif
