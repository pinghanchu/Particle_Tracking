/* Forward declarations for various image processing routines collected from
image_processing.c, segmentation.c and peakfitting.c */

#ifndef IMAGE_PROC_C
#define IMAGE_PROC_C

#include <optv/parameters.h>

int peak_fit_new (unsigned char	*img, char par_file[],
    int xmin, int xmax, int ymin, int ymax, target pix[], int nr, 
    control_par *cpar);

void simple_connectivity(unsigned char *img0,
    unsigned char *img, char par_file[],
    int xmin, int xmax, int ymin, int ymax,
    target pix[], int nr, int *num, control_par *cpar);

void targ_rec (unsigned char *img0, unsigned char *img, char par_file[], 
    int xmin, int xmax, int ymin, int ymax,
    target pix[], int nr, int *num, control_par *cpar);

void highpass(unsigned char *img, unsigned char *img_hp, int dim_lp,
    int filter_hp, control_par *cpar);

#endif

