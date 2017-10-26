/* Forward declarations for calibration (camera orientation) functions found in
orientation.c */

#ifndef ORIENTATION_H
#define ORIENTATION_H

#include <optv/parameters.h>
#include <optv/vec_utils.h>
#include <optv/calibration.h>

void prepare_eval(control_par *cpar, int *n_fix);
void prepare_eval_shake(control_par *cpar);

void num_deriv_exterior(int cam, Exterior ext, Interior I0, Glass G0, ap_52 ap0,
    mm_np mm, double dpos, double dang, vec3d pos,
    double x_ders[6], double y_ders[6]);

int orient_v3(Exterior Ex0, Interior I0, Glass G0, ap_52 ap0,
    mm_np mm, int nfix, coord_3d fix[], coord_2d crd[], 
    Exterior *Ex, Interior *I, Glass *G, ap_52 *ap, int nr,
    double resid_x[], double resid_y[], int pixnr[], int *num_used);
int raw_orient_v3(Exterior Ex0, Interior I, Glass G0, ap_52 ap,
    mm_np mm, int nfix, coord_3d fix[], coord_2d crd[], 
    Exterior *Ex, Glass *G, int nr, int only_show);
void orient_v5(control_par *cpar, int nfix, Calibration cal[]);

#endif
