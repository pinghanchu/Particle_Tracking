/*  Unit tests for functions related to the frame buffer. Uses the Check
    framework: http://check.sourceforge.net/
    
    To run it, type "make check" when in the top C directory, src_c/
    If that doesn't run the tests, use the Check tutorial:
    http://check.sourceforge.net/doc/check_html/check_3.html
*/

#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <optv/tracking_frame_buf.h>
#include <optv/calibration.h>
#include <optv/parameters.h>
#include <optv/trafo.h>
#include "../src_c/typedefs.h"

extern int init_proc_c(); /* Until main parameters are moved to parameters.h */
extern void img_coord(double X, double Y, double Z,
    Exterior Ex, Interior I, Glass G, ap_52 ap, mm_np mm, double *x, double *y);
extern n_tupel **correspondences(frame *frm, Calibration **calib, volume_par *vpar);
#include "../src_c/globals.h"

START_TEST(test_correspondences)
{
    frame frm;
    Calibration *calib[4];
    volume_par *vpar;
    n_tupel **corres_lists;
    mm_np media_par = {1, 1., {1., 0., 0.}, {1., 0., 0.}, 1., 1.};
    
    char ori_tmpl[] = "cal/sym_cam%d.tif.ori";
    char ori_name[25]; 
    
    int num_cams = 4, cam, cpt_vert, cpt_horz, cpt_ix;
    target *targ;
    
    int subset_size, num_corres;
    n_tupel *corres;
    
    ck_abort_msg("Known failure: j/p2 in find_candidate_plus breaks this.");
    chdir("testing_fodder/");
    init_proc_c();
    
    int i,j;
    /* Four cameras on 4 quadrants looking down into a calibration target.
       Calibration taken from an actual experimental setup */
    frame_init(&frm, num_cams, 16);
    for (cam = 0; cam < num_cams; cam++) {
        sprintf(ori_name, ori_tmpl, cam + 1);
        calib[cam] = read_calibration(ori_name, "cal/cam1.tif.addpar", NULL);
        
        /* Until calibration globals are removed: */
        Ex[cam] = calib[cam]->ext_par;
        I[cam] = calib[cam]->int_par;
        G[cam] = calib[cam]->glass_par;
        ap[cam] = calib[cam]->added_par;
        
        frm.num_targets[cam] = 16;
        
        /* Construct a scene representing a calibration target, generate
           tergets for it, then use them to reconstruct correspondences. */
        for (cpt_horz = 0; cpt_horz < 4; cpt_horz++) {
            for (cpt_vert = 0; cpt_vert < 4; cpt_vert++) {
                cpt_ix = cpt_horz*4 + cpt_vert;
                if (cam % 2) cpt_ix = 15 - cpt_ix; /* Avoid symmetric case */
                
                targ = &(frm.targets[cam][cpt_ix]);
                targ->pnr = cpt_ix;
                
                img_coord(cpt_vert * 10, cpt_horz * 10, 0, 
                    calib[cam]->ext_par, calib[cam]->int_par,
                    calib[cam]->glass_par, calib[cam]->added_par, media_par,
                    &(targ->x), &(targ->y));
                metric_to_pixel(&(targ->x), &(targ->y), targ->x, targ->y, cpar);
                
                /* These values work in check_epi, so used here too */
                targ->n = 25;
                targ->nx = targ->ny = 5;
                targ->sumg = 10;
            }
        }
    }
    vpar = read_volume_par("parameters/criteria.par");
    corres_lists = correspondences(&frm, calib, vpar);
    fail_unless(corres_lists[0] == NULL);
    
    for (subset_size = 2; subset_size <= num_cams; subset_size++) {
        corres = corres_lists[subset_size - 1];
        num_corres = 0;
        
        while (corres->corr >= 0) {
            num_corres++;
            corres++;
        }
        fail_unless(num_corres == ((subset_size == 4) ? 16 : 0));
    }
}
END_TEST

Suite* corres_suite(void) {
    Suite *s = suite_create ("Correspondences");

    TCase *tc = tcase_create ("Calibration target correspondences");
    tcase_add_test(tc, test_correspondences);
    suite_add_tcase (s, tc);

    return s;
}

int main(void) {
    int number_failed;
    Suite *s = corres_suite();
    SRunner *sr = srunner_create (s);
    srunner_run_all (sr, CK_ENV);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

