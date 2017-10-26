/* Finding correspondences between particles in a set of images. 

References:
[1] http://en.wikipedia.org/wiki/Tree_traversal#Depth-first
*/

#include <stdlib.h>
#include <stdio.h>

#include "typedefs.h"
#include "globals.h"
#include "epi.h"
#include "btree.h"

#include <optv/calibration.h>
#include <optv/parameters.h>
#include <optv/trafo.h>

/*  advance_nd_iterator() increments by 1 an iterator over n dimensions. The
    iterator is an array of counters, one for each dimension. Least significant
    dimension is on the left. 

    Arguments:
    int nd - number of dimensions for this iterator.
    int *current - pointer to the iterator.
    int *count - pointer to array of the same size holding the maximal value
        for each dimension.
    
    Returns:
    1 if the iterator has finished and was reset to zeroes, 0 otherwise.
*/
int advance_nd_iterator(int nd, int *current, int *count) {
    int d, carry = 1;
    
    for (d = 0; d < nd; d++) {
        current[d] += carry;
        if (current[d] >= count[d]) {
            current[d] = 0;
            carry = 1;
        } else {
            carry = 0;
            break;
        }
    }
    return carry;
}

/* correspondences() generates a list of correspondence groups. Each
   correspondence group is the set of target numbers identifying each one
   target in one image from the set, such that each target in the cgroup is an
   image of the same 3D particle.

   Arguments:
   frame *frm - holds all needed image data.
   calibration **calib - holds camera calibration information needed for getting
       real metric target coordinates corrected from those identified in image.
       Array of num_cams pointers to calibration objects.
   volume_par *vpar - search volume parameters, contains the search volume for
       epipolar lines and the measure of correspondence needed for a candidate
       target.
   
   Returns:
   n_tupel **corres_lists - a num_cams list of corres-lists, each is a list of 
        4-item arrays containing the index of target for each particle in the
        list in each of up to 4 images, in order. Each list is terminated by an
        n_tupel with negative .corr and garbage p[]. The top-level array is by
        number of particles in correspondence, so that the first element is
        always NULL and is there for convenience only.
*/
n_tupel **correspondences(frame *frm, Calibration **calib, volume_par *vpar) {
    int chfield = 0; /* temporary replacement of obsolete global */
    int cam, part; /* loop counters */
    double img_x, img_y; /* image center */
    coord_2d **corrected;
    
    int part_img, epi_img;
    double epi_start[2], epi_end[2]; /* each for x,y coordinates of a point */
    target *targ; /* shortcut to working target. */
    candidate cand[maxcand];
    int num_cands, cix, adj_row, adj_ix;
    int *accum_parts;  /* accumulated number of particles through each image */
    int *cam_set, *cam_subset, *subset_count, *subset_iter;
    candidate *adjacency; /* Adjacency matrix. */
    
    int subset, subset_size, pow_set_size;
    double corr, dist; /* accumulators for set correspondences calculation. */
    int is_clique, *num_matches;
    btree_node **match_trees;
    n_tupel *match;
    
    int stack_top = 0; /* For iterative tree traversal. */
    btree_node *node;
    btree_node **stack = NULL;
    n_tupel **corres_lists, *flat_top;
    int **part_used, use_corres;
    int list_len;
    
    /* Hidden first element stores 0, to save some branches below. */
    accum_parts = (int *) malloc(frm->num_cams * sizeof(int) + 1);
    accum_parts[0] = 0;
    accum_parts++;
    
    cam_set = (int *) malloc(frm->num_cams * sizeof(int));
    
    /* We work on distortion-corrected image coordinates of particles.
       This loop does the correction. It also recycles the iteration on
       frm->num_cams to allocate some arrays needed later and do some related
       preparation. */
    corrected = (coord_2d **) malloc(frm->num_cams * sizeof(coord_2d *));
    for (cam = 0; cam < frm->num_cams; cam++) {
        corrected[cam] = (coord_2d *) malloc(
            frm->num_targets[cam] * sizeof(coord_2d));
        
        accum_parts[cam] = frm->num_targets[cam];
        if (cam > 0) accum_parts[cam] += accum_parts[cam - 1];
        cam_set[cam] = 1 << cam;
        
        for (part = 0; part < frm->num_targets[cam]; part++) {
            pixel_to_metric(&corrected[cam][part].x, &corrected[cam][part].y,
                frm->targets[cam][part].x, frm->targets[cam][part].y, cpar);
            
            img_x = corrected[cam][part].x - calib[cam]->int_par.xh;
            img_y = corrected[cam][part].y - calib[cam]->int_par.yh;
            
            correct_brown_affin (img_x, img_y, calib[cam]->added_par,
               &corrected[cam][part].x, &corrected[cam][part].y);
            
            corrected[cam][part].pnr = frm->targets[cam][part].pnr;
        }
        
        /* This is expected by find_candidate_plus() */
        quicksort_coord2d_x(corrected[cam], frm->num_targets[cam]);
    }

    /* Build adjacency matrix of particles in all images to each-other, using a
       list of candidates on the epipolar line from particle in one image to a
       paired image. A link exists if a cell's .corr attribute != 0.
    */
    adjacency = (candidate *) calloc( /* calloc() zeroes memory in advance */
        accum_parts[frm->num_cams - 1] * accum_parts[frm->num_cams - 1],
        sizeof(candidate));
    
    for (part_img = 0; part_img < frm->num_cams - 1; part_img++) {
        for (epi_img = part_img + 1; epi_img < frm->num_cams; epi_img++) {
            /* Note that we always take the particle from the first image
               and candidates from the second. The reverse may give slightly
               different results. A better algorithm would check both, but
               ambition must wait for now. */
            
            for (part = 0; part < frm->num_targets[part_img]; part++) {
                /* Find epipolar line on corrected image */
                epi_mm(epi_img, corrected[part_img][part].x, corrected[part_img][part].y,
                    calib[part_img]->ext_par, calib[part_img]->int_par,
                    calib[part_img]->glass_par, calib[epi_img]->ext_par, 
                    calib[epi_img]->int_par, calib[epi_img]->glass_par, *(cpar->mm), vpar,
                    &(epi_start[0]), &(epi_start[1]), &(epi_end[0]), &(epi_end[1]));
                
                /* Find candidates close to epipolar line */
                targ = &(frm->targets[part_img][corrected[part_img][part].pnr]);
                find_candidate_plus (corrected[epi_img], frm->targets[epi_img],
                    frm->num_targets[epi_img], epi_start[0], epi_start[1],
                    epi_end[0], epi_end[1], targ->n, targ->nx, targ->ny,
                    targ->sumg, cand, &num_cands, epi_img, vpar, cpar);
                
                /* Copy candidate information to the adjacency matrix. */
                if (num_cands > maxcand) num_cands = maxcand;
                adj_row = (accum_parts[part_img - 1] + corrected[part_img][part].pnr) *
                    accum_parts[frm->num_cams - 1];
                
                for (cix = 0; cix < num_cands; cix++) {
                    adj_ix =  adj_row + accum_parts[epi_img - 1] + cand[cix].pnr;

                    adjacency[adj_ix].pnr = cand[cix].pnr;
                    adjacency[adj_ix].corr = cand[cix].corr;
                    adjacency[adj_ix].tol = cand[cix].tol;
                }
            }
        }
    }
    
    /* Image coordinates not needed beyond this point. */
    for (cam = 0; cam < frm->num_cams; cam++) {
        free(corrected[cam]);
    }
    free(corrected);
    
    /* Match candidates between sets of (n > 1) images. Iteration is over the
       subset of the power-set of 0 .. num_cams - 1, i.e. all combinations of
       2, ..., num_cams members of the camera indices set (unordered
       combinations).
    */
    match_trees = (btree_node **) calloc(frm->num_cams, sizeof(btree_node *));
    num_matches = (int *) calloc(frm->num_cams, sizeof(int));
    cam_subset = (int *) malloc(frm->num_cams * sizeof(int));
    subset_count = (int *) malloc(frm->num_cams * sizeof(int));
    subset_iter = (int *) malloc(frm->num_cams * sizeof(int));
    
    pow_set_size = (1 << frm->num_cams) - 1;
    for (subset = 1; subset <= pow_set_size; subset++) {
        /* Prepare the subset for iterating over its particle groups. */
        subset_size = 0;
        for (cam = 0; cam < frm->num_cams; cam++) {
            if (subset & cam_set[cam]) {
                cam_subset[subset_size] = cam;
                subset_iter[subset_size] = 0;
                subset_count[subset_size] = frm->num_targets[cam];
                subset_size++;
            }
        }
        if (subset_size < 2) continue;
        
        do { /* For each combination of particles, one from each image in set */
            is_clique = 1;
            corr = dist = 0;
            /* Check that the adjacency matrix of the subset represents a
               complete subgraph by checking that the upper triangle is all
               ones. */
            for (part_img = 0; part_img < subset_size - 1; part_img++) {
                adj_row = (accum_parts[cam_subset[part_img] - 1] + \
                    subset_iter[part_img]) * accum_parts[frm->num_cams - 1];
                
                for (epi_img = part_img + 1; epi_img < subset_size; epi_img++) {
                    adj_ix =  adj_row + accum_parts[cam_subset[epi_img] - 1] + \
                        subset_iter[epi_img];
                    if (adjacency[adj_ix].corr == 0) {
                        is_clique = 0;
                        break;
                    }
                    corr += adjacency[adj_ix].corr;
                    dist += adjacency[adj_ix].tol;
                }
                if (!is_clique) break;
            }
            if ((!is_clique) || (corr/dist < vpar->corrmin)) continue;
            corr /= dist;
            
            /* record a match */
            match = (n_tupel *) malloc(sizeof(n_tupel));
            for (cam = 0; cam < 4; cam++) match->p[cam] = -2;
            /* n_tupel is hard-coded to max. 4 cameras for now. */
            
            match->corr = corr;
            for (cam = 0; cam < subset_size; cam++) {
                match->p[cam_subset[cam]] = subset_iter[cam];
            }
            btree_insert(&match_trees[subset_size - 1], (void *) match, corr);
            num_matches[subset_size - 1]++;
        } while (!advance_nd_iterator(subset_size, subset_iter, subset_count));
    }
    
    free(adjacency);
    free(--accum_parts);
    
    /* Flatten trees while filtering already-used particles.
       Using in-order iterative btree traversal. [1] */
    
    corres_lists = (n_tupel **) calloc(frm->num_cams, sizeof(n_tupel *));
    part_used = (int **) malloc(frm->num_cams * sizeof(int *));
    for (cam = 0; cam < frm->num_cams; cam++) {
        part_used[cam] = (int *) calloc(frm->num_targets[cam], sizeof(int));
    }
    
    for (subset_size = frm->num_cams; subset_size > 1; subset_size--) {
        corres_lists[subset_size - 1] = (n_tupel *) calloc(
            num_matches[subset_size - 1] + 1, sizeof(n_tupel));
        stack = (btree_node **) realloc(stack,
            num_matches[subset_size - 1] * sizeof(btree_node *));
        flat_top = corres_lists[subset_size - 1];
        corres_lists[subset_size - 1][0].corr = -1;
        
        stack_top = 0;
        node = match_trees[subset_size - 1];
        
        /* Reverse in-order btree traversal [1] */
        while ((stack_top != 0) || (node != NULL)) {
            if (node != NULL) {
                stack[stack_top++] = node;
                node = node->right;
            } else {
                node = stack[--stack_top];
                
                /* Insert to output array only if particle wasn't taken already */
                match = (n_tupel *)(node->val);
                use_corres = 1;
                for (cam = 0; cam < frm->num_cams; cam++) {
                    if ((match->p[cam] >= 0) && part_used[cam][match->p[cam]])
                    {
                        use_corres = 0;
                        break;
                    }
                }
                
                if (use_corres) {
                    for (cam = 0; cam < 4; cam++) 
                        if (match->p[cam] >= 0)
                            part_used[cam][match->p[cam]] = 1;
                    
                    *flat_top = *((n_tupel *) (node->val));
                    flat_top++;
                }
                
                node = node->left;
            }
        }
        
        /* Release unused memory and terminate lists properly */
        btree_free(match_trees[subset_size - 1]);
        list_len = flat_top - corres_lists[subset_size - 1];
        corres_lists[subset_size - 1] = (n_tupel *) realloc(
            corres_lists[subset_size - 1], (list_len + 1) * sizeof(n_tupel));
        corres_lists[subset_size - 1][list_len].corr = -1;
    }
    
    free(match_trees);
    free(stack);
    free(num_matches);
    free(cam_set);
    free(cam_subset);
    free(subset_count);
    free(subset_iter);
    
    return corres_lists;
}

