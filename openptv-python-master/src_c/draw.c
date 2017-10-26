/****************************************************************************

Author/Copyright:      	Jochen Willneff

Address:	      	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	end of 97
	
Description:	       	different drawing function to interact
                        with Tcl/Tk
	
Routines contained:     drawcross, drawvector, draw_pnr, mark_detections
                        mark_correspondences, mark_corr, mark_track_c

****************************************************************************/
#include "ptv.h"
#include "imgcoord.h"
#include <optv/parameters.h>

/* draws crosses for detected points in a displayed image 

   Arguments:
   int i - frame number
   control_par *cpar - scene parameters
*/
int trajectories_c(int i, control_par *cpar) 
{
  int   k, intx1, inty1, intx2, inty2;
  int  anz1, anz2, m, j;
  FILE *fp1;
  char val[256];
  vector *line1, *line2;
  double color;
  coord_2d p1[4], p2[4];
  sequence_par *seq_par;

      if (i < 10)             sprintf (val, "res/ptv_is.%1d", i);
      else if (i < 100)       sprintf (val, "res/ptv_is.%2d",  i);
      else       sprintf (val, "res/ptv_is.%3d",  i);
 
      fp1 = fopen (val, "r");
      
    if (cpar->num_cams < 3){
        seq_par = read_sequence_par("parameters/sequence.par", 4);
    } else {
        seq_par = read_sequence_par("parameters/sequence.par", cpar->num_cams);
    }
  color = ((double)(i - seq_par->first)) / ((double)(seq_par->last - 2 - seq_par->first));
      fscanf (fp1,"%d\n", &anz1);
      
      line1 = (vector *) malloc (anz1 * sizeof (vector));
      for (j=0;j<anz1;j++) {
	fscanf (fp1, "%d\n", &line1[j].p);
	fscanf (fp1, "%d\n", &line1[j].n);
	fscanf (fp1, "%lf\n", &line1[j].x1);
	fscanf (fp1, "%lf\n", &line1[j].y1);
	fscanf (fp1, "%lf\n", &line1[j].z1);
      }

      strcpy(val, "");     
      fclose (fp1);

      /* read next time step */     
      if (i+1 < 10)             sprintf (val, "res/ptv_is.%1d", i+1);
      else if (i+1 < 100)       sprintf (val, "res/ptv_is.%2d",  i+1);
      else       sprintf (val, "res/ptv_is.%3d",  i+1);
      
      fp1 = fopen (val, "r");      
      fscanf (fp1,"%d\n", &anz2);
      line2 = (vector *) calloc (anz2, sizeof (vector));
      
      for (j=0;j<anz2;j++) {
	fscanf (fp1, "%d\n", &line2[j].p);
	fscanf (fp1, "%d\n", &line2[j].n);
	fscanf (fp1, "%lf\n", &line2[j].x1);
	fscanf (fp1, "%lf\n", &line2[j].y1);
	fscanf (fp1, "%lf\n", &line2[j].z1);
      }
      fclose (fp1);
      m1_tr=0;
      for(j=0;j<anz1;j++) { 	
	m = line1[j].n;

	if (m >= 0)  {	  
        
	  for (k=0; k < cpar->num_cams; k++)
	    {
	      img_coord (k, line1[j].x1, line1[j].y1, line1[j].z1, 
            Ex[k],I[k], G[k], ap[k], *(cpar->mm), &p1[k].x, &p1[k].y);
	      metric_to_pixel (&p1[k].x, &p1[k].y, p1[k].x, p1[k].y, cpar);
	      
	      img_coord (k, line2[m].x1, line2[m].y1, line2[m].z1, 
            Ex[k],I[k], G[k], ap[k], *(cpar->mm), &p2[k].x, &p2[k].y);
	      metric_to_pixel(&p2[k].x, &p2[k].y, p2[k].x, p2[k].y, cpar); 
	      
	      if ( fabs( p2[k].x-zoom_x[k]) < cpar->imx/(2*zoom_f[k])
		   && ( fabs(p2[k].y-zoom_y[k]) < cpar->imy/(2*zoom_f[k])) )
		{	
        
		  intx1 = (int)(cpar->imx/2 + zoom_f[k]*(p1[k].x - zoom_x[k]));
		  inty1 = (int)(cpar->imy/2 + zoom_f[k]*(p1[k].y - zoom_y[k]));
		  intx2 = (int)(cpar->imx/2 + zoom_f[k]*(p2[k].x - zoom_x[k]));
		  inty2 = (int)(cpar->imy/2 + zoom_f[k]*(p2[k].y - zoom_y[k]));
            
        intx1_tr[k][m1_tr]=intx1;
        inty1_tr[k][m1_tr]=inty1;
        intx2_tr[k][m1_tr]=intx2;
        inty2_tr[k][m1_tr]=inty2;
		}	
       
	    }
         m1_tr++;
	}
	}


      strcpy(val, "");
      free(line1); free(line2);
  
  return 0;
}
