/****************************************************************************

Routine:	       	sortgrid.c

Author/Copyright:      	Hans-Gerd Maas

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	22.6.88

Description:	       	reads objects, detected by detection etc.,
		       	sorts them with respect to the 13 wires of the grid,
		       	detects missing points
		       	and writes the result to  a new file
			
		       	does not work in each imaginable case !
****************************************************************************/

#include <stdio.h>

#include "ptv.h"
#include "tools.h"

#include <optv/parameters.h>
#include <optv/trafo.h>

void just_plot (Exterior Ex, Interior I, Glass G, ap_52 ap, 
    int nfix, coord_3d fix[], int n_img, control_par *cpar) {
  int	       	i;
  double       	xp, yp;
  
  ncal_points[n_img]=nfix;
  
  /* reproject all calibration plate points into pixel space
     and search a detected target nearby */
  
  for (i=0; i<nfix; i++) {
      img_coord (n_img, fix[i].x, fix[i].y, fix[i].z,  Ex, I, G, ap, 
          *(cpar->mm), &xp, &yp);
      metric_to_pixel (&xp, &yp, xp, yp, cpar);
      
      /* draw projected points for check purpuses */
      x_calib[n_img][i] = (int)xp;
      y_calib[n_img][i] = (int)yp;
  }
}


int nearest_neighbour_pix (pix, num, x, y, eps)
target 	pix[];
int    	num;
double 	x, y, eps;
{
  register int	j;
  int	       	pnr = -999;
  double       	d, dmin=1e20, xmin, xmax, ymin, ymax;

  xmin = x - eps;  xmax = x + eps;  ymin = y - eps;  ymax = y + eps;

  for (j=0; j<num; j++)		    			/* candidate search */
    {
      if (pix[j].y>ymin && pix[j].y<ymax && pix[j].x>xmin && pix[j].x<xmax)
	{
	  d = sqrt ((x-pix[j].x)*(x-pix[j].x) + (y-pix[j].y)*(y-pix[j].y));
	  if (d < dmin)
	    {
	      dmin = d; pnr = j;
	    }
	}
    }
  return (pnr);
}


void sortgrid_man (Exterior Ex, Interior I, Glass G, ap_52 ap, 
    int nfix, coord_3d fix[], int num, target pix[], int n_img, 
    control_par *cpar) {
  int	       	i, j;
  double       	xp, yp, eps=10.0;
  target       	old[1024];
  FILE *fpp;
  
  /* copy and re-initialize pixel data before sorting */
  for (i=0; i<num; i++)	old[i] = pix[i];
  for (i=0; i<nfix; i++)
    {
      pix[i].pnr = -999;  pix[i].x = -999;  pix[i].y = -999;
      pix[i].n = 0; pix[i].nx = 0; pix[i].ny = 0;
      pix[i].sumg = 0;
    }
  
  
  fpp = fopen ("parameters/sortgrid.par", "r");
  if (fpp) {
    fscanf (fpp, "%lf", &eps);
    printf ("Sortgrid search radius: %.1f pixel (from sortgrid.par)\n",eps);
    fclose (fpp);
  }
  else {
    printf ("parameters/sortgrid.par does not exist, ");
    printf ("using default search radius 10 pixel\n");
  }
  
  
  /* reproject all calibration plate points into pixel space
     and search a detected target nearby */
  
  for (i=0; i<nfix; i++) {
      img_coord (n_img, fix[i].x, fix[i].y, fix[i].z,  Ex, I, G, ap, 
        *(cpar->mm), &xp, &yp);
      metric_to_pixel (&xp, &yp, xp, yp, cpar);
      
      /* draw projected points for check purpuses */
      x_calib[n_img][i] = (int)xp;
      y_calib[n_img][i] = (int)yp;

      if (xp > -eps  &&  yp > -eps  &&  xp < cpar->imx + eps  &&  yp < cpar->imy + eps) {
	  j = nearest_neighbour_pix (old, num, xp, yp, eps);
	  
	  if (j != -999)
	    {
	      pix[i] = old[j];  pix[i].pnr = fix[i].pnr;
	      z_calib[n_img][i]=fix[i].pnr; // Alex, 18.05.11
	      printf("z_calib[%d][%d]=%d\n",n_img,i,z_calib[n_img][i]);
	    }
	}
    }
}
