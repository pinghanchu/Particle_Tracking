#include "ptv.h"
#include <optv/parameters.h>
#include <optv/trafo.h>
#include "epi.h"
#include "tools.h"

/* Variables used here and defined extern in 'globals.h' */
int rclick_intx1[4], rclick_inty1[4], rclick_intx2[4], rclick_inty2[4],\
    rclick_points_x1[4][10000], rclick_points_y1[4][10000], rclick_count[4],\
    rclick_points_intx1, rclick_points_inty1;

int mouse_proc_c (int click_x, int click_y, int kind, int num_image, volume_par *vpar, \
	control_par *cpar){
  int     i, j, n, zf;
  double  x, y;
  double  xa12, xb12, ya12, yb12;
  int     k, pt1, intx1, inty1, count, intx2, inty2, pt2;
  candidate cand[maxcand];
  
  
  printf("entered mouse_proc_c \n");
 
  if (zoom_f[0] == 1) {zf = 2;} else { zf = zoom_f[0];}

  n=num_image;
  if (examine)	zf *= 2;
  
  switch (kind) 
    {

/* -------------------------- MIDDLE MOUSE BUTTON ---------------------------------- */

   

    case 3: /* generate epipolar line segments */
      
      /* get geometric coordinates of nearest point in img[n] */
      x = (float) (click_x - cpar->imx/2)/zoom_f[n] + zoom_x[n];
      y = (float) (click_y - cpar->imy/2)/zoom_f[n] + zoom_y[n];

      pixel_to_metric (&x, &y, x,y, cpar);
      x -= I[n].xh;	y -= I[n].yh;
      correct_brown_affin (x, y, ap[n], &x, &y);
      
      k = nearest_neighbour_geo (geo[n], num[n], x, y, 0.05);
      
      if (k == -999){	  
	      printf  ("No point near click coord! Click again! \n"); 
	      return -1;
	      }
	      
	      pt1 = geo[n][k].pnr;

      intx1 = (int) ( cpar->imx/2 + zoom_f[n] * (pix[n][pt1].x-zoom_x[n]));
      inty1 = (int) ( cpar->imy/2 + zoom_f[n] * (pix[n][pt1].y-zoom_y[n]));
      rclick_points_intx1=intx1;
      rclick_points_inty1=inty1;
    
      //drawcross (interp, intx1, inty1, cr_sz+2, n, "BlueViolet");

      printf ( "pt1,nx,ny,n,sumg: %d %d %d %d %d\n", pt1, pix[n][pt1].nx, pix[n][pt1].ny,
	       pix[n][pt1].n, pix[n][pt1].sumg);  
	       	       	       
	       for (i = 0; i < cpar->num_cams; i++) if (i != n) {
		   /* calculate epipolar band in img[i] */
		   epi_mm (i, geo[n][k].x,geo[n][k].y,
			   Ex[n],I[n], G[n], Ex[i],I[i], G[i], *(cpar->mm), vpar,
			   &xa12, &ya12, &xb12, &yb12);
		   
		   /* search candidate in img[i] */
		   printf("\ncandidates in img: %d\n", i);
		   find_candidate_plus_msg (geo[i], pix[i], num[i],
					    xa12, ya12, xb12, yb12,
					    pix[n][pt1].n, pix[n][pt1].nx, pix[n][pt1].ny,
					    pix[n][pt1].sumg, cand, &count, i, vpar, cpar);

		   distort_brown_affin (xa12,ya12, ap[i], &xa12,&ya12);
		   distort_brown_affin (xb12,yb12, ap[i], &xb12,&yb12);
		   xa12 += I[i].xh;	ya12 += I[i].yh;
		   xb12 += I[i].xh;	yb12 += I[i].yh;
            
		   metric_to_pixel(&xa12, &ya12, xa12, ya12, cpar);
		   metric_to_pixel(&xb12, &yb12, xb12, yb12, cpar);
            
		   intx1 = (int) ( cpar->imx/2 + zoom_f[i] * (xa12 - zoom_x[i]));
		   inty1 = (int) ( cpar->imy/2 + zoom_f[i] * (ya12 - zoom_y[i]));
		   intx2 = (int) ( cpar->imx/2 + zoom_f[i] * (xb12 - zoom_x[i]));
		   inty2 = (int) ( cpar->imy/2 + zoom_f[i] * (yb12 - zoom_y[i]));
           
            rclick_intx1[i]=intx1;
            rclick_inty1[i]=inty1;
            rclick_intx2[i]=intx2;
            rclick_inty2[i]=inty2;


		  // drawvector ( interp, intx1, inty1, intx2, inty2, 1, i, val);
            rclick_count[i]=count;
                   for (j=0; j<count; j++)
                     {
                       pt2 = cand[j].pnr;
                       intx2 = (int) ( cpar->imx/2 + zoom_f[i] * (pix[i][pt2].x - zoom_x[i]));
                       inty2 = (int) ( cpar->imy/2 + zoom_f[i] * (pix[i][pt2].y - zoom_y[i]));
                         rclick_points_x1[i][j]=intx2;
                         rclick_points_y1[i][j]=inty2;
                       //drawcross (interp, intx2, inty2, cr_sz+2, i, "orange");
                     }
   
		   
		 }

	       break;

	       	       
    case 4: /* delete points, which should not be used for orientation */

      
      j = kill_in_list (n, num[n], click_x, click_y);
      if (j != -1)
	{
	  num[n] -= 1;
	  printf ("point %d deleted", j);  
	}
      else {
	  printf ("no point near click coord !");  
      }
      break;

    }
    printf("finished mouse_proc_c \n");
  return 0;
  
}
