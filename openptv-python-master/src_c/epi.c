#include <math.h>
#include <stdio.h>
#include <string.h>

#include "epi.h"
#include "globals.h"

#include <optv/lsqadj.h>
#include <optv/ray_tracing.h>

int dumbbell_pyptv;

double epi_line (xl, yl, Ex1, I1, G1, Ex2, I2, G2)

double    xl, yl;
Exterior  Ex1, Ex2;
Interior  I1, I2;
Glass     G1, G2;

{
  int i,j;
  double m2;
  double vect1[3], vect2[3], vect3[3], nk[3], n2[3],
    p1l[3], K2[3], k2[3], D2t[3][3];
  void crossprod ();

  /* base O1 -> O2 */
  vect1[0] = Ex2.x0 - Ex1.x0;
  vect1[1] = Ex2.y0 - Ex1.y0;
  vect1[2] = Ex2.z0 - Ex1.z0;

  /* coordinates of arbitrary point P1 in image1 */
  p1l[0] = xl;  p1l[1] = yl;	p1l[2] = - I1.cc;

  /* beam O1 -> P in space */
  matmul(vect2, (double *) Ex1.dm, p1l, 3, 3, 1, 3, 3);

  /* normale to epipolar plane */
  crossprod (vect1,vect2,nk);

  /* normale to image2 */
  vect3[0] = 0;	vect3[1] = 0;	vect3[2] = - I2.cc;

  /* normale to image 2, in space */
  matmul(n2, (double *) Ex2.dm, vect3, 3, 3, 1, 3, 3);

  /* epipolar line in image2, in space */
  crossprod (nk,n2,K2);

  /* epipolar line in image2 */
  for (i=0; i<3; i++)  for (j=0; j<3; j++)  D2t[i][j] = Ex2.dm[j][i];
  matmul(k2, (double *) D2t, K2, 3, 3, 1, 3, 3);
  m2 = k2[1] / k2[0];
  return (m2);
}


/* epi_mm() finds epipolar line search window.
   
   Arguments:
   int cam - number of the camera on which the epipolar line is drawn. Only
      needed if the LUT structure is used.
   double x1, y1 - image-plane coordinate of point in 1st image.
   Exterior Ex1, Interior I1, Glass G1 - calibration of the camera from which
      the line emanates
   Exterior Ex2, Interior I2, Glass G2 - calibration of the camera on which 
      the epipolar line is drawn.
   mm_np mmp - multimed param. (layers)
   volume_par *vpar - parameters of the observed volume
   
   Output parameters:
   double *xmin, *ymin, *xmax, *ymax - output search window.
*/
int epi_mm (int cam, double x1, double y1, Exterior Ex1, Interior I1, Glass G1,
    Exterior Ex2, Interior I2, Glass G2, mm_np mmp, volume_par *vpar, 
    double *xmin, double *ymin, double *xmax, double *ymax)
{
  /*  ray tracing gives the point of exit and the direction
      cosines at the waterside of the glass;
      min. and max. depth give window in object space,
      which can be transformed into _2 image
      (use img_xy_mm because of comparison with img_geo)  */

  double a[3], xa,ya,xb,yb;
  double X1[3], X[3];
  double Zmin, Zmax;
  Calibration cal1;
  memcpy(&(cal1.ext_par), &Ex1, sizeof(Exterior));\
  memcpy(&(cal1.int_par), &I1, sizeof(Interior));\
  memcpy(&(cal1.glass_par), &G1, sizeof(Glass));\

  ray_tracing(x1,y1, &(cal1), mmp, X1, a);

  /* calculate min and max depth for position (valid only for one setup) */
  Zmin = vpar->Zmin_lay[0]
    + (X1[0] - vpar->X_lay[0]) * (vpar->Zmin_lay[1] - vpar->Zmin_lay[0]) / 
    (vpar->X_lay[1] - vpar->X_lay[0]);
  Zmax = vpar->Zmax_lay[0]
    + (X1[0] - vpar->X_lay[0]) * (vpar->Zmax_lay[1] - vpar->Zmax_lay[0]) / 
    (vpar->X_lay[1] - vpar->X_lay[0]);

  X[2] = Zmin; \
  X[1] = X1[1] + (X[2] - X1[2]) * a[1]/a[2]; \
  X[0] = X1[0] + (X[2] - X1[2]) * a[0]/a[2]; \
  img_xy_mm_geo(cam, X[0], X[1], X[2], Ex2, I2, G2, mmp, &xa, &ya);

  X[2] = Zmax; \
  X[1] = X1[1] + (X[2] - X1[2]) * a[1]/a[2]; \
  X[0] = X1[0] + (X[2] - X1[2]) * a[0]/a[2]; \
  img_xy_mm_geo(cam, X[0], X[1], X[2], Ex2, I2, G2, mmp, &xb, &yb);

  /*  ==> window given by xa,ya,xb,yb  */

  *xmin = xa;  *ymin = ya;  *xmax = xb;  *ymax = yb;

  return (0);
}

int epi_mm_2D (x1, y1, Ex1, I1, G1, mmp, vpar, xp,yp,zp)

double     x1, y1;	  	/* input coord */
Exterior   Ex1;           	/* orientation data */
Interior   I1;	      	/* orientation data */
Glass      G1;	      	/* glass data */
mm_np	   mmp;		        /* multimed param. (layers) */
volume_par *vpar;
double *xp, *yp, *zp;
{
  /*  ray tracing gives the point of exit and the direction
      cosines at the waterside of the glass;
      min. and max. depth give window in object space,
      which can be transformed into _2 image
      (use img_xy_mm because of comparison with img_geo)  */

  double a[3], xa,ya,xb,yb;
  double X1[3], X[3];
  double Zmin, Zmax;
  Calibration cal1;
  memcpy(&(cal1.ext_par), &Ex1, sizeof(Exterior));\
  memcpy(&(cal1.int_par), &I1, sizeof(Interior));\
  memcpy(&(cal1.glass_par), &G1, sizeof(Glass));\
  
  ray_tracing(x1,y1, &(cal1), mmp, X1, a);

  /* calculate min and max depth for position (valid only for one setup) */
  Zmin = vpar->Zmin_lay[0]
    + (X1[0] - vpar->X_lay[0]) * (vpar->Zmin_lay[1] - vpar->Zmin_lay[0]) / 
    (vpar->X_lay[1] - vpar->X_lay[0]);
  Zmax = vpar->Zmax_lay[0]
    + (X1[0] - vpar->X_lay[0]) * (vpar->Zmax_lay[1] - vpar->Zmax_lay[0]) /
    (vpar->X_lay[1] - vpar->X_lay[0]);

  X[2] = 0.5*(Zmin+Zmax);   
  X[0] = X1[0] + (X[2] - X1[2]) * a[0]/a[2];   
  X[1] = X1[1] + (X[2] - X1[2]) * a[1]/a[2];
  
  *xp = X[0];
  *yp = X[1];
  *zp = X[2];

  return (0);
}

void find_candidate_plus (crd, pix, num, xa,ya,xb,yb, n, nx, ny, sumg,
						  cand, count, nr, vpar, cpar)
/*  binarized search in a x-sorted coord-set, exploits shape information  */

coord_2d	crd[];
target		pix[];
int    		num, *count;
double		xa, ya, xb, yb;
int    		n, nx, ny, sumg;
candidate	cand[];
int	       	nr;	       	/* image number for ap etc. */
volume_par *vpar;
control_par *cpar;

{
  register int	j;
  int dummy;
  int	       	j0, dj, p2;
  double      	m, b, d, temp, qn, qnx, qny, qsumg, corr;
  double       	xmin, xmax, ymin, ymax,particle_size;
  int dumbbell=0;
  double tol_band_width;
  
//Beat Mai 2010 for dumbbell

    if (dumbbell_pyptv==1) dumbbell=1;

  if (dumbbell==0){
	
	/* I overwrite this change using the old version of tol_band_width to be just the 
	value in millimeters given by the user in the main parameters/observation volume */
	
	  /////here is new Beat version of April 2010
	  /*
	  if (nx>ny) particle_size=nx;
	  else       particle_size=ny;
	  tol_band_width = vpar->eps0*0.5*(pix_x + pix_y)*particle_size;
	  */
	  
	  tol_band_width = vpar->eps0;
  }
  else{
      tol_band_width = vpar->eps0;
  }
  if(tol_band_width <= 0){
       tol_band_width = 1e-6;
  }

  /* define sensor format for search interrupt */
  xmin = (-1) * cpar->pix_x * cpar->imx/2;
  xmax = cpar->pix_x * cpar->imx/2;
  
  ymin = (-1) * cpar->pix_y * cpar->imy/2;
  ymax = cpar->pix_y * cpar->imy/2;
  
  xmin -= I[nr].xh;	ymin -= I[nr].yh;
  xmax -= I[nr].xh;	ymax -= I[nr].yh;
  correct_brown_affin (xmin,ymin, ap[nr], &xmin,&ymin);
  correct_brown_affin (xmax,ymax, ap[nr], &xmax,&ymax);

  for (j=0; j<4; j++)
    {
      cand[j].pnr = -999;  cand[j].tol = -999;  cand[j].corr = -999;
    }


  /* line equation: y = m*x + b */
  if (xa == xb)	xa += 1e-10;
  m = (yb-ya)/(xb-xa);  b = ya - m*xa;

  if (xa > xb)
    {
      temp = xa;  xa = xb;  xb = temp;
    }
  if (ya > yb)
    {
      temp = ya;  ya = yb;  yb = temp;
    }

  if ( (xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax))  /* sensor area */
    {
      /* binarized search for start point of candidate search */
      for (j0=num/2, dj=num/4; dj>1; dj/=2)
	{
	  if (crd[j0].x < (xa - tol_band_width))  j0 += dj;
	  else  j0 -= dj;
	}
      j0 -= 12;  if (j0 < 0)  j0 = 0;		       	/* due to trunc */

      for (j=j0, *count=0; j<num; j++)			/* candidate search */
	{
	  if (crd[j].x > xb+tol_band_width)  return;		/* finish search */

	  if ((crd[j].y > ya-tol_band_width) && (crd[j].y < yb+tol_band_width))
	    {
	      if ((crd[j].x > xa-tol_band_width) && (crd[j].x < xb+tol_band_width))
		{
		  d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));
          
		  /* Beat: modified in April 2010 to allow for better treatment of 
		  //different sized traced particles, in particular colloids and tracers
		  if ( d < eps ){
          */
          /////here is new Beat version of April 2010
		   //if (nx>ny) particle_size=nx;
		   //else       particle_size=ny;
		   if ( d < tol_band_width ){
		   ///////end of new Beat version

		      p2 = crd[j].pnr;
		      if (n  < pix[p2].n)      	qn  = (double) n/pix[p2].n;
		      else		       	qn  = (double) pix[p2].n/n;
		      if (nx < pix[p2].nx)	qnx = (double) nx/pix[p2].nx;
		      else		       	qnx = (double) pix[p2].nx/nx;
		      if (ny < pix[p2].ny)	qny = (double) ny/pix[p2].ny;
		      else		       	qny = (double) pix[p2].ny/ny;
		      if (sumg < pix[p2].sumg)
			        qsumg = (double) sumg/pix[p2].sumg;
		      else	qsumg = (double) pix[p2].sumg/sumg;

		      // empirical correlation coefficient
			  // from shape and brightness parameters 
		      corr = (4*qsumg + 2*qn + qnx + qny);
		      // create a tendency to prefer those matches
			  // with brighter targets 
		      corr *= ((double) (sumg + pix[p2].sumg));

		      if (qn >= vpar->cn && qnx >= vpar->cnx && \
                 qny >= vpar->cny && qsumg > vpar->csumg) {

				 if ( *count < maxcand) {
			        cand[*count].pnr = j;
			        cand[*count].tol = d;
			        cand[*count].corr = corr;
			        (*count)++;
		         } else {
			        dummy=(int)maxcand;
			        printf("in find_candidate_plus: count > maxcand\n");}
			     }
		      }
		   }
           
           
	    }
	}
    }

  else  *count = -1;	 	/* out of sensor area */
}




void find_candidate_plus_msg (crd, pix, num, xa,ya,xb,yb, n, nx, ny, sumg,
							  cand, count, i12, vpar, cpar)

/*  binarized search in a x-sorted coord-set, exploits shape information  */
/*  gives messages (in examination)  */

coord_2d	crd[];
target		pix[];
int    		num, *count, i12;
double		xa, ya, xb, yb;
int    		n, nx, ny, sumg;
volume_par *vpar;
control_par *cpar;
/*
candidate	cand[3];
*/
candidate	cand[];

{
  register int	j;
  int	       	j0, dj, p2;
  double        m, b, d, temp, qn, qnx, qny, qsumg, corr;
  double       	xmin, xmax, ymin, ymax;
  double tol_band_width,particle_size;

  /* define sensor format for search interrupt */
  xmin = (-1) * cpar->pix_x * cpar->imx/2;
  xmax = cpar->pix_x * cpar->imx/2;
  
  ymin = (-1) * cpar->pix_y * cpar->imy/2;
  ymax = cpar->pix_y * cpar->imy/2;

  xmin -= I[i12].xh;	ymin -= I[i12].yh;
  xmax -= I[i12].xh;	ymax -= I[i12].yh;
  correct_brown_affin (xmin,ymin, ap[i12], &xmin,&ymin);
  correct_brown_affin (xmax,ymax, ap[i12], &xmax,&ymax);

  /* see the note in the find_candidate_plus about overwriting the tol_band_width definition */
  tol_band_width = vpar->eps0;

  for (j=0; j<4; j++)
    {
      cand[j].pnr = -999;  cand[j].tol = 999;
    }
  m = (yb-ya)/(xb-xa);  b = ya - m*xa;   /* line equation: y = m*x + b */

  if (xa > xb)
    {
      temp = xa;  xa = xb;  xb = temp;
    }
  if (ya > yb)
    {
      temp = ya;  ya = yb;  yb = temp;
    }

  if ( (xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax)) /* sensor area */
    {
      /* binarized search for start point of candidate search */
      for (j0=num/2, dj=num/4; dj>1; dj/=2)
	{
	  if (crd[j0].x < (xa - tol_band_width))  j0 += dj;
	  else  j0 -= dj;
	}
      j0 -= 12;  if (j0 < 0)  j0 = 0;  	/* due to trunc */

      for (j=j0, *count=0; j<num; j++) 	/* candidate search */
	{
	  if (crd[j].x > xb+tol_band_width)  return;      	/* finish search */

	  if ((crd[j].y > ya-tol_band_width) && (crd[j].y < yb+tol_band_width))
	    {
	      if ((crd[j].x > xa-tol_band_width) && (crd[j].x < xb+tol_band_width))
		{
		  d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));
          if ( d < tol_band_width ){
		      p2 = crd[j].pnr;
		      if (n  < pix[p2].n)      	qn  = (double) n/pix[p2].n;
		      else		       	qn  = (double) pix[p2].n/n;
		      if (nx < pix[p2].nx)	qnx = (double) nx/pix[p2].nx;
		      else		       	qnx = (double) pix[p2].nx/nx;
		      if (ny < pix[p2].ny)	qny = (double) ny/pix[p2].ny;
		      else		       	qny = (double) pix[p2].ny/ny;
		      if (sumg < pix[p2].sumg)
			qsumg = (double) sumg/pix[p2].sumg;
		      else	qsumg = (double) pix[p2].sumg/sumg;


		      /* empirical correlation coefficient
			 from shape and brightness parameters */
		      corr = (4*qsumg + 2*qn + qnx + qny);
		      /* create a tendency to prefer those matches
			 with brighter targets */
		      corr *= ((double) (sumg + pix[p2].sumg));

            if (qn >= vpar->cn && qnx >= vpar->cnx && qny >= vpar->cny && 
                qsumg > vpar->csumg)
			{
			  if (*count>=maxcand)
			    { printf("More candidates than (maxcand): %d\n",*count); return; }
			  cand[*count].pnr = p2;
			  cand[*count].tol = d;
 			  cand[*count].corr = corr;
			  (*count)++;
			  printf ("candidates: pnt: %d corr:%3.0f epi-dist:%3.1f \n", p2, corr, d*1000);
			}
		    }
		}
	    }
	}
      if (*count == 0)  puts ("- - -");
    }
  else  *count = -1;		       	       /* out of sensor area */

}





void crossprod (a,b,c)

double  a[3], b[3], c[3];

{
	c[0] = a[1] * b[2]  -  a[2] * b[1];
	c[1] = a[2] * b[0]  -  a[0] * b[2];
	c[2] = a[0] * b[1]  -  a[1] * b[0];
}
