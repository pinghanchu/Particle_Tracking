/****************************************************************************

Routine:	       	pointpos.c

Author/Copyright:      	Hans-Gerd Maas

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	July 1988, June 1989
	
Description:	       	point positioning with least squares adjustment
		       	(multimedia, 2 or 3 channels)
	
Routines contained:		

****************************************************************************/
#include "ptv.h"
#include <optv/parameters.h>
#include <optv/ray_tracing.h>

void dist_to_ray(x, y, Ex, I, G, ap, mm, Xp,Yp,Zp, dist)

Exterior	Ex;
Interior	I;
Glass   	G;
ap_52		ap;
mm_np		mm;
double		x, y,Xp,Yp,Zp, *dist;
{
  double    gX[4],gY[4],gZ[4],a[4],b[4],c[4];
  double    x01,x02,x03,x12,x13,x23;
  double    y01,y02,y03,y12,y13,y23;
  double    z01,z02,z03,z12,z13,z23;
  
    
  	
  *dist=1;

}

void det_lsq_3d (Calibration *cal, mm_np mm, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double *Xp, double *Yp, double *Zp, int num_cams) {
	    int     i,count_inner=0,n,m, flag[4] = {0., 0., 0., 0.};
	    double  d_inner=0.,x,y;
	    double X[4][3], a[4][3];
        double dist,dist_error,X_pos[6],Y_pos[6],Z_pos[6],XX,YY,ZZ,si0,sqX,sqY,sqZ;
	    
		if(x1>-999){
			flag[0]=1;
			x = x1 - I[0].xh;
	        y = y1 - I[0].yh;
	        //correct_brown_affin (x, y, ap[0], &x, &y);
		    ray_tracing(x,y, &(cal[0]), mm, X[0], a[0]);
		}		
		if(x2>-999){
			flag[1]=1;
			x = x2 - I[1].xh;
	        y = y2 - I[1].yh;
	        //correct_brown_affin (x, y, ap[1], &x, &y);
		    ray_tracing(x,y, &(cal[1]), mm, X[1], a[1]);
		}		
		if(x3>-999){
			flag[2]=1;
			x = x3 - I[2].xh;
	        y = y3 - I[2].yh;
	        //correct_brown_affin (x, y, ap[2], &x, &y);
		    ray_tracing(x,y, &(cal[2]), mm, X[2], a[2]);
		}		
		if(x4>-999){
			flag[3]=1;
			x = x4 - I[3].xh;
	        y = y4 - I[3].yh;
	        //correct_brown_affin (x, y, ap[3], &x, &y);
		    ray_tracing(x,y, &(cal[3]), mm, X[3], a[3]);
		}

		count_inner=0;
		for (n = 0; n < num_cams; n++){
			for(m = n+1; m < num_cams; m++){
				if(flag[n]==1 && flag[m]==1){
                    mid_point(X[n][0], X[n][1], X[n][2],
                        a[n][0], a[n][1], a[n][2], 
                        X[m][0], X[m][1], X[m][2],
                        a[m][0], a[m][1], a[m][2], &dist,&XX,&YY,&ZZ);
                    d_inner += dist;
					X_pos[count_inner]=XX;Y_pos[count_inner]=YY;Z_pos[count_inner]=ZZ;
					count_inner++;
				}
			}
		}
        d_inner/=(double)count_inner;		
		XX=0.;YY=0.;ZZ=0.;
		for(i=0;i<count_inner;i++){
           XX+=X_pos[i]; 
		   YY+=Y_pos[i];
		   ZZ+=Z_pos[i];
		}
		XX/=(double)count_inner;YY/=(double)count_inner;ZZ/=(double)count_inner;
		//end of new det_lsq
		*Xp=XX;
		*Yp=YY;
		*Zp=ZZ;

		//statistics
		si0=0.;sqX=0.;sqY=0.;sqZ=0.;
		for(i=0;i<count_inner;i++){
           si0+=pow(X_pos[i]-XX,2.)+pow(Y_pos[i]-YY,2.)+pow(Z_pos[i]-ZZ,2.);
           sqX+=pow(X_pos[i]-XX,2.); 
		   sqY+=pow(Y_pos[i]-YY,2.);
		   sqZ+=pow(Z_pos[i]-ZZ,2.);		   
		}
		si0/=(double)count_inner;sqX/=(double)count_inner;sqY/=(double)count_inner;sqZ/=(double)count_inner;
		
		mean_sigma0 += pow(si0,0.5);
        rmsX += pow(sqX,0.5);
        rmsY += pow(sqY,0.5);
        rmsZ += pow(sqZ,0.5);
		//end of statistics

}

