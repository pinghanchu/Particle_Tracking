/****************************************************************************

Routine:	       	multimed.c

Author/Copyright:      	Hans-Gerd Maas

Address:	      	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
	       	       	CH - 8093 Zurich

Creation Date:          20.4.88
	
Description:			
1 point, exterior orientation, multimedia parameters   	==> X,Y=shifts
		  (special case 3 media, medium 2 = plane parallel)
							
Routines contained:		-

****************************************************************************/

#include "ptv.h"
#include <optv/parameters.h>
#include <optv/ray_tracing.h>
#include <optv/trafo.h>

double get_mmf_from_mmLUT ();
extern int lut_inited; /* temp. until optv multimed integ */

/* multimed_nlay_v2()
   
   Arguments:
   ...
   int cam - the camera number, for getting LUT data. If no LUT data is to be
        used, then this variable is ignored.
 */
void  multimed_nlay_v2 (Exterior ex, Exterior ex_o, mm_np mm, 
    double X, double Y, double Z, double *Xq, double *Yq, int cam)
{
  
  //Beat Lüthi, Nov 2007 comment actually only Xq is affected since all Y and Yq are always zero
  int		i, it=0;
  double	beta1, beta2[32], beta3, r, rbeta, rdiff, rq, mmf;
  
  // interpolation in mmLUT, if selected (requires some global variables) 
  if (lut_inited) {
      mmf = get_mmf_from_mmLUT (cam, X, Y, Z);
      
      if (mmf > 0)
	{
	  *Xq = ex.x0 + (X-ex.x0) * mmf;
	  *Yq = ex.y0 + (Y-ex.y0) * mmf;
	  return;
	}
    }
  
  // iterative procedure (if mmLUT does not exist or has no entry) 
  r = sqrt ((X-ex.x0)*(X-ex.x0)+(Y-ex.y0)*(Y-ex.y0));
  rq = r;
  
  do
    {
      beta1 = atan (rq/(ex.z0-Z));
      for (i=0; i<mm.nlay; i++)	beta2[i] = asin (sin(beta1) * mm.n1/mm.n2[i]);
      beta3 = asin (sin(beta1) * mm.n1/mm.n3);
      
      rbeta = (ex.z0-mm.d[0]) * tan(beta1) - Z * tan(beta3);
      for (i=0; i<mm.nlay; i++)	rbeta += (mm.d[i] * tan(beta2[i]));
      rdiff = r - rbeta;
      rq += rdiff;
      it++;
    }
  while (((rdiff > 0.001) || (rdiff < -0.001))  &&  it < 40);
  
  if (it >= 40)
    {
      *Xq = X; *Yq = Y;
      puts ("Multimed_nlay stopped after 40. Iteration");	return;
    }
    
  if (r != 0)
    {
      *Xq = ex.x0 + (X-ex.x0) * rq/r;
      *Yq = ex.y0 + (Y-ex.y0) * rq/r;
    }
  else
    {
      *Xq = X;
      *Yq = Y;
    }
	
}

void back_trans_Point_back(X_t,Y_t,Z_t,mm,G,cross_p,cross_c,X,Y,Z)

double X_t,Y_t,Z_t,cross_p[],cross_c[];
mm_np mm;
Glass G;
double *X,*Y,*Z;
{
    double nVe,nGl;
	nGl=sqrt(pow(G.vec_x,2.)+pow(G.vec_y,2.)+pow(G.vec_z,2.));
	nVe=sqrt( pow(cross_p[0]-(cross_c[0]-mm.d[0]*G.vec_x/nGl),2.)
		     +pow(cross_p[1]-(cross_c[1]-mm.d[0]*G.vec_y/nGl),2.)
			 +pow(cross_p[2]-(cross_c[2]-mm.d[0]*G.vec_z/nGl),2.));
	

	*X=cross_c[0]+X_t*(cross_p[0]-(cross_c[0]-mm.d[0]*G.vec_x/nGl))/nVe+Z_t*G.vec_x/nGl;
	*Y=cross_c[1]+X_t*(cross_p[1]-(cross_c[1]-mm.d[0]*G.vec_y/nGl))/nVe+Z_t*G.vec_y/nGl;
	*Z=cross_c[2]+X_t*(cross_p[2]-(cross_c[2]-mm.d[0]*G.vec_z/nGl))/nVe+Z_t*G.vec_z/nGl;
}



void trans_Cam_Point_back(ex,mm,gl,X,Y,Z,ex_t,X_t,Y_t,Z_t,cross_p,cross_c)

Exterior	ex;
mm_np		mm;
Glass		gl;
double		X, Y, Z;
Exterior	*ex_t;
double		*X_t, *Y_t, *Z_t;
double      *cross_p, *cross_c;
{
  //--Beat Lüthi June 07: I change the stuff to a system perpendicular to the interface
  double dummy;
  double nGl;
  
  nGl=sqrt(gl.vec_x*gl.vec_x+gl.vec_y*gl.vec_y+gl.vec_z*gl.vec_z);
  
  dummy=ex.x0*gl.vec_x/nGl+ex.y0*gl.vec_y/nGl+ex.z0*gl.vec_z/nGl-nGl-mm.d[0];
  cross_c[0]=ex.x0-dummy*gl.vec_x/nGl;
  cross_c[1]=ex.y0-dummy*gl.vec_y/nGl;
  cross_c[2]=ex.z0-dummy*gl.vec_z/nGl;
  ex_t->x0=0.;
  ex_t->y0=0.;
  ex_t->z0=dummy;

  dummy=nGl-X*gl.vec_x/nGl+Y*gl.vec_y/nGl+Z*gl.vec_z/nGl;
  cross_p[0]=X+dummy*gl.vec_x/nGl;
  cross_p[1]=Y+dummy*gl.vec_y/nGl;
  cross_p[2]=Z+dummy*gl.vec_z/nGl;
  *Z_t=-dummy-mm.d[0];
  dummy=sqrt( pow(cross_p[0]-(cross_c[0]-mm.d[0]*gl.vec_x/nGl),2.)
		     +pow(cross_p[1]-(cross_c[1]-mm.d[0]*gl.vec_y/nGl),2.)
			 +pow(cross_p[2]-(cross_c[2]-mm.d[0]*gl.vec_z/nGl),2.));
  *X_t=dummy;
  *Y_t=0;
      
}

void trans_Cam_Point(ex,mm,gl,X,Y,Z,ex_t,X_t,Y_t,Z_t,cross_p,cross_c)

Exterior	ex;
mm_np		mm;
Glass		gl;
double		X, Y, Z;
Exterior	*ex_t;
double		*X_t, *Y_t, *Z_t;
double      *cross_p, *cross_c;
{
  //--Beat Lüthi June 07: I change the stuff to a system perpendicular to the interface
  double dist_cam_glas,dist_point_glas,dist_o_glas; //glas inside at water 
  
  // Origin to glass.
  dist_o_glas=sqrt(gl.vec_x*gl.vec_x+gl.vec_y*gl.vec_y+gl.vec_z*gl.vec_z);
  // camera to glass. Assumes origin and camera are on opposite sides of glass.
  dist_cam_glas   = ex.x0*gl.vec_x/dist_o_glas+ex.y0*gl.vec_y/dist_o_glas+ex.z0*gl.vec_z/dist_o_glas-dist_o_glas-mm.d[0];
  dist_point_glas = X    *gl.vec_x/dist_o_glas+Y    *gl.vec_y/dist_o_glas+Z    *gl.vec_z/dist_o_glas-dist_o_glas; 

  cross_c[0]=ex.x0-dist_cam_glas*gl.vec_x/dist_o_glas;
  cross_c[1]=ex.y0-dist_cam_glas*gl.vec_y/dist_o_glas;
  cross_c[2]=ex.z0-dist_cam_glas*gl.vec_z/dist_o_glas;
  cross_p[0]=X    -dist_point_glas*gl.vec_x/dist_o_glas;
  cross_p[1]=Y    -dist_point_glas*gl.vec_y/dist_o_glas;
  cross_p[2]=Z    -dist_point_glas*gl.vec_z/dist_o_glas;

  ex_t->x0=0.;
  ex_t->y0=0.;
  ex_t->z0=dist_cam_glas+mm.d[0];

  *X_t=sqrt( pow(cross_p[0]-(cross_c[0]-mm.d[0]*gl.vec_x/dist_o_glas),2.)
		    +pow(cross_p[1]-(cross_c[1]-mm.d[0]*gl.vec_y/dist_o_glas),2.)
			+pow(cross_p[2]-(cross_c[2]-mm.d[0]*gl.vec_z/dist_o_glas),2.));
  *Y_t=0;
  *Z_t=dist_point_glas;
      
}

void back_trans_Point(X_t,Y_t,Z_t,mm,G,cross_p,cross_c,X,Y,Z)

double X_t,Y_t,Z_t,cross_p[],cross_c[];
mm_np mm;
Glass G;
double *X,*Y,*Z;
{
    
    double nVe,nGl;
	nGl=sqrt(pow(G.vec_x,2.)+pow(G.vec_y,2.)+pow(G.vec_z,2.));
	nVe=sqrt( pow(cross_p[0]-(cross_c[0]-mm.d[0]*G.vec_x/nGl),2.)
		     +pow(cross_p[1]-(cross_c[1]-mm.d[0]*G.vec_y/nGl),2.)
			 +pow(cross_p[2]-(cross_c[2]-mm.d[0]*G.vec_z/nGl),2.));

	*X = cross_c[0]-mm.d[0]*G.vec_x/nGl + Z_t*G.vec_x/nGl;
	*Y = cross_c[1]-mm.d[0]*G.vec_y/nGl + Z_t*G.vec_y/nGl;
	*Z = cross_c[2]-mm.d[0]*G.vec_z/nGl + Z_t*G.vec_z/nGl;

    if (nVe > 0) {  
        // We need this for when the cam-point line is exactly perp. to glass.
        // The degenerate case will have nVe == 0 and produce NaNs on the
        // following calculations.
        *X += X_t*(cross_p[0] - (cross_c[0] - mm.d[0]*G.vec_x/nGl))/nVe;
        *Y += X_t*(cross_p[1] - (cross_c[1] - mm.d[0]*G.vec_y/nGl))/nVe;
        *Z += X_t*(cross_p[2] - (cross_c[2] - mm.d[0]*G.vec_z/nGl))/nVe;
    }
}

/* multimed_r_nlay_v2() alculates and returns the radial shift.
   
   Arguments:
   ...
   int cam - the camera number, for getting LUT data. If no LUT data is to be
        used, then this variable is ignored.
 */
double multimed_r_nlay_v2 (Exterior ex, Exterior ex_o, mm_np mm,
    double X, double Y, double Z, int cam)
{
  int  	i, it=0;
  double beta1, beta2[32], beta3, r, rbeta, rdiff, rq, mmf;
  double ocf=1.0; // over compensation factor for faster convergence 

  double absError=100;
  int counter=0;
  double dir_water_x=-X;
  double dir_water_z=ex.z0-Z;
  double dist=pow(dir_water_x*dir_water_x+dir_water_z*dir_water_z,0.5);
  double xInInterFace,comp_parallel,comp_perpendicular,dir_air_x,dir_air_z,error_x,error_z;
  dir_water_x=dir_water_x/dist;
  dir_water_z=dir_water_z/dist;
  
  
  // 1-medium case 
  if (mm.n1==1 && mm.nlay == 1 && mm.n2[0]==1 && mm.n3==1) return (1.0);
  
  
  // interpolation in mmLUT, if selected (requires some global variables) 
  if (lut_inited) {
    mmf = get_mmf_from_mmLUT(cam, X, Y, Z);
    if (mmf > 0) return (mmf);
  }
 
  // iterative procedure 
  r = sqrt ((X-ex.x0)*(X-ex.x0)+(Y-ex.y0)*(Y-ex.y0));
  rq = r;
  
  do
    {
      beta1 = atan (rq/(ex.z0-Z));
      for (i=0; i<mm.nlay; i++)	beta2[i] = asin (sin(beta1) * mm.n1/mm.n2[i]);
      beta3 = asin (sin(beta1) * mm.n1/mm.n3);
      
      rbeta = (ex.z0-mm.d[0]) * tan(beta1) - Z * tan(beta3);
      for (i=0; i<mm.nlay; i++)	rbeta += (mm.d[i] * tan(beta2[i]));
      rdiff = r - rbeta;
      rdiff *= ocf;
      rq += rdiff;
      it++;
    }
  while (((rdiff > 0.001) || (rdiff < -0.001))  &&  it < 40);
  
  if (it >= 40)
    {
      puts ("Multimed_r_nlay_v2 stopped after 40. Iteration");	return (1.0);
    }
  
  if (r != 0)	return (rq/r);	else return (1.0);
}



void init_mmLUT (int i_cam, Calibration *cal, control_par *cpar)
{
  register int	i,j, nr, nz;
  double       	X[3], R, X1[3], Zmin, Rmax=0,Zmax, a[3];
  double       	x,y, *Ri,*Zi;
  double       	rw = 2; //was 2, has unit [mm]??, Beat Lüthi Aug 1, 2007, is ok
  Exterior	    Ex, Ex_t;
  Interior I;
  Glass G;
  ap_52 ap;
  double       	X_t,Y_t,Z_t,cross_p[3],cross_c[3],Zmin_t,Zmax_t;
  
  ap = cal->added_par;
  Ex = cal->ext_par;
  I = cal->int_par;
  G = cal->glass_par;
  
  volume_par *vpar;
    
  /* find extrema of imaged object volume */
  /* ==================================== */
  
  vpar = read_volume_par("parameters/criteria.par");
  Zmin = (vpar->Zmin_lay[1] < vpar->Zmin_lay[0]) ? vpar->Zmin_lay[1] : vpar->Zmin_lay[0];
  Zmax = (vpar->Zmax_lay[1] < vpar->Zmax_lay[0]) ? vpar->Zmax_lay[1] : vpar->Zmax_lay[0];
  
  Zmin -= fmod (Zmin, rw);
  Zmax += (rw - fmod (Zmax, rw));
  Zmin_t=Zmin;
  Zmax_t=Zmax;
  
  /* intersect with image vertices rays */
  // this is to cut down duplicate code without refactoring the function 
  //totally which would happen later.
#define UPDATE_R(p2mx, p2my) \
  pixel_to_metric (&x, &y, (p2mx), (p2my), cpar); \
  x = x - I.xh; \
  y = y - I.yh; \
  correct_brown_affin (x, y, ap, &x,&y); \
  ray_tracing(x,y, cal, *(cpar->mm), X1, a); \
  X[2] = Zmin; \
  X[1] = X1[1] + (X[2] - X1[2]) * a[1]/a[2]; \
  X[0] = X1[0] + (X[2] - X1[2]) * a[0]/a[2]; \
  \
  trans_Cam_Point(Ex, *(cpar->mm), G, X[0], X[1], X[2], &Ex_t, &X_t, &Y_t, &Z_t,\
        &cross_p,&cross_c);\
  if(Z_t<Zmin_t)Zmin_t=Z_t;\
  if(Z_t>Zmax_t)Zmax_t=Z_t;\
  \
  R = sqrt (  (X_t - Ex_t.x0)*(X_t - Ex_t.x0)\
	      + (Y_t - Ex_t.y0)*(Y_t - Ex_t.y0));\
  \
  if (R > Rmax)	Rmax = R;\
  X[2] = Zmax;\
  X[1] = X1[1] + (X[2] - X1[2]) * a[1]/a[2];\
  X[0] = X1[0] + (X[2] - X1[2]) * a[0]/a[2];\
  \
  trans_Cam_Point(Ex, *(cpar->mm), G, X[0], X[1], X[2], &Ex_t, &X_t, &Y_t, &Z_t,\
        &cross_p, &cross_c);\
  if(Z_t<Zmin_t)Zmin_t=Z_t;\
  if(Z_t>Zmax_t)Zmax_t=Z_t;\
  \
  R = sqrt (  (X_t - Ex_t.x0)*(X_t - Ex_t.x0)\
	      + (Y_t - Ex_t.y0)*(Y_t - Ex_t.y0));
  
  //--00
  UPDATE_R(0., 0.)

  //--0y
  if (R > Rmax)	Rmax = R;
  UPDATE_R(0., (double) cpar->imy)
  
  //--x0
  if (R > Rmax)	Rmax = R;  
  UPDATE_R((double) cpar->imx, 0.)

  //--xy
  if (R > Rmax)	Rmax = R;
  UPDATE_R((double) cpar->imx, (double) cpar->imy)
  
  if (R > Rmax)	Rmax = R;
  
  /* round values (-> enlarge) */
  Rmax += (rw - fmod (Rmax, rw));
    
  /* get # of rasterlines in r,z */
  
  nr = (int)(Rmax/rw + 1);
  nz = (int)((Zmax_t-Zmin_t)/rw + 1);
 
  /* create twodimensional mmLUT structure */
  trans_Cam_Point(Ex, *(cpar->mm), G, X[0], X[1], X[2], &Ex_t, &X_t, &Y_t, &Z_t,\
        &cross_p, &cross_c);\

  mmLUT[i_cam].origin.x = Ex_t.x0;
  mmLUT[i_cam].origin.y = Ex_t.y0;
  mmLUT[i_cam].origin.z = Zmin_t;
  mmLUT[i_cam].nr = nr;
  mmLUT[i_cam].nz = nz;
  mmLUT[i_cam].rw = rw;
  mmLUT[i_cam].data = (double *) malloc (nr*nz * sizeof (double));
   
  /* fill mmLUT structure */
  /* ==================== */
  
  Ri = (double *) malloc (nr * sizeof (double));
  for (i=0; i<nr; i++)	Ri[i] = i*rw;
  Zi = (double *) malloc (nz * sizeof (double));
  for (i=0; i<nz; i++)	Zi[i] = Zmin_t + i*rw;
  
  for (i=0; i<nr; i++)	for (j=0; j<nz; j++) {
    trans_Cam_Point(Ex, *(cpar->mm), G, X[0], X[1], X[2], &Ex_t, &X_t, &Y_t, &Z_t,\
        &cross_p, &cross_c);\
    mmLUT[i_cam].data[i*nz + j] = multimed_r_nlay_v2 (Ex_t, Ex, *(cpar->mm), 
        Ri[i] + Ex_t.x0, Ex_t.y0, Zi[j], i_cam);
  }
}



double get_mmf_from_mmLUT (i_cam, X,Y,Z)
int		i_cam;
double	X,Y,Z;
{
  int		i, ir,iz, nr,nz, v4[4];
  double	R, sr,sz, rw, mmf=1;
  
  rw =  mmLUT[i_cam].rw;
  
  Z -= mmLUT[i_cam].origin.z; sz = Z/rw; iz = (int) sz;	sz -= iz;
  X -= mmLUT[i_cam].origin.x;
  Y -= mmLUT[i_cam].origin.y;
  R = sqrt (X*X + Y*Y);	sr = R/rw; ir = (int) sr; sr -= ir;
    
  nz =  mmLUT[i_cam].nz;
  nr =  mmLUT[i_cam].nr;
    
  /* check whether point is inside camera's object volume */
  if (ir > nr)				return (0);
  if (iz < 0  ||  iz > nz)	return (0);
  
  /* bilinear interpolation in r/z box */
  /* ================================= */
  
  /* get vertices of box */
  v4[0] = ir*nz + iz;
  v4[1] = ir*nz + (iz+1);
  v4[2] = (ir+1)*nz + iz;
  v4[3] = (ir+1)*nz + (iz+1);
  
  /* 2. check wther point is inside camera's object volume */
  /* important for epipolar line computation */
  for (i=0; i<4; i++)
    if (v4[i] < 0  ||  v4[i] > nr*nz)	return (0);
  
  /* interpolate */
  mmf = mmLUT[i_cam].data[v4[0]] * (1-sr)*(1-sz)
    + mmLUT[i_cam].data[v4[1]] * (1-sr)*sz
    + mmLUT[i_cam].data[v4[2]] * sr*(1-sz)
    + mmLUT[i_cam].data[v4[3]] * sr*sz;
  
  return (mmf);
}



void volumedimension (Calibration cal[4], double *xmax, double *xmin, 
    double *ymax, double *ymin, double *zmax, double *zmin, control_par *cpar)
{
  int	i_cam;
  Exterior Ex;
  double X[3], R, X1[3], Zmin, Rmax = 0, Zmax, a[3];
  double x,y;
  
  volume_par *vpar;
    
  /* find extrema of imaged object volume */
  /* ==================================== */
  
  
  vpar = read_volume_par("parameters/criteria.par");
  Zmin = (vpar->Zmin_lay[1] < vpar->Zmin_lay[0]) ? vpar->Zmin_lay[1] : vpar->Zmin_lay[0];
  Zmax = (vpar->Zmax_lay[1] < vpar->Zmax_lay[0]) ? vpar->Zmax_lay[1] : vpar->Zmax_lay[0];

  *zmin=Zmin;
  *zmax=Zmax;

  // this is to cut down duplicate code without refactoring the function 
  //totally which would happen later.
#define UPDATE_R_VD(p2mx, p2my) \
  pixel_to_metric (&x, &y, (p2mx), (p2my), cpar); \
  x = x - cal[i_cam].int_par.xh; \
  y = y - cal[i_cam].int_par.yh; \
  correct_brown_affin (x, y, cal[i_cam].added_par, &x,&y); \
  ray_tracing(x,y, cal + i_cam, *(cpar->mm), X1, a); \
  X[2] = Zmin; \
  X[1] = X1[1] + (X[2] - X1[2]) * a[1]/a[2]; \
  X[0] = X1[0] + (X[2] - X1[2]) * a[0]/a[2]; \
  \
  R = sqrt (  (X[0] - Ex.x0)*(X[0] - Ex.x0)\
	      + (X[1] - Ex.y0)*(X[1] - Ex.y0));\
  \
  if ( X[0] > *xmax) *xmax=X[0];\
  if ( X[0] < *xmin) *xmin=X[0];\
  if ( X[1] > *ymax) *ymax=X[1];\
  if ( X[1] < *ymin) *ymin=X[1];\
  \
  if (R > Rmax)	Rmax = R;\
  X[2] = Zmax;\
  X[1] = X1[1] + (X[2] - X1[2]) * a[1]/a[2];\
  X[0] = X1[0] + (X[2] - X1[2]) * a[0]/a[2];\
  \
  R = sqrt (  (X[0] - Ex.x0)*(X[0] - Ex.x0)\
	      + (X[1] - Ex.y0)*(X[1] - Ex.y0));\
  if ( X[0] > *xmax) *xmax=X[0];\
  if ( X[0] < *xmin) *xmin=X[0];\
  if ( X[1] > *ymax) *ymax=X[1];\
  if ( X[1] < *ymin) *ymin=X[1];\

  for (i_cam = 0; i_cam < cpar->num_cams; i_cam++) {  
    Ex = cal[i_cam].ext_par;
    
    //--00
    UPDATE_R_VD(0., 0.)

    //--0y
    if (R > Rmax)	Rmax = R;
    UPDATE_R_VD(0., (double) cpar->imy)
  
    //--x0
    if (R > Rmax)	Rmax = R;  
    UPDATE_R_VD((double) cpar->imx, 0.)

    //--xy
    if (R > Rmax)	Rmax = R;
    UPDATE_R_VD((double) cpar->imx, (double) cpar->imy)
  }
}
