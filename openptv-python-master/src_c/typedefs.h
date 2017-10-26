#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <optv/tracking_frame_buf.h>
#include <optv/calibration.h>

#define sqr(x) x*x
#define maxcand 200 //Beat changed it on 090325
#define maxtrack 64 //Beat changed it on 090325

#define M 20000 //Beat changed it on 090325

#ifndef  M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

typedef struct
{
  int pnr;
  double x, y, z;
    }
coord_3d;

typedef struct
{
  int pnr;
  double x, y;
}
coord_2d;

typedef struct
{
  double n1,n2,n3,d;
  int    lut;
}
mm_3p;

typedef struct
{
  int 	pos, status;
  short	xmin, xmax, ymin, ymax;
  int   n, sumg;
  double  x, y;
  int   unr, touch[4], n_touch;	/* unified with target unr, touching ... */
}
peak;

typedef struct
{
  short	       	x,y;
  unsigned char	g;
  short	       	tnr;
}
targpix;

typedef struct
{
  int  	left, right;
  double  tol, corr;
}
conjugate;

typedef struct
{
  int     p[4];
  double  corr;
}
n_tupel;

typedef struct
{
  int    	p1;	       	/* point number of master point */
  int    	n;	       	/* # of candidates */
  int    	p2[maxcand];	/* point numbers of candidates */
  double	corr[maxcand];	/* feature based correlation coefficient */
  double	dist[maxcand];	/* distance perpendicular to epipolar line */
}
correspond;	       	/* correspondence candidates */

typedef struct
{
  coord_3d	origin;
  int          	nr,nz;
  double       	rw;
  double       	*data;
}
mm_LUT;


typedef struct /* struct for what was found to corres */
{
 int ftnr, freq, whichcam[4];
}
foundpix;

typedef struct
{
  int multi, h[maxtrack], freq[maxtrack];
  double quali[maxtrack];
}
currlist;

typedef struct
{
  int p,n;
  double x1, y1, z1;
  int type;
}
vector;

typedef struct
{
  int z1, c[maxcand], n[maxcand];
  double quali[maxcand];
}
prevlist;


#define MAX_FILENAME_LEN 1024
#define FILENAME_IN "res/rt_is"
#define FILENAME_OUT "res/ptv_is"
#define COORDPARAFILE "res/coord.par"
#define STATUSFILE "res/track.out"

#endif
