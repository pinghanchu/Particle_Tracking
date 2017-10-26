/****************************************************************************

Author/Copyright:      	Jochen Willneff

Address:	      	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	end of 99 ...going on

Description:	       	different search routines and track tools

Routines contained:     pix_in_next, candsearch_in_pix, searchposition,
                        predict, readseqtrackcrit, readtrackdata,
                        searchquader

****************************************************************************/
#include "ptv.h"
#include "imgcoord.h"
#include <optv/parameters.h>
#include <optv/trafo.h>

int pix_in_next (next, num, x, y, dl, dr, du, dd, found)
target  next[];
int     num;
double  x, y, dl, dr, du, dd;
int found[POSI];
{
  /* search of POSI near candidates in targetlist */

  int  j, j0, dj;
  int  zaehler=0;
  double  d, xmin, xmax, ymin, ymax;
  double dcand[POSI];
  int cand[POSI];
  xmin = x - dl;  xmax = x + dr;  ymin = y - du;  ymax = y + dd;


  /* binarized search for start point of candidate search */
  for (j0=num/2, dj=num/4; dj>1; dj/=2)
    {
      if (next[j0].y < ymin) j0 += dj;
      else j0 -= dj;
    }

  j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to trunc */
  for (j=j0; j<num; j++)		       	/* candidate search */
    {
      if (next[j].y > ymax)  break;	       	/* finish search */
      if (next[j].x > xmin  &&  next[j].x < xmax && next[j].y > ymin  &&  next[j].y < ymax && next[j].tnr>0)
	{
	  d = sqrt ((x-next[j].x)*(x-next[j].x) + (y-next[j].y)*(y-next[j].y));
	  cand[zaehler]=next[j].tnr; dcand[zaehler]=d;
	  zaehler++;
	}
    }

  return (zaehler);
}


int candsearch_in_pix (target next[], int num, double x, double y, double dl,
    double dr, double du, double dd, int p[4], control_par *cpar) {
  /* search of four near candidates in targetlist */

  int  	  j, j0, dj, pnr = -999;
  int  zaehler=0, p1, p2, p3, p4;
  double  d, dmin=1e20, xmin, xmax, ymin, ymax;
  double d1, d2, d3, d4;
  xmin = x - dl;  xmax = x + dr;  ymin = y - du;  ymax = y + dd;

  if(xmin<0.0) xmin=0.0;
  if(xmax > cpar->imx)
        xmax = cpar->imx;
  if(ymin<0.0) ymin=0.0;
  if(ymax > cpar->imy)
    ymax = cpar->imy;

  if(x<0.0) x=0.0;
  if(x > cpar->imx)
    x = cpar->imx;
  if(y<0.0) y=0.0;
  if(y > cpar->imy)
    y = cpar->imy;

  p1 = p2 = p3 = p4 = -999;
  d1 = d2 = d3 = d4 = dmin;

  if (x >= 0.0 && x <= cpar->imx ) { if (y >= 0.0 && y <= cpar->imy ) {

  /* binarized search for start point of candidate search */
  for (j0=num/2, dj=num/4; dj>1; dj/=2)
    {
      if (next[j0].y < ymin) j0 += dj;
      else j0 -= dj;
    }

  j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to trunc */
  for (j=j0; j<num; j++)		       	/* candidate search */
    {
      if (next[j].tnr != -1 ) {
	if (next[j].y > ymax )  break;	       	/* finish search */
	if (next[j].x > xmin  &&  next[j].x < xmax && next[j].y > ymin  &&  next[j].y < ymax )
	  {
	    d = sqrt ((x-next[j].x)*(x-next[j].x) + (y-next[j].y)*(y-next[j].y));

	    if (d < dmin) { dmin = d; pnr = j;
	    }
	    if ( d < d1 )
	      {
		p4=p3; p3=p2; p2=p1; p1=j;
		d4=d3; d3=d2; d2=d1; d1=d;
	      }
	    else if ( d1 < d &&  d < d2 )
	      {
		p4=p3; p3=p2; p2=j;
		d4=d3; d3=d2; d2=d;
	      }
	    else if ( d2 < d && d < d3 )
	      {
		p4=p3; p3=j;
		d4=d3; d3=d;
	      }
	    else if ( d3 < d && d < d4 )
	      {
		p4=j;
		d4=d;
	      }
	  }
      }
    }

  p[0]=p1;
  p[1]=p2;
  p[2]=p3;
  p[3]=p4;
  for (j=0; j<4; j++) if ( p[j] != -999 ) {zaehler++; }
  } }
  return (zaehler);
}


int candsearch_in_pixrest(target  next[], int num, double x, double y,
    double dl, double dr, double du, double dd, int p[4], control_par *cpar) {
  /* search of four near candidates in targetlist */

  int  	  j, j0, dj;
  int  zaehler=0, p1, p2, p3, p4;
  double  d, dmin=1e20, xmin, xmax, ymin, ymax;
  xmin = x - dl;  xmax = x + dr;  ymin = y - du;  ymax = y + dd;

  if(xmin<0.0) xmin=0.0;
  if(xmax > cpar->imx) 
    xmax = cpar->imx;
  if(ymin<0.0) ymin=0.0;
  if(ymax > cpar->imy)
    ymax = cpar->imy;

  if(x<0.0) x=0.0;
  if(x > cpar->imx)
    x = cpar->imx;
  if(y<0.0) y=0.0;
  if(y > cpar->imy)
    y = cpar->imy;

  p1 = p2 = p3 = p4 = -999;

  /* binarized search for start point of candidate search */
  for (j0=num/2, dj=num/4; dj>1; dj/=2)
    {
      if (next[j0].y < ymin) j0 += dj;
      else j0 -= dj;
    }

  j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to trunc */
  for (j=j0; j<num; j++)		       	/* candidate search */
    {
      if (next[j].tnr == -1 ) {
	if (next[j].y > ymax )  break;	       	/* finish search */
	if (next[j].x > xmin  &&  next[j].x < xmax && next[j].y > ymin  &&  next[j].y < ymax )
	  {
	    d = sqrt ((x-next[j].x)*(x-next[j].x) + (y-next[j].y)*(y-next[j].y));
	    if (d < dmin) { dmin = d; p1 = j; }
	  }
      }
    }

  p[0]=p1;
  p[1]=p2;
  p[2]=p3;
  p[3]=p4;
  for (j=0; j<4; j++) if ( p[j] != -999 ) {zaehler++; }
  return (zaehler);
}



void predict (x1, y1, x2, y2, x3, y3)
double x1, y1, x2, y2, *x3, *y3;
{
  *x3 = 2*x2 - x1;
  *y3 = 2*y2 - y1;
  return;
}

void searchquader(X, Y, Z, xr, xl, yd, yu, tpar, cpar)
double X, Y, Z, xr[4], xl[4], yd[4], yu[4];
track_par *tpar;
control_par *cpar;
{
  int k, i;
  double x, y, xz, yz;
  coord_3d quader[8], point;

  /* project quader in image space to define search area */
  for (k=0; k<8; k++)
    {
      quader[k].pnr=k;
    }
  /* calculation of quader points */
  point.pnr=0; point.x=X; point.y=Y; point.z=Z;

    quader[0].x = X + tpar->dvxmin;
    quader[0].y = Y + tpar->dvymin;
    quader[0].z = Z + tpar->dvzmin; /* --- */
    
    quader[1].x = X + tpar->dvxmax;
    quader[1].y = Y + tpar->dvymin;
    quader[1].z = Z + tpar->dvzmin; /* +-- */
    
    quader[2].x = X + tpar->dvxmin;
    quader[2].y = Y + tpar->dvymax;
    quader[2].z = Z + tpar->dvzmin; /* -+- */
    
    quader[3].x = X + tpar->dvxmin;
    quader[3].y = Y + tpar->dvymin;
    quader[3].z = Z + tpar->dvzmax; /* --+ */ //changed by Beat and JuliAn Nov 2008
    
    quader[4].x = X + tpar->dvxmax;
    quader[4].y = Y + tpar->dvymax;
    quader[4].z = Z + tpar->dvzmin; /* ++- */
    
    quader[5].x = X + tpar->dvxmax;
    quader[5].y = Y + tpar->dvymin;
    quader[5].z = Z + tpar->dvzmax; /* +-+ */
    
    quader[6].x = X + tpar->dvxmin;
    quader[6].y = Y + tpar->dvymax;
    quader[6].z = Z + tpar->dvzmax; /* -++ */
    
    quader[7].x = X + tpar->dvxmax;
    quader[7].y = Y + tpar->dvymax;
    quader[7].z = Z + tpar->dvzmax; /* +++ */

  /* calculation of search area */
  for (i = 0; i < cpar->num_cams; i++)
    {
      xr[i]=0;
      xl[i] = cpar->imx;
      yd[i]=0;
      yu[i] = cpar->imy;
      img_coord (i, point.x, point.y, point.z, 
        Ex[i], I[i], G[i], ap[i], *(cpar->mm), &xz,&yz);
      metric_to_pixel (&xz, &yz, xz, yz, cpar);

      for (k=0; k<8; k++) {
        img_coord (i, quader[k].x, quader[k].y, quader[k].z, 
            Ex[i], I[i], G[i], ap[i], *(cpar->mm), &x,&y);
	  metric_to_pixel (&x, &y, x, y, cpar);

	  if (x <xl[i] ) xl[i]=x;
	  if (y <yu[i] ) yu[i]=y;
	  if (x >xr[i] ) xr[i]=x;
	  if (y >yd[i] ) yd[i]=y;
	}
      if (xl[i] < 0 ) xl[i]=0;
      if (yu[i] < 0 ) yu[i]=0;
      if (xr[i] > cpar->imx)
        xr[i] = cpar->imx;
      if (yd[i] > cpar->imy)
        yd[i] = cpar->imy;

      xr[i]=xr[i]-xz;
      xl[i]=xz-xl[i];
      yd[i]=yd[i]-yz;
      yu[i]=yz-yu[i];
    }
  return;

}



void sortwhatfound (foundpix item[16], int *zaehler, int num_cams)
{
  int i,j,m, different;
  foundpix temp;

  different=0;

  /* where what was found */
  for (i=0; i<16; i++)
    for (j=0; j<4; j++)
      for (m=0; m<4; m++)
	if(item[i].ftnr == item[4*j+m].ftnr)
	  {
	    item[i].whichcam[j]=1;
	  }

  /* how often was ftnr found */
  for (i=0; i<16; i++)
    for (j=0; j < num_cams; j++)
      if (item[i].whichcam[j] == 1 && item[i].ftnr !=-1) item[i].freq++;

  /* sort freq */
  for (i=1; i<16; ++i)  for (j=16-1; j>=i; --j)
    {
      if ( item[j-1].freq < item[j].freq )
	{
	  temp = *(item+j-1); *(item+j-1) = *(item+j); *(item+j) = temp;
	}
    }

  for (i=0; i<16; i++)
    for (j=i+1; j<16; j++)
      {
	if (item[i].ftnr == item[j].ftnr || item[j].freq <2)
	  {
	    item[j].freq=0;
	    item[j].ftnr=-1;
	  }
      }

  /* sort freq */
  for (i=1; i<16; ++i)  for (j=16-1; j>=i; --j)
    {
      if ( item[j-1].freq < item[j].freq )
	{
	  temp = *(item+j-1); *(item+j-1) = *(item+j); *(item+j) = temp;
	}
    }
  for (i=0; i<16; ++i) if(item[i].freq != 0) different++;
  *zaehler=different;

}

