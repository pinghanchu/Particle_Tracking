#include <stdlib.h>
#include "ptv.h"
#include "tools.h"

/* declaration */
void qs_coord2d_x ();
void qs_target_y ();
void qs_coord2d_pnr ();
void qs_con ();

FILE *fopen_r (filename)
char *filename;
/*	tries to open a file;
	gives a message, if it cannot open it
	and waits until it has been created 	 */
{
  FILE	*fpr;
  int  	count;

  fpr = fopen (filename, "r");
  if ( ! fpr)
    {
      printf ("could not open %s, please create this file\n", filename);

      /* wait until file can be opened */
      while ( ! fpr)	fpr = fopen (filename, "r");

      /* wait until file really created */
      for (count=0; count<100000; count++);
    }

  return (fpr);
}

void compose_name_plus_nr (basename, str, nr, filename)
char	basename[256], str[256], filename[256];
int		nr;
{
  char	nr_ch[256];

//  if (nr < 10)		sprintf (nr_ch, "00%1d", nr);
//  else if (nr < 100)	sprintf (nr_ch, "0%2d",  nr);
  if (nr < 10)		sprintf (nr_ch, "%1d", nr);
  else if (nr < 100)	sprintf (nr_ch, "%2d",  nr);
  else	sprintf (nr_ch, "%3d",  nr);

  sprintf (filename, "%s%s%s", basename, str, nr_ch);
}

void compose_name_plus_nr_str (basename, str, nr, filename)
char	basename[256], str[256], filename[256];
int		nr;
{
  char	nr_ch[256];

  if (nr < 10)		sprintf (nr_ch, "%1d", nr);
  else if (nr < 100)	sprintf (nr_ch, "%2d",  nr);
  else 	sprintf (nr_ch, "%3d",  nr);

  sprintf (filename, "%s%s%s", basename, nr_ch, str);
}

/* find nearest neighbours */

int kill_in_list ( nr, num, ms_x, ms_y)
int    	nr, num;
int    	ms_x, ms_y;
{
  int 	i, imin = 9999, intx, inty;
  double  x, y, d, dmin = 9999;

  if (zoom_f[nr] > 1)
    {
      printf ("cannot delete point from zoomed image");
      return (0);
    }

  for (i=0; i<num; i++)
    {
      x = (double) ms_x - pix[nr][i].x;
      y = (double) ms_y - pix[nr][i].y;
      d = sqrt (x*x + y*y);
      if (d < dmin)
	{
	  dmin = d; imin = i;
	}
    }
  if (dmin > 10)	return (-1);	       	/*  limit: 10 pixel  */
  intx = (int) pix[nr][imin].x;
  inty = (int) pix[nr][imin].y;

  for (i=imin; i<num; i++)  pix[nr][i] = pix[nr][i+1];

  return (imin);
}



int nearest_neighbour_geo (crd, num, x, y, eps)
coord_2d  crd[];
int       num;
double 	  x, y, eps;
{
  register int	j;
  int	       	j0, dj, pnr = -999;
  double       	d, dmin=1e20, xmin, xmax, ymin, ymax;

  xmin = x - eps;  xmax = x + eps;  ymin = y - eps;  ymax = y + eps;

  /* binarized search for start point of candidate search */
  for (j0=num/2, dj=num/4; dj>1; dj/=2)
    {
      if (crd[j0].x < xmin)  j0 += dj;
      else  j0 -= dj;
    }
  j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to trunc */

  for (j=j0; j<num; j++)		       	/* candidate search */
    {
      if (crd[j].x > xmax)  break;	       	/* finish search */

      if (crd[j].y > ymin  &&  crd[j].y < ymax)
	{
	  d = sqrt ((x-crd[j].x)*(x-crd[j].x) + (y-crd[j].y)*(y-crd[j].y));
	  if (d < dmin)
	    {
	      dmin = d; pnr = j;
	    }
	}
    }
  return (pnr);
}



/***********************************************************************/
/***********************************************************************/

/* sorting routines */

/***********************************************************************/

/* bubble sorts */
void bubble_conlist (item, count)
correspond	*item;
int    		count;
{
	int			i,j;
	correspond	temp;

	for (i=1; i<count; ++i)  for (j=count-1; j>=i; --j)
	{
		if (item[j-1].corr > item[j].corr)
		{
			temp = *(item+j-1);  *(item+j-1) = *(item+j);  *(item+j) = temp;
		}
	}
}

/***********************************************************************/
/***********************************************************************/


/* quicksort algorithms for several issues */

/***********************************************************************/

/* quicksort of 2d coordinates in x-order */

void quicksort_coord2d_x (crd, num)
coord_2d	*crd;
int			num;
{
	qs_coord2d_x (crd, 0, num-1);
}



void qs_coord2d_x (crd, left, right)
coord_2d	*crd;
int			left, right;
{
	register int	i, j;
	double			xm;
	coord_2d		temp;

	i = left;	j = right;	xm = crd[(left+right)/2].x;

	do
	{
		while (crd[i].x < xm  &&  i<right)	i++;
		while (xm < crd[j].x  &&  j>left)	j--;

		if (i <= j)
		{
			temp = crd[i];
			crd[i] = crd[j];
			crd[j] = temp;
			i++;	j--;
		}
	}
	while (i <= j);

	if (left < j)	qs_coord2d_x (crd, left, j);
	if (i < right)	qs_coord2d_x (crd, i, right);
}


/***********************************************************************/

/* quicksort of 2d coordinates in pnr-order */

void quicksort_coord2d_pnr (crd, num)
coord_2d	*crd;
int	       	num;
{
  qs_coord2d_pnr (crd, 0, num-1);
}



void qs_coord2d_pnr (crd, left, right)
coord_2d	*crd;
int    		left, right;
{
  register int	i, j;
  double       	pnrm;
  coord_2d     	temp;

  i = left;	j = right;	pnrm = crd[(left+right)/2].pnr;

  do
    {
      while (crd[i].pnr < pnrm  &&  i<right)	i++;
      while (pnrm < crd[j].pnr  &&  j>left)	j--;

      if (i <= j)
	{
	  temp = crd[i];
	  crd[i] = crd[j];
	  crd[j] = temp;
	  i++;	j--;
	}
    }
  while (i <= j);

  if (left < j)	qs_coord2d_pnr (crd, left, j);
  if (i < right)	qs_coord2d_pnr (crd, i, right);
}




/***********************************************************************/


/* quicksort of targets in y-order */

void quicksort_target_y (pix, num)
target 	*pix;
int    	num;
{
  qs_target_y (pix, 0, num-1);
}




void qs_target_y (pix, left, right)
target 	*pix;
int    	left, right;
{
  register int	i, j;
  double			ym;
  target			temp;

  i = left;	j = right;	ym = pix[(left+right)/2].y;

  do
    {
      while (pix[i].y < ym  &&  i<right)	i++;
      while (ym < pix[j].y  &&  j>left)	j--;

      if (i <= j)
	{
	  temp = pix[i];
	  pix[i] = pix[j];
	  pix[j] = temp;
	  i++;	j--;
	}
    }
  while (i <= j);

  if (left < j)	qs_target_y (pix, left, j);
  if (i < right)	qs_target_y (pix, i, right);
}



/***********************************************************************/


/* quicksort for list of correspondences in order of match quality */
/* 4 camera version */

void quicksort_con (con, num)
n_tupel	*con;
int    	num;
{ 
  if (num > 0) qs_con (con, 0, num-1);
}


void qs_con (con, left, right)
n_tupel	*con;
int    	left, right;
{
  register int	i, j;
  double       	xm;
  n_tupel      	temp;

  i = left;	j = right;	xm = con[(left+right)/2].corr;

  do
    {
      while (con[i].corr > xm  &&  i<right)	i++;
      while (xm > con[j].corr  &&  j>left)	j--;

      if (i <= j)
	{
	  temp = con[i];
	  con[i] = con[j];
	  con[j] = temp;
	  i++;	j--;
	}
    }
  while (i <= j);

  if (left < j)	qs_con (con, left, j);
  if (i < right)	qs_con (con, i, right);
}


/***SORTING ALGORIHTMUS****/

void sort(int n, float a[], int b[])
{
  int flag = 0, i, itemp;
  float ftemp;

  do {
    flag = 0;
    for(i=0; i<(n-1); i++)
      if(a[i] > a[i+1]) {
	ftemp =  a[i];
	itemp =  b[i];
	a[i] = a[i+1];
	b[i] = b[i+1];
	a[i+1] = ftemp;
	b[i+1] = itemp;
        flag = 1;
      }
  }while(flag);
}

/***********************************************************************/



