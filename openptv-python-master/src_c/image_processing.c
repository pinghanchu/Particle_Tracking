/****************************************************************************

Routine:	       	image_processing.c

Author/Copyright:      	Hans-Gerd Maas

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	1988
	
Description:	       	different image processing routines ...
	
Routines contained:  
		       	lowpass_n:	n*n local average,, fast
		       			computation time independent from n
		       	enhance:	enhances gray value spectrum to 0..255,
		      			some extreme gray values are cut off
		       	mark_img:	reads image and pixel coordinate set,
			       		marks the image in a certain (monocomp)
			       		color and writes image

****************************************************************************/

#include "ptv.h"

/*------------------------------------------------------------------------
	Subtract mask, Matthias Oswald, Juli 08
  ------------------------------------------------------------------------*/
void subtract_mask (img, img_mask, img_new) 

unsigned char	*img, *img_mask, *img_new;

{
	register unsigned char 	*ptr1, *ptr2, *ptr3;
	int i;
	
	for (i=0, ptr1=img, ptr2=img_mask, ptr3=img_new; i<imgsize; ptr1++, ptr2++, ptr3++, i++)
    {
      if (*ptr2 == 0)  *ptr3 = 0;
      else  *ptr3 = *ptr1;
    }
 }
