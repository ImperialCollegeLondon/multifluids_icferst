/*
 * Ypref.c
 *
 *>>	Designed to calculate mechanical aperture as local variation
 *>> 	in length due to fracture opening
 *
 *  Created on: 31 Oct 2017
 *      Author: ao1414
 *
 *      By Asiri Indika Obeysekara (a.obeysekara14@imperial.ac.uk)
 *      IMPERIAL COLLEGE LONDON
 *
 *      As part of development in Solidity Project FEMDEM solver
 *      Fluid-Solid coupling for fluid-driven fracturing project
 *
 * This code is provided by copyright holders under the GNU Lesser
 * General Public License (LGPL). It is open source code; you can
 * redistribute it and/or modify it under the terms of the GNU Lesser
 * General Public License version 3.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
 * the GNU Lesser General Public License for more details,
 * http://www.gnu.org/licenses/lgpl.txt.
 *
 *
 *
 *
 *
 */



#include "Yproto.h"



/**********************************************************************/

static void Yapm(INT nelem, INT nelemi, INT **i2elto, INT *i1elty, DBL *d1ncix, DBL *d1nciy, DBL *d1nciz,
		 DBL *d1nccx, DBL *d1nccy, DBL *d1nccz, DBL**d2xap, DBL dthick)
{
  INT i, j, k;
  INT ielem, jelem;
  INT i0, i1, i2, i3, i4, i5;

  /* length change in each node of joint */
  DBL d1x, d1y, d1z;
  DBL d2x, d2y, d2z;
  DBL d3x, d3y, d3z;

  DBL *d1xap1, *d1xap2, *d1xap3;

  d1xap1 = DBL1NULL;
  d1xap1 = TalDBL1 (nelem);

  d1xap2 = DBL1NULL;
  d1xap2 = TalDBL1 (nelem);

  d1xap3 = DBL1NULL;
  d1xap3 = TalDBL1 (nelem);

  FILE *fp;
  fp = fopen ("mechanical_opening.txt", "w+");

  for (ielem = nelemi; ielem < nelem; ielem++)  //only joints (so use nelemi) - Totally wrong
    {
      i0 = i2elto[ielem][0];
      i1 = i2elto[ielem][1];
      i2 = i2elto[ielem][2];
      i3 = i2elto[ielem][3];
      i4 = i2elto[ielem][4];
      i5 = i2elto[ielem][5];

      fprintf (fp, "-------- %i ------%i--- \n\n", ielem, i1elty[ielem]);

      /* new fractures */
      if (i1elty[ielem] > 2)
	{
	/*  d1x = d1ncix[in0] - d1nccx[in0];
	  d1y = d1ncix[in0] - d1nccx[in0];
	  d1z = d1ncix[in0] - d1nccx[in0];

	  d2x = d1ncix[in1] - d1nccx[in1];
	  d2y = d1ncix[in1] - d1nccx[in1];
	  d2z = d1ncix[in1] - d1nccx[in1];

	  d3x = d1ncix[in2] - d1nccx[in2];
	  d3y = d1ncix[in2] - d1nccx[in2];
	  d3z = d1ncix[in2] - d1nccx[in2]; */


	       d1x=RP5*(d1nccx[i0]-d1nccx[i5]);
	       d1y=RP5*(d1nccy[i0]-d1nccy[i5]);
	       d1z=RP5*(d1nccz[i0]-d1nccz[i5]);

	       d2x=RP5*(d1nccx[i1]-d1nccx[i4]);
	       d2y=RP5*(d1nccy[i1]-d1nccy[i4]);
	       d2z=RP5*(d1nccz[i1]-d1nccz[i4]);

	       d3x=RP5*(d1nccx[i2]-d1nccx[i3]);
	       d3y=RP5*(d1nccy[i2]-d1nccy[i3]);
	       d3z=RP5*(d1nccz[i2]-d1nccz[i3]);



	  fprintf (fp, "opening is: 1x=[%e] 1y=[%e] 1z=[%e] \n", d1x, d1y, d1z);
	  fprintf (fp, "opening is: 2x=[%e] 2y=[%e] 2z=[%e] \n", d2x, d2y, d2z);
	  fprintf (fp, "opening is: 3x=[%e] 3y=[%e] 3z=[%e] \n", d3x, d3y, d3z);

	  V3DLen(d1xap1[ielem], d1x, d1y, d1z);
	  V3DLen(d1xap2[ielem], d2x, d2y, d2z);
	  V3DLen(d1xap3[ielem], d3x, d3y, d3z);

	  fprintf (fp, "opening length 1[%e] - 2[%e] - 3[%e] \n\n",
		   d1xap1[ielem], d1xap2[ielem], d1xap3[ielem]);
	  fprintf (fp, "----------------- \n\n");

	  if (d1xap1[ielem] > dthick)
	    {
	      d2xap[ielem][0] = d1xap1[ielem];
	    }
	  else if (dthick>=d1xap1[ielem] > EPSILON)
	    {
	      d2xap[ielem][0] = dthick;
	    }
	  else
	    {
	      d2xap[ielem][0] = EPSILON;
	    }

	  if (d1xap2[ielem] > dthick)
	    {
	      d2xap[ielem][1] = d1xap2[ielem];
	    }
	  else if (dthick>=d1xap2[ielem] > EPSILON)
	    {
	      d2xap[ielem][1] = dthick;
	    }
	  else
	    {
	      d2xap[ielem][1] = EPSILON;
	    }

	  if (d1xap3[ielem] > dthick)
	    {
	      d2xap[ielem][2] = d1xap3[ielem];
	    }
	  else if (dthick>=d1xap3[ielem] > EPSILON)
	    {
	      d2xap[ielem][2] = dthick;
	    }
	  else
	    {
	      d2xap[ielem][2] = EPSILON;
	    }

	}
    }

  fclose (fp);


  FREE(d1xap1);
  FREE(d1xap2);
  FREE(d1xap3);

}






/*********************PUBLIC*****************************/

void Yap (YDE yde, YDN ydn, YDX ydx)
{

  /* calcualting mechanial aperture */
  Yapm(yde->nelem, yde-> nelemi, yde->i2elto, yde->i1elty, ydn->d2nci[0], ydn->d2nci[1],
       ydn->d2nci[2],ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2ncc[2], yde->d2xap, ydx->dthick);

}
