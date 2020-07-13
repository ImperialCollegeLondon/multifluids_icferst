/*
 * Y3Drscale.c
 *
 *>>	Takes the initial data, and new aperture values to assign on
 *>>	for both dc nodes and surface nodes. Such that the ring-mesh
 *>>	can be scaled (normal to tetrhedral surface) according to
 *>>	aperture data.
 *
 *  Created on: 2 Oct 2017
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



/****************** PUBLIC ************************************/
void Yscale(YDE yde, YDX ydx, YDN ydn, YDXN ydxn)

{
  INT in,inopo;
  INT ielem;
  INT i, j, jono, jn;

  FILE *fp;
  FILE *fpp;
  FILE *fppp;

  DBL *d1ap;
  DBL *d1xap;

  DBL inap=1e-8; //EPSILON;

  d1ap=DBL1NULL;
  d1ap=TalDBL1(ydn->nnopo);
  /* initiliase to 0 in loop so that nodes are not overwritten */
  for (i=0; i<ydn->nnopo; i++)
    {
      d1ap[i]=R0;
    }


  d1xap=DBL1NULL;
  d1xap=TalDBL1(ydn->mnopo);


  /* assigning value based on aperture and joint type to all nodes of a face*/
  for (ielem = yde->nelemi; ielem < yde->nelem; ielem++)
    {
      for (i = 0; i < 3; i++)
	{
	  in = yde->i2elto[ielem][i];  //each ielem should have unique 4 nodes

	  /* main loop to assign apertures on nodes */
	  if (yde->i1elty[ielem] == 1) // && d1ap[in] <= inap) 	//boundary and unbroken
	    {
	      d1ap[in] = inap;     				//minimum element thickness for the ring
	    }
	  else if (yde->i1elty[ielem] == 2)
	    {
	      /* assigning intial aperture based on constant value */
	      d1ap[in] = ydx->dthick;
	    }

	  else if (yde->i1elty[ielem] > 2)
	    {
	      /* new aperture from change in position */
	      d1ap[in] = yde->d2xap[ielem][i];
	    }


	}
    }




  // looop over all the nodes around a field - > discont to cont
  for(inopo=0;inopo<ydxn->ncmno;inopo++)
      {
        d1xap[inopo]=R0;
        jono=0;

        for(j=0;j<100;j++)
  	{
  	  jn=ydxn->i2xbjno[inopo][j];

  	  if(jn>=0)
  	    {

  	      d1xap[inopo]=d1xap[inopo]+d1ap[jn];
  	      jono+=1;
  	    }
  	}
        if(jono>0)
  	{
  	  d1xap[inopo]=d1xap[inopo]/jono;
          ydxn->d1shap[inopo]=d1xap[inopo];
  	}
      }




/*
  // 2 PRINT
  fpp = fopen ("joints_top_2.txt", "w+");
  for (i = yde->nelemi; i < yde->nelem; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  in = yde->i2elto[i][j];
	  fprintf (fpp, " %i -%i- %i--- %e ----%i \n", i, j, in, d1ap[in],
		   yde->i1elty[i]);
	}
    }
  fclose (fpp);

// 3 PRINT
  fp = fopen ("joints_top_3.txt", "w+");
  for (i = 0; i < ydxn->ncmno; i++)
    {
      fprintf (fp, " %i ----- %e ----%e \n", i, d1xap[i], ydxn->d1shap[i]);
    }
  fclose (fp);
*/




  FREE(d1xap);
  FREE(d1ap);
  }
/*************************** END PUBLIC ********************************************/


