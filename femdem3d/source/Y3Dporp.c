/*
 * 	Y3Dporp.c
 *
 *>>	Designed to take pore-pressure from fluid model and convert from
 *>> 	continuous mesh to discontinuous mesh, in order to modify the cauchy
 *>> 	stress tensor in Y3Dfd.c
 *
 *  Created on: 13 Nov 2017
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

/*************************** PUBLIC ********************************************/
void Ypore(yde, ydn, ydx, d1flpore /** calculating poro fluid pressure **/)

YDE yde;YDN ydn;YDX ydx; DBL *d1flpore;
{
	INT i0, i1, i, j, k, ielem, in, jn, kn;

	FILE *fp;
	DBL *d1tpore;


	fp=fopen("porp.txt", "w+");

	d1tpore = TalDBL1(ydn->nnopo);

	/* for all nodes in the full mesh */
	for (i = 0; i < ydn->nnopo; i++) {
		d1tpore[i] = R0;

		//map from continous mesh (d1flpore) to [temp] discontinous mesh (d1tpore)
		j=ydn->i1ntoC2D[i];
		fprintf(fp, "dnode:%i  ---- cnode=%i --- pressure=%e \n\n", i,j, d1tpore[i]);

		d1tpore[i]=d1flpore[j];
	}
	fclose(fp);

	// map from [temp] nodes to elements
	for (ielem = 0; ielem < yde->nelem; ielem++) {
		for (i = 0; i < 3; i++)      // Nodes to elements
		{
			in = yde->i2elto[i][ielem];
			yde->d1npore[ielem] = d1tpore[in];
		}
	}

	FREE(d1tpore);
}
/*************************** END PUBLIC ********************************************/
