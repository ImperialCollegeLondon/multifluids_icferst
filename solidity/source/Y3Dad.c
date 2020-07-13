/* This Y3Dad module was developed for the topological diagnosis of 
 * interconnected fracture patterns and the initial configuration
 * of fracture apertures
 * (1) the topological analysis is based on a ternary-tree search
 *     algorithm
 * (2) the initial apertures can be assigned correlated to trace
 *     lengths or determined empirically by fracture roughenss 
 *     coefficient
 */	// written by Qinghua Lei


/**********************************************************************/
/**********************************************************************/


#include "Yproto.h"
static void Yad3JointDataExtraction(	  /* joint data extraction */
  INT  njoint, INT **i2elto, INT  *i1elty, INT	*i1jtid,
  DBL *d1ncix, DBL  *d1nciy, DBL  *d1nciz, DBL **d2jnrm,
  INT *i1fj2j, DBL **d2fjcx, DBL **d2fjcy, DBL **d2fjcz,
  DBL *d1jtar, DBL **d2fjnm
  )
{ INT ijoint,ifrjt,inode,ijtid,isupward,itmp;
  DBL u[3],v[3],n[3],x[6],y[6],z[6],area,dtmp;

  /* extract joint element information */
  ifrjt=0;
  for(ijoint=0;ijoint<njoint;ijoint++)
  { ijtid=i1jtid[ijoint];
	// coordinates
	for(inode=0;inode<6;inode++)
	{ x[inode]=d1ncix[i2elto[ijtid][inode]];
	  y[inode]=d1nciy[i2elto[ijtid][inode]];
	  z[inode]=d1nciz[i2elto[ijtid][inode]];
	}
	// area and unit normal of joint element
	u[0]=x[1]-x[0]; u[1]=y[1]-y[0]; u[2]=z[1]-z[0];
	v[0]=x[2]-x[0]; v[1]=y[2]-y[0]; v[2]=z[2]-z[0];
	V3DCro(n[0],n[1],n[2],u[0],u[1],u[2],v[0],v[1],v[2]);
	V3DLen(area,n[0],n[1],n[2]); area/=R2;
	V3DNor(dtmp,n[0],n[1],n[2]);
	// resort the node order (upward normal counter-clockwise)
	if(DABS(n[2])>EPSILON)
	{ if(n[2]>R0) isupward=1;
	  else isupward=0;
	}
	else
	{ if(DABS(n[1])>EPSILON)
	  { if(n[1]>R0) isupward=1;
		else isupward=0;
	  }
	  else
	  { if(n[0]>R0) isupward=1;
		else isupward=0;
	} }
	if(isupward==0)
	{ if(i1elty[ijtid]>=1)
	  { itmp=i2elto[ijtid][0];
            i2elto[ijtid][0]=i2elto[ijtid][1];
            i2elto[ijtid][1]=itmp;
            itmp=i2elto[ijtid][5];
            i2elto[ijtid][5]=i2elto[ijtid][4];
            i2elto[ijtid][4]=itmp;
	  }
	  n[0]=-n[0]; n[1]=-n[1]; n[2]=-n[2];
	}
	// store joint element normal
	d2jnrm[ijoint][0]=n[0];
	d2jnrm[ijoint][1]=n[1];
	d2jnrm[ijoint][2]=n[2];
	d1jtar[ijoint]=area;
	// store fracture joint element data
	if(i1elty[ijtid]>1)
	{ for(inode=0;inode<6;inode++)
	  { d2fjcx[ifrjt][inode]=d1ncix[i2elto[ijtid][inode]];
		d2fjcy[ifrjt][inode]=d1nciy[i2elto[ijtid][inode]];
		d2fjcz[ifrjt][inode]=d1nciz[i2elto[ijtid][inode]];
	  }
	  d2fjnm[ifrjt][0]=n[0];
	  d2fjnm[ifrjt][1]=n[1];
	  d2fjnm[ifrjt][2]=n[2];
	  i1fj2j[ifrjt]=ijoint;
	  ifrjt++;
  } }
}

static void Yad3JointConnectivity(	  /* joint connectivity analysis */
  INT    nfrjt, DBL **d2fjcx, DBL **d2fjcy, DBL **d2fjcz,
  DBL **d2fjnm, DBL **d2fjbx, INT **i2fjcn)
{ DBL xmin,xmax,ymin,ymax,zmin,zmax;
  DBL vecdot;
  INT inode,jnode,knode,ijoint,jjoint;
  INT i1,i2,j1,j2;

  /* bounding box */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { xmin=BEPSILON; xmax=-BEPSILON;
	ymin=BEPSILON; ymax=-BEPSILON;
	zmin=BEPSILON; zmax=-BEPSILON;
	for(inode=0;inode<6;inode++)
	{ xmin=MINIM(d2fjcx[ijoint][inode],xmin);
	  xmax=MAXIM(d2fjcx[ijoint][inode],xmax);
	  ymin=MINIM(d2fjcy[ijoint][inode],ymin);
	  ymax=MAXIM(d2fjcy[ijoint][inode],ymax);
	  zmin=MINIM(d2fjcz[ijoint][inode],zmin);
	  zmax=MAXIM(d2fjcz[ijoint][inode],zmax);
	}
	d2fjbx[ijoint][0]=xmin; d2fjbx[ijoint][1]=xmax;
	d2fjbx[ijoint][2]=ymin; d2fjbx[ijoint][3]=ymax;
	d2fjbx[ijoint][4]=zmin; d2fjbx[ijoint][5]=zmax;
  }
  /* initialise connectivity */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { i2fjcn[ijoint][0]=-1;
	i2fjcn[ijoint][1]=-1;
	i2fjcn[ijoint][2]=-1;
  }
  for(ijoint=0;ijoint<nfrjt-1;ijoint++)
  { for(jjoint=ijoint+1;jjoint<nfrjt;jjoint++)
	{ // bounding box test
	  if(d2fjbx[ijoint][0]>d2fjbx[jjoint][1]+EPSILON
	  || d2fjbx[ijoint][1]<d2fjbx[jjoint][0]-EPSILON
	  || d2fjbx[ijoint][2]>d2fjbx[jjoint][3]+EPSILON
	  || d2fjbx[ijoint][3]<d2fjbx[jjoint][2]-EPSILON
	  || d2fjbx[ijoint][4]>d2fjbx[jjoint][5]+EPSILON
	  || d2fjbx[ijoint][5]<d2fjbx[jjoint][4]-EPSILON)
	  { continue;
	  }
	  // connectivity analysis
	  for(inode=0;inode<3;inode++)
	  { i1=inode+1; if(i1>2) i1=0;
		i2=i1+1; if(i2>2) i2=0;
		for(jnode=0;jnode<3;jnode++)
		{ j1=jnode+1; if(j1>2) j1=0;
		  j2=j1+1; if(j2>2) j2=0;
		  // check overlapping
		  if((DABS(d2fjcx[ijoint][i1]-d2fjcx[jjoint][j1])<EPSILON
			&&DABS(d2fjcy[ijoint][i1]-d2fjcy[jjoint][j1])<EPSILON
			&&DABS(d2fjcz[ijoint][i1]-d2fjcz[jjoint][j1])<EPSILON
			&&DABS(d2fjcx[ijoint][i2]-d2fjcx[jjoint][j2])<EPSILON
			&&DABS(d2fjcy[ijoint][i2]-d2fjcy[jjoint][j2])<EPSILON
			&&DABS(d2fjcz[ijoint][i2]-d2fjcz[jjoint][j2])<EPSILON)
		   ||(DABS(d2fjcx[ijoint][i2]-d2fjcx[jjoint][j1])<EPSILON
			&&DABS(d2fjcy[ijoint][i2]-d2fjcy[jjoint][j1])<EPSILON
			&&DABS(d2fjcz[ijoint][i2]-d2fjcz[jjoint][j1])<EPSILON
			&&DABS(d2fjcx[ijoint][i1]-d2fjcx[jjoint][j2])<EPSILON
			&&DABS(d2fjcy[ijoint][i1]-d2fjcy[jjoint][j2])<EPSILON
			&&DABS(d2fjcz[ijoint][i1]-d2fjcz[jjoint][j2])<EPSILON))
		  { // connectivity of the 1st joint
			if(i2fjcn[ijoint][inode]==-1)  // store neighbour
			{ if(i2fjcn[jjoint][jnode]!=-1&&i2fjcn[jjoint][jnode]!=ijoint)
			  { i2fjcn[ijoint][inode]=-2;
			  }
			  else
			  { i2fjcn[ijoint][inode]=jjoint;
			} }
			else if(i2fjcn[ijoint][inode]>=0)  // fracture intersection
			{ for(knode=0;knode<2;knode++)
			  { if(i2fjcn[i2fjcn[ijoint][inode]][knode]==ijoint)
				{ i2fjcn[i2fjcn[ijoint][inode]][knode]=-2;
				  break;
			  } }
			  i2fjcn[ijoint][inode]=-2;
			}
			// connectivity of the 2nd joint
			if(i2fjcn[jjoint][jnode]==-1)  // store neighbour
			{ if(i2fjcn[ijoint][inode]!=-1&&i2fjcn[ijoint][inode]!=jjoint)
			  { i2fjcn[jjoint][jnode]=-2;
			  }
			  else
			  { i2fjcn[jjoint][jnode]=ijoint;
			} }
			else if(i2fjcn[jjoint][jnode]>=0)  // fracture intersection
			{ for(knode=0;knode<2;knode++)
			  { if(i2fjcn[i2fjcn[jjoint][jnode]][knode]==jjoint)
				{ i2fjcn[i2fjcn[jjoint][jnode]][knode]=-2;
				  break;
			  } }
			  i2fjcn[jjoint][jnode]=-2;
			}
  } } } } }
  /* coplanarity check */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { for(inode=0;inode<3;inode++)
	{ if(i2fjcn[ijoint][inode]>=0)
	  { jjoint=i2fjcn[ijoint][inode];
		V3DDot(vecdot,d2fjnm[ijoint][0],d2fjnm[ijoint][1],d2fjnm[ijoint][2],
		  d2fjnm[jjoint][0],d2fjnm[jjoint][1],d2fjnm[jjoint][2]);
		if(DABS(vecdot)<R1-ANGLETOLER)
		{ i2fjcn[ijoint][inode]=-1;
		  for(jnode=0;jnode<2;jnode++)
		  { if(i2fjcn[jjoint][jnode]==ijoint)
			i2fjcn[jjoint][jnode]=-1;
  } } } } }
}

static void Yad3CountChildren(	  /* count real children of a tree-node */
  INT **i2jtcn, INT ijoint, INT *i1fj2s, INT *nreal)
{ INT inode,jjoint;
  (*nreal)=0;
  if(ijoint!=-1)
  { for(inode=0;inode<3;inode++)
	{ jjoint=i2jtcn[ijoint][inode];
	  if(jjoint>=0 && i1fj2s[jjoint]==-1)
	  { i1fj2s[jjoint]=-2;
		(*nreal)++;
  } } }
}

static void Yad3FindChildren(	/* find real children of a tree-node */
  INT **i2jtcn, INT ijoint, INT sectindx,
  INT nreal, INT *i1fj2s, INT *i1child)
{ INT inode, jjoint, ichild;
  ichild=0;
  for(inode=0;inode<3;inode++)
  { jjoint=i2jtcn[ijoint][inode];
	if(jjoint>=0 && i1fj2s[jjoint]<0)
	{ i1child[ichild]=jjoint;
	  i1fj2s[jjoint]=sectindx;
	  ichild++;
  } }
}

static void Yad3LinkJoints2Sections(	  /* identify isolated sections */
  INT **i2fjcn, INT nfrjt, INT *i1fj2s, INT *nsect)
{ INT ijoint,jjoint;
  INT *i1child,*i1parent,*i1subchild;
  INT iparent,nparent,ichild,nchild,ireal,*nreal;

  /* initialisation */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { i1fj2s[ijoint]=-1;
  }
  (*nsect)=0;

  /* ternary tree search */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { if(i1fj2s[ijoint]==-1)
	{ // a new section
	  i1fj2s[ijoint]=*nsect;
	  nchild=3; nparent=1;
	  i1parent=TalINT1(nparent);
	  *i1parent=ijoint;
	  while(nchild>0)
	  { // count number of real children
		nchild=0;
		nreal=TalINT1(nparent);
		for(iparent=0;iparent<nparent;iparent++)
		{ jjoint=i1parent[iparent];
		  Yad3CountChildren(i2fjcn,jjoint,i1fj2s,&nreal[iparent]);
		  nchild+=nreal[iparent];
		}
		// find and mark real children
		if(nchild>0)
		{ ichild=0;
		  i1child=TalINT1(nchild);
		  for(iparent=0;iparent<nparent;iparent++)
		  { if(nreal[iparent]>0)
			{ jjoint=i1parent[iparent];
			  i1subchild=TalINT1(nreal[iparent]);
			  Yad3FindChildren(i2fjcn,jjoint,*nsect,nreal[iparent],i1fj2s,i1subchild);
			  for(ireal=0;ireal<nreal[iparent];ireal++)
			  { i1child[ichild]=i1subchild[ireal];
				ichild++;
			  }
			  FREE(i1subchild);
		  } }
		  nparent=nchild;
		  i1parent=TalINT1(nparent);
		  i1parent=i1child;
	  } }
	  (*nsect)++;
	  FREE(i1parent);
  } }
}

static void Yad3SectionDataExtraction(	  /* extract isolated section data */
  DBL *d1jtar, INT *i1fj2j, INT **i2fjcn, INT   nfrjt, INT *i1fj2s,
  INT nsect, INT *i1scnj, DBL *d1scar, INT *i1scbj)
{ INT inode,ifrjt,isect;
  
  /* initialisation */
  for(isect=0;isect<nsect;isect++)
  { i1scnj[isect]=0;
	d1scar[isect]=R0;
  }
  for(ifrjt=0;ifrjt<nfrjt;ifrjt++)
  { i1scbj[ifrjt]=0;
  }

  /* area calculation */
  for(isect=0;isect<nsect;isect++)
  { for(ifrjt=0;ifrjt<nfrjt;ifrjt++)
	{ if(i1fj2s[ifrjt]==isect)
	  { d1scar[isect]+=d1jtar[i1fj2j[ifrjt]];
		i1scnj[isect]++;
  } } }

  /* boundary joint identification */
  for(ifrjt=0;ifrjt<nfrjt;ifrjt++)
  { for(inode=0;inode<3;inode++)
	{ if(i2fjcn[ifrjt][inode]<0)
	  { i1scbj[ifrjt]=1;
	  }
  } }
}

static void Yad3LinkSections2Fractures(	  /* group sections to form fractures */
  DBL **d2fjcx, DBL **d2fjcy, DBL **d2fjcz, INT **i2fjcn, DBL **d2fjnm,
  DBL **d2fjbx, INT    nfrjt, INT  *i1fj2s, INT  *i1scnj, INT    nsect,
  INT  *i1scbj, INT  *i1sc2f, INT   *nfrac)
{ DBL **d2bdbx;			// bounding boxes of sections
  DBL xmin,xmax,ymin,ymax,zmin,zmax;
  INT isect,jsect,ksect,ijoint,jjoint,inode,jnode,lowID,uppID;
  INT i1,i2,j1,j2;
  DBL vecdot;
  INT connflag,isolfrac;

  /* initialisation */
  for(isect=0;isect<nsect;isect++)
  { i1sc2f[isect]=-1;
  }
  (*nfrac)=0;
  d2bdbx=TalDBL2(nsect,6);

  /* only one section */
  if(nsect==1)
  { i1sc2f[0]=0;
	*nfrac=1;
	return;
  }

  /* bounding box */
  for(isect=0;isect<nsect;isect++)
  { xmin=BEPSILON; xmax=-BEPSILON;
	ymin=BEPSILON; ymax=-BEPSILON;
	zmin=BEPSILON; zmax=-BEPSILON;
	for(ijoint=0;ijoint<nfrjt;ijoint++)
	{ for(inode=0;inode<3;inode++)
	  { xmin=MINIM(d2fjcx[ijoint][inode],xmin);
		xmax=MAXIM(d2fjcx[ijoint][inode],xmax);
		ymin=MINIM(d2fjcy[ijoint][inode],ymin);
		ymax=MAXIM(d2fjcy[ijoint][inode],ymax);
		zmin=MINIM(d2fjcz[ijoint][inode],zmin);
		zmax=MAXIM(d2fjcz[ijoint][inode],zmax);
	} }
	d2bdbx[isect][0]=xmin; d2bdbx[isect][1]=xmax;
	d2bdbx[isect][2]=ymin; d2bdbx[isect][3]=ymax;
	d2bdbx[isect][4]=zmin; d2bdbx[isect][5]=zmax;
  }

  /* section connectivity */
  for(isect=0;isect<nsect-1;isect++)
  { isolfrac=TRUE;
	for(jsect=isect+1;jsect<nsect;jsect++)
	{ // bounding box test
	  if(d2bdbx[isect][0]>d2bdbx[jsect][1]+EPSILON
	  || d2bdbx[isect][1]<d2bdbx[jsect][0]-EPSILON
	  || d2bdbx[isect][2]>d2bdbx[jsect][3]+EPSILON
	  || d2bdbx[isect][3]<d2bdbx[jsect][2]-EPSILON
	  || d2bdbx[isect][4]>d2bdbx[jsect][5]+EPSILON
	  || d2bdbx[isect][5]<d2bdbx[jsect][4]-EPSILON)
	  { continue;
	  }
	  // connected state indicator
	  connflag=FALSE;
	  // recognise connected joint element pairs
	  for(ijoint=0;ijoint<nfrjt;ijoint++)
	  { // not boundary element of the ith section
		if(i1fj2s[ijoint]!=isect||i1scbj[ijoint]==0)
		{ continue;
		}
		for(jjoint=0;jjoint<nfrjt;jjoint++)
		{ // not boundary element of the jth section
		  if(i1fj2s[jjoint]!=jsect||i1scbj[jjoint]==0)
		  { continue;
		  }
		  // bounding box test
		  if(d2fjbx[ijoint][0]>d2fjbx[jjoint][1]+EPSILON
		  || d2fjbx[ijoint][1]<d2fjbx[jjoint][0]-EPSILON
		  || d2fjbx[ijoint][2]>d2fjbx[jjoint][3]+EPSILON
		  || d2fjbx[ijoint][3]<d2fjbx[jjoint][2]-EPSILON
		  || d2fjbx[ijoint][4]>d2fjbx[jjoint][5]+EPSILON
		  || d2fjbx[ijoint][5]<d2fjbx[jjoint][4]-EPSILON)
		  { continue;
		  }
		  // check coplanarity
		  V3DDot(vecdot,d2fjnm[ijoint][0],d2fjnm[ijoint][1],d2fjnm[ijoint][2],
			d2fjnm[jjoint][0],d2fjnm[jjoint][1],d2fjnm[jjoint][2]);
		  if(DABS(vecdot)<R1-ANGLETOLER) continue;
		  // check connectivity
		  for(inode=0;inode<3;inode++)
		  { if(i2fjcn[ijoint][inode]>=0) continue;
			i1=inode+1; if(i1>2) i1=0;
			i2=i1+1; if(i2>2) i2=0;
			for(jnode=0;jnode<3;jnode++)
			{ if(i2fjcn[jjoint][jnode]>=0) continue;
			  j1=jnode+1; if(j1>2) j1=0;
			  j2=j1+1; if(j2>2) j2=0;
			  if((DABS(d2fjcx[ijoint][i1]-d2fjcx[jjoint][j1])<EPSILON
				&&DABS(d2fjcy[ijoint][i1]-d2fjcy[jjoint][j1])<EPSILON
				&&DABS(d2fjcz[ijoint][i1]-d2fjcz[jjoint][j1])<EPSILON
				&&DABS(d2fjcx[ijoint][i2]-d2fjcx[jjoint][j2])<EPSILON
				&&DABS(d2fjcy[ijoint][i2]-d2fjcy[jjoint][j2])<EPSILON
				&&DABS(d2fjcz[ijoint][i2]-d2fjcz[jjoint][j2])<EPSILON)
			   ||(DABS(d2fjcx[ijoint][i2]-d2fjcx[jjoint][j1])<EPSILON
				&&DABS(d2fjcy[ijoint][i2]-d2fjcy[jjoint][j1])<EPSILON
				&&DABS(d2fjcz[ijoint][i2]-d2fjcz[jjoint][j1])<EPSILON
				&&DABS(d2fjcx[ijoint][i1]-d2fjcx[jjoint][j2])<EPSILON
				&&DABS(d2fjcy[ijoint][i1]-d2fjcy[jjoint][j2])<EPSILON
				&&DABS(d2fjcz[ijoint][i1]-d2fjcz[jjoint][j2])<EPSILON))
			  { connflag=TRUE;
				isolfrac=FALSE;
				// update joint connectivity data
				i2fjcn[ijoint][inode]=jjoint;
				i2fjcn[jjoint][jnode]=ijoint;
				break;
			} }
			if(connflag==TRUE) break;
		  }
		  if(connflag==TRUE) break;
		}
		if(connflag==TRUE) break;
	  }
	  if(connflag==TRUE)
	  { // the two sections are connected
		if(i1sc2f[isect]==-1 && i1sc2f[jsect]==-1)
		{ i1sc2f[isect]=*nfrac;
		  i1sc2f[jsect]=*nfrac;
		  (*nfrac)++;
		}
		else if(i1sc2f[isect]==-1)
		{ i1sc2f[isect]=i1sc2f[jsect];
		}
		else if(i1sc2f[jsect]==-1)
		{ i1sc2f[jsect]=i1sc2f[isect];
		}
		else if(i1sc2f[isect]!=i1sc2f[jsect])
		{ lowID=MINIM(i1sc2f[isect],i1sc2f[jsect]);
		  uppID=MAXIM(i1sc2f[isect],i1sc2f[jsect]);
		  i1sc2f[isect]=lowID; i1sc2f[jsect]=lowID;
		  for(ksect=0;ksect<nsect;ksect++)
		  { if(i1sc2f[ksect]==uppID)
			{ i1sc2f[ksect]=lowID;
			}
			else if(i1sc2f[ksect]>uppID)
			{ i1sc2f[ksect]-=1;
		  }	}
		  (*nfrac)--;
	} } }
	// it is an isolated fracture
	if(isolfrac==TRUE && i1sc2f[isect]<0)
	{ i1sc2f[isect]=*nfrac;
	  (*nfrac)++;
	}
	// the last fracture is isolated
	if(connflag==FALSE && isect==nsect-2 && i1sc2f[isect+1]<0)
	{ i1sc2f[isect+1]=*nfrac;
	  (*nfrac)++;
  } }

  /* T type intersection */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { for(inode=0;inode<2;inode++)
	{ if(i2fjcn[ijoint][inode]==-2)
	  { i2fjcn[ijoint][inode]=-1;
  } } }

  /* memory clean */
  FREE(d2bdbx);
}

static void Yad3FractureDataExtraction(	  /* extract fracture data */
  DBL *d1scar, INT nsect, INT *i1sc2f, INT nfrac, DBL *d1frar)
{ INT isect, ifrac;
  
  /* initialisation */
  for(ifrac=0;ifrac<nfrac;ifrac++)
  { d1frar[ifrac] = 0;
  }

  /* area calculation */
  for(isect=0;isect<nsect;isect++)
  { d1frar[i1sc2f[isect]]+=d1scar[isect];
  }
}

static void Yad3LinkJoints2Fractures(	/* link joint elements to fractures */
  INT *i1fj2s, INT   nfrjt, INT *i1sc2f, INT nsect, INT nfrac,
  INT *i1fj2f, INT *i1frnj)
{ INT ifrjt,ifrac;	  // loop variables
  
  for(ifrjt=0;ifrjt<nfrjt;ifrjt++)
  { i1fj2f[ifrjt]=i1sc2f[i1fj2s[ifrjt]];
  }

  /* count number of joints in each fracture */
  for(ifrac=0;ifrac<nfrac;ifrac++)
  { i1frnj[ifrac]=0;
	for(ifrjt=0;ifrjt<nfrjt;ifrjt++)
	{ if(i1fj2f[ifrjt]==ifrac)
	  { i1frnj[ifrac]++;
  } } }
}

static void Yad3EmpFractureJointConfiguration(	/* configure fracture joints empirically */
  INT    nfrjt, INT *i1jtid, INT *i1fj2j, INT *i1fj2s, DBL *d1scar, INT *i1fj2f,
  INT  *i1elpr, DBL *d1pjrc, DBL *d1pjcs, DBL *d1pjsl,
  DBL  *d1jkni, DBL *d1jknc, DBL *d1jksc, DBL *d1jnst, DBL *d1jsst,
  DBL  *d1japi, DBL *d1japc, DBL *d1japh, DBL *d1japr,
  DBL **d2jsdc, DBL *d1jdlc, DBL *d1jsdp, DBL *d1jefl,
  DBL  *d1jjrc, DBL *d1jjcs, DBL *d1jphi)
{ INT ifrjt,ifrac,ijoint;	  // loop variable
  DBL a0;					  // initial aperture
  DBL ar;					  // residual aperture
  DBL vm;					  // maximum normal closure
  INT propID;				  // joint element property id
  DBL JRC0,JCS0,L0,JCSn,JRCn; // joint parameters

  for(ifrjt=0;ifrjt<nfrjt;ifrjt++)
  { // index of joint element
	ijoint=i1fj2j[ifrjt];

	// effective fracture length (i.e. block size)
	// d1jefl[ijoint]=sqrt(d1scar[i1fj2s[ifrjt]]);
	d1jefl[ijoint]=d1scar[i1fj2s[ifrjt]]/0.1;

	// mechanical properties
	propID=i1elpr[i1jtid[ijoint]];
	JRC0=d1pjrc[propID];
	JCS0=d1pjcs[propID];
	L0=d1pjsl[propID];

	// scale-dependent JRC and JCS
	JRCn=JRC0*pow(d1jefl[ijoint]/L0,-0.02*JRC0);
	JCSn=JCS0*pow(d1jefl[ijoint]/L0,-0.03*JRC0);
	// JRCn=JRC0; JCSn=JCS0;
	d1jjrc[ijoint]=JRCn;
	d1jjcs[ijoint]=JCSn;
	d1jphi[ijoint]=R0;

	// initial aperture
	a0=JRC0/50/1e3;

	// closure
	vm=(-0.1032-0.0074*JRCn+1.1350*pow(JCSn/a0*1e-3,-0.2510))*1e-3;
	if(vm<R0) vm=a0*R9/R10;

	// residual aperture
	ar=a0-vm;

	// initial and current aperture
	d1japi[ijoint]=a0;
	d1japc[ijoint]=a0;

	// hydraulic aperture
	d1japh[ijoint]=pow(a0*1e6,R2)/pow(JRCn,2.5)/1e6;

	// residual aperture
	d1japr[ijoint]=ar;

	// initial and current joint normal stiffness
	d1jkni[ijoint]=(-7.15+1.75*JRCn+0.02*JCSn/d1japi[ijoint]*1e-3)*1e9;
	d1jknc[ijoint]=d1jkni[ijoint];

	// current joint shear stiffness
	d1jksc[ijoint]=R0;

	// normal and shear stress
	d1jnst[ijoint]=R0; d1jsst[ijoint]=R0;

	// current and peak shear displacement, and dilation
	d2jsdc[ijoint][0]=R0; d2jsdc[ijoint][1]=R0; d2jsdc[ijoint][2]=R0;
	d1jdlc[ijoint]=R0;
	d1jsdp[ijoint]=d1jefl[ijoint]/500.0*pow(JRCn/d1jefl[ijoint],0.33);
  }
}

static void Yad3UnbrokenJointConfiguration(	/* configure unbroken joints */
  INT  njoint, INT  *i1jtid, INT *i1elty, INT *i1elpr,
  DBL *d1pjrc, DBL  *d1pjcs, DBL *d1pjsl,
  DBL *d1jkni, DBL  *d1jknc, DBL *d1jksc, DBL *d1jnst,
  DBL *d1jsst, DBL  *d1japi, DBL *d1japc, DBL *d1japh,
  DBL *d1japr, DBL **d2jsdc, DBL *d1jdlc, DBL *d1jsdp,
  DBL *d1jefl, DBL  *d1jjrc, DBL *d1jjcs, DBL *d1jphi
  )
{ INT ijoint,ielem,propID;
  DBL JRC0,JCS0;

  for(ijoint=0;ijoint<njoint;ijoint++)
  { ielem=i1jtid[ijoint];
	if(i1elty[ielem]==0)
	{ propID=i1elpr[ielem];
	  JRC0=d1pjrc[propID];
	  JCS0=d1pjcs[propID];
	  d1jefl[ijoint]=d1pjsl[propID];
	  d1jjrc[ijoint]=JRC0;
	  d1jjcs[ijoint]=JCS0;
	  d1jphi[ijoint]=R0;
	  d1japi[ijoint]=d1jjrc[ijoint]/50.0*1e-3;
	  d1japc[ijoint]=d1japi[ijoint];
	  d1japh[ijoint]=R0;
	  d1japr[ijoint]=MAXIM(d1japi[ijoint]/R10,(-0.296-0.0056*d1jjrc[ijoint]+
		2.241*pow(d1jjcs[ijoint]/d1japi[ijoint]*1e-3,-0.245))*1e-3);
	  d1jkni[ijoint]=(-7.15+1.75*d1jjrc[ijoint]+
		0.02*d1jjcs[ijoint]/d1japi[ijoint]*1e-3)*1e9;
	  d1jknc[ijoint]=d1jkni[ijoint];
	  d1jksc[ijoint]=R0; d1jnst[ijoint]=R0;
	  d1jsst[ijoint]=R0; d1jdlc[ijoint]=R0;
	  d2jsdc[ijoint][0]=R0; d2jsdc[ijoint][1]=R0; d2jsdc[ijoint][2]=R0;
	  d1jsdp[ijoint]=d1jefl[ijoint]/500.0*pow(d1pjrc[propID]/d1jefl[ijoint],0.33);
  } }
}

static void Yad3JointConfiguration(	  /* configure joint element properties */
  INT   njoint, INT **i2elto, INT  *i1elpr, INT *i1elty,
  INT  *i1jtid, DBL  *d1jkni, DBL  *d1jknc, DBL *d1jksc,
  DBL  *d1jnst, DBL  *d1jsst, DBL  *d1japi, DBL *d1japc,
  DBL  *d1japh, DBL  *d1japr, DBL **d2jsdc, DBL *d1jdlc,
  DBL  *d1jsdp, DBL  *d1jefl, DBL  *d1jjrc, DBL *d1jjcs,
  DBL  *d1jphi, DBL **d2jnrm, DBL  *d1jtar,
  DBL  *d1ncix, DBL  *d1nciy, DBL  *d1nciz,
  DBL  *d1pefr, DBL  *d1pegt, DBL  *d1pela,
  DBL  *d1pemu, DBL  *d1pjrc, DBL  *d1pjcs, DBL *d1pjsl)
{ INT   ijoint;		// joint loop variable
  INT    nfrjt;		// number of fracture joints
  INT  *i1fj2j;		// fracture joint to joint element indicator
  DBL **d2fjcx;		// fracture joint x coordinates
  DBL **d2fjcy;		// fracture joint y coordinates
  DBL **d2fjcz;		// fracture joint z coordinates
  INT **i2fjcn;		// fracture joint connectivity matrix
  DBL **d2fjbx;		// fracture joint bounding boxes
  DBL **d2fjnm;		// fracture joint upward unit normal
  INT    nsect;		// number of isolated sections
  INT  *i1fj2s;		// fracture joint to isolated section indicator
  INT  *i1scnj;		// number of joints in each isolated section
  INT  *i1scbj;		// isolated section boundary joints
  DBL  *d1scar;		// isolated section area
  INT  *i1sc2f;		// isolated section to fracture indicator
  INT    nfrac;		// number of fractures
  INT  *i1fj2f;		// fracture joint to fracture indicator
  DBL  *d1frar;		// area of each fracture
  INT  *i1frhd;		// fracture head joint
  INT  *i1frnx;		// fracture next joint
  INT  *i1frnj;		// number of joints in each fracture
  INT  ifrac,i;

  /******************************************/
  /*	  Extract joint element info		*/
  /******************************************/
  /* count fracture joints */
  nfrjt=0;
  for(ijoint=0;ijoint<njoint;ijoint++)
  { if(i1elty[i1jtid[ijoint]]>1)
	{ nfrjt++;
  } }

  /* joint data extraction */
  i1fj2j=TalINT1(nfrjt);
  d2fjcx=TalDBL2(nfrjt,6);
  d2fjcy=TalDBL2(nfrjt,6);
  d2fjcz=TalDBL2(nfrjt,6);
  d2fjnm=TalDBL2(nfrjt,3);
  Yad3JointDataExtraction(
	njoint,i2elto,i1elty,i1jtid,
	d1ncix,d1nciy,d1nciz,d2jnrm,
	i1fj2j,d2fjcx,d2fjcy,d2fjcz,
	d1jtar,d2fjnm);

  /******************************************/
  /*	Configure unbroken joint elements	*/
  /******************************************/
  Yad3UnbrokenJointConfiguration(
	njoint,i1jtid,i1elty,i1elpr,
	d1pjrc,d1pjcs,d1pjsl,
	d1jkni,d1jknc,d1jksc,d1jnst,
	d1jsst,d1japi,d1japc,d1japh,
	d1japr,d2jsdc,d1jdlc,d1jsdp,
	d1jefl,d1jjrc,d1jjcs,d1jphi);

  /******************************************/
  /*	Configure fracture joint elements	*/
  /******************************************/
  if(nfrjt>0)
  { /* joint connectivity analysis */
	d2fjbx=TalDBL2(nfrjt,6);
	i2fjcn=TalINT2(nfrjt,3);
	Yad3JointConnectivity(
	  nfrjt ,d2fjcx,d2fjcy,d2fjcz,
	  d2fjnm,d2fjbx,i2fjcn);

	/* joint and isolated section linkage */
	i1fj2s=TalINT1(nfrjt);
	Yad3LinkJoints2Sections(
	  i2fjcn,nfrjt,i1fj2s,&nsect);

	/* isolated section data extraction */
	i1scnj=TalINT1(nsect);
	d1scar=TalDBL1(nsect);
	i1scbj=TalINT1(nfrjt);
	Yad3SectionDataExtraction(
	  d1jtar,i1fj2j,i2fjcn,nfrjt ,i1fj2s,
	  nsect,i1scnj,d1scar,i1scbj);

	/* section and fracture linkage */
	i1sc2f=TalINT1(nsect);
	Yad3LinkSections2Fractures(
	  d2fjcx,d2fjcy,d2fjcz,i2fjcn,d2fjnm,
	  d2fjbx, nfrjt,i1fj2s,i1scnj, nsect,
	  i1scbj,i1sc2f,&nfrac);

	/* fracture data extraction */
	d1frar=TalDBL1(nfrac);
	Yad3FractureDataExtraction(
	  d1scar,nsect,i1sc2f,nfrac,d1frar);

	/* joint element and fracture linkage */
	i1fj2f=TalINT1(nfrjt); i1frhd=TalINT1(nfrac);
	i1frnx=TalINT1(nfrjt); i1frnj=TalINT1(nfrac);
	Yad3LinkJoints2Fractures(
	  i1fj2s,nfrjt ,i1sc2f,nsect,nfrac,
	  i1fj2f,i1frnj);

	/* roughness-correlated */
	Yad3EmpFractureJointConfiguration(
	  nfrjt ,i1jtid,i1fj2j,i1fj2s,d1scar,i1fj2f,
	  i1elpr,d1pjrc,d1pjcs,d1pjsl,
	  d1jkni,d1jknc,d1jksc,d1jnst,d1jsst,
	  d1japi,d1japc,d1japh,d1japr,
	  d2jsdc,d1jdlc,d1jsdp,d1jefl,
	  d1jjrc,d1jjcs,d1jphi);

	/* memory clean */
	FREE(i2fjcn);
	FREE(d2fjbx);
	FREE(i1fj2s);
	FREE(i1scnj);
	FREE(i1scbj);
	FREE(d1scar);
	FREE(i1sc2f);
	FREE(d1frar);
	FREE(i1fj2f);
	FREE(i1frnj);
  }
  FREE(i1fj2j);
  FREE(d2fjcx);
  FREE(d2fjcy);
  FREE(d2fjcz);
  FREE(d2fjnm);
}

/*********************PUBLIC**********************/
void Yad(YDC ydc, YDE yde, YDJ ydj, YDN ydn, YDP ydp)
{ INT ielem,ijoint;		// loop variables
  INT njoint;			// number of joint elements
  INT ntrele;			// number of triangle elements
  
  /* initialisation of joint data structure */
  njoint=0;
  for(ielem=0;ielem<yde->nelem;ielem++)
  { if(yde->i1elty[ielem]!=-1)
	{ njoint++;
  } }
  ntrele=(yde->nelem)-njoint;
  ydj->njoint=njoint;
  ydj->i1jtid=TalINT1(njoint);
  ydj->d1jkni=TalDBL1(njoint);
  ydj->d1jknc=TalDBL1(njoint);
  ydj->d1jksc=TalDBL1(njoint);
  ydj->d1jnst=TalDBL1(njoint);
  ydj->d1jsst=TalDBL1(njoint);
  ydj->d1japi=TalDBL1(njoint);
  ydj->d1japc=TalDBL1(njoint);
  ydj->d1japh=TalDBL1(njoint);
  ydj->d1japr=TalDBL1(njoint);
  ydj->d2jsdc=TalDBL2(njoint,3);
  ydj->d1jdlc=TalDBL1(njoint);
  ydj->d1jsdp=TalDBL1(njoint);
  ydj->d1jefl=TalDBL1(njoint);
  ydj->d1jjrc=TalDBL1(njoint);
  ydj->d1jjcs=TalDBL1(njoint);
  ydj->d1jphi=TalDBL1(njoint);
  ydj->d1jfmd=TalDBL1(njoint);
  ydj->d2jnrm=TalDBL2(njoint,3);
  ydj->d1jtar=TalDBL1(njoint);
  for(ijoint=0;ijoint<njoint;ijoint++)
  { ydj->i1jtid[ijoint]=ntrele+ijoint;
	ydj->d1jkni[ijoint]=BEPSILON;
	ydj->d1jknc[ijoint]=BEPSILON;
	ydj->d1jksc[ijoint]=R0;
	ydj->d1jnst[ijoint]=R0;
	ydj->d1jsst[ijoint]=R0;
	ydj->d1japi[ijoint]=R0;
	ydj->d1japc[ijoint]=R0;
	ydj->d1japh[ijoint]=R0;
	ydj->d1japr[ijoint]=R0;
	ydj->d2jsdc[ijoint][0]=R0;
	ydj->d2jsdc[ijoint][1]=R0;
	ydj->d2jsdc[ijoint][2]=R0;
	ydj->d1jdlc[ijoint]=R0;
	ydj->d1jsdp[ijoint]=R0;
	ydj->d1jefl[ijoint]=R0;
	ydj->d1jjrc[ijoint]=R0;
	ydj->d1jjcs[ijoint]=R0;
	ydj->d1jphi[ijoint]=R0;
	ydj->d1jfmd[ijoint]=R0;
	ydj->d2jnrm[ijoint][0]=R0;
	ydj->d2jnrm[ijoint][1]=R0;
	ydj->d2jnrm[ijoint][2]=R0;
	ydj->d1jtar[ijoint]=R0;
  }

  /* configuration of joint element properties */
  Yad3JointConfiguration(	/* configure joint properties */
	njoint,yde->i2elto,yde->i1elpr,yde->i1elty,
	ydj->i1jtid,ydj->d1jkni,ydj->d1jknc,ydj->d1jksc,
	ydj->d1jnst,ydj->d1jsst,ydj->d1japi,ydj->d1japc,
	ydj->d1japh,ydj->d1japr,ydj->d2jsdc,ydj->d1jdlc,
	ydj->d1jsdp,ydj->d1jefl,ydj->d1jjrc,ydj->d1jjcs,
	ydj->d1jphi,ydj->d2jnrm,ydj->d1jtar,
	ydn->d2nci[0],ydn->d2nci[1],ydn->d2nci[2],
	ydp->d1pefr,ydp->d1pegf,ydp->d1pela,
	ydp->d1pemu,ydp->d1pjrc,ydp->d1pjcs,ydp->d1pjsl);
}
