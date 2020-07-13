/*! \file Y3D.c
 *  \brief main for Y3D model
 *
 *  \todo
 *  finish documentation
 *
 *  Copyright (C) 2008, Queen Mary University of London (QMUL) & 
 *  Imperial College of Science, Technology and Medicine (ICSTM).
 *  All rights reserved. Implemented by Prof Antonio Munjiza & 
 *  Dr Jiansheng Xiang.
 *
 *  This code is part of the Virtual Geoscience Workbench (VGW) 
 *  developed jointly by ICSTM and QMUL through two related parallel 
 *  projects at ICSTM and QMUL respectively funded by EPSRC. 
 *
 *  This code is provided by copyright holders under the GNU Lesser 
 *  General Public License (LGPL). It is open source code; you can 
 *  redistribute it and/or modify it under the terms of the GNU Lesser 
 *  General Public License version 3.  
 *  
 *  This code is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty 
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
 *  the GNU Lesser General Public License for more details,
 *  http://www.gnu.org/licenses/lgpl-3.0.txt. 
 *
 *  You should have received a copy of the GNU Lesser General Public 
 *  License along with this code; if not, write to:
 *
 *  Dr Jiansheng Xiang   <j.xiang@imperial.ac.uk>    \n
 *  Prof Antonio Munjiza <a.munjiza@qmul.ac.uk>      \n
 *  Dr John-Paul Latham  <j.p.latham@imperial.ac.uk> \n
 *
 */


/**********************************************************************/
/**********************************************************************/


#include "Yproto.h"

/*Z for development only */
#define SHOW_STATS(t) \
{ \
  double sum; \
  t[0]/=ydc->mcstep; t[1]/=ydc->mcstep; t[2]/=ydc->mcstep; t[3]/=ydc->mcstep; t[4]/=ydc->mcstep; t[5]/=ydc->mcstep;\
  sum = t[0] + t[1] + t[2] + t[3] + t[4] + t[5]; \
  fprintf(stderr,"average time per module (units=%f s)\n", 1.0/CLOCKS_PER_SEC); \
  fprintf(stderr,"Ymd:  %08.3f (%07.3f%%)\n", t[0], 100.0*t[0]/sum); \
  fprintf(stderr,"Yfd:  %08.3f (%07.3f%%)\n", t[1], 100.0*t[1]/sum); \
  fprintf(stderr,"Ycd:  %08.3f (%07.3f%%)\n", t[2], 100.0*t[2]/sum); \
  fprintf(stderr,"Yid:  %08.3f (%07.3f%%)\n", t[3], 100.0*t[3]/sum); \
  fprintf(stderr,"Ysd:  %08.3f (%07.3f%%)\n", t[4], 100.0*t[4]/sum); \
  fprintf(stderr,"Yod:  %08.3f (%07.3f%%)\n", t[5], 100.0*t[5]/sum); \
}


/**********************************************************************/
/* PRIVATE                                                            */
/**********************************************************************/


static void Ycopyright(void) 
{
  CHRw(stdout, "******************************************************");
  CHRw(stdout, "****************\n");
  CHRw(stdout, "** Copyright (C) 2008,");
  CHRw(stdout, "                                              **\n");
  CHRw(stdout, "** Queen Mary University of London (QMUL) & Imperial "); 
  CHRw(stdout, "College        **\n");
  CHRw(stdout, "** of Science, Technology and Medicine (ICSTM). "); 
  CHRw(stdout, "All rights reserved.**\n");
  CHRw(stdout, "** Implemented for you by Prof Antonio Munjiza & ");
  CHRw(stdout, "Dr Jiansheng Xiang **\n");
  CHRwcr(stdout);
  CHRwcr(stdout);
  CHRwcr(stdout);
  CHRw(stdout, "* This code is part of the Virtual Geoscience ");
  CHRw(stdout, "Workbench (VGW) developed\n");
  CHRw(stdout, "* jointly by ICSTM and QMUL through two related ");
  CHRw(stdout, "parallel projects at \n");
  CHRw(stdout, "* ICSTM and QMUL respectively funded by EPSRC. \n");
  CHRw(stdout, "*  \n");
  CHRw(stdout, "* This code is provided by copyright holders under");
  CHRw(stdout, " the GNU lesser \n");
  CHRw(stdout, "* General Public License (LGPL)."); 
  CHRw(stdout, " It is open source code; you can \n");
  CHRw(stdout, "* redistribute it and/or modify it under the terms ");
  CHRw(stdout, "of the GNU Lesser\n");
  CHRw(stdout, "* General Public License version 3.  \n");
  CHRw(stdout, "*  \n");
  CHRw(stdout, "* This code is distributed in the hope that it will ");
  CHRw(stdout, "be useful, \n");
  CHRw(stdout, "* but WITHOUT ANY WARRANTY; without even the implied ");
  CHRw(stdout, "warranty \n");
  CHRw(stdout, "* of MERCHANTABILITY or FITNESS FOR A PARTICULAR ");
  CHRw(stdout, "PURPOSE. See \n");
  CHRw(stdout, "* the GNU Lesser General Public License");
  CHRw(stdout, " for more details,\n");
  CHRw(stdout, "* http://www.gnu.org/licenses/lgpl-3.0.txt. \n");
  CHRw(stdout, "*  \n");
  CHRw(stdout, "* You should have received a copy of the GNU Lesser ");
  CHRw(stdout, "General Public\n");
  CHRw(stdout, "* License along with this code; if not, write to \n");
  CHRw(stdout, "* Dr Jiansheng Xiang Prof Antonio Munjiza ");
  CHRw(stdout, "or Dr John-Paul Latham \n");
  CHRw(stdout, "* j.xiang@imperial.ac.uk a.munjiza@qmul.ac.uk \n");
  CHRw(stdout, "* or j.p.latham@imperial.ac.uk  \n");
  CHRw(stdout, "* *************************************************");
  CHRw(stdout, "****************** *\n");
}


/**********************************************************************/
/* PUBLIC                                                             */
/**********************************************************************/


/*! \brief main program Y3D
 *  \param[in] argc number of command line arguments
 *  \param[in] argv array of command line arguments
 *  \par Details:
 *  
 *  The program takes one argument which is the file name of the 
 *  model definition file (*.Y3D).
 *  
 *  The model input is parsed using Yrd(), and runs are defined and 
 *  executed until there is no more input in the file.
 *
 *  Each run calls the following functions in sequence for each time 
 *  step:
 *
 *  \b Yfd() // calculate nodal forces      \n
 *  \b Ycd() // calculate contact detection \n
 *  \b Yid() // detect contact interaction  \n
 *  \b Ysd() // solve equations             \n
 *  \b Yod() // output results              \n
 *
 *
 *
 */
main(argc, argv)
  INT argc; char **argv;
{
  /* display copyright */
  Ycopyright();

  if(argv[1]!=NULL)
  { 
    TIMER t;                 /*Z t is simulation start time           */
    CHR c1name[300];         /*2 c1name is name of the input file     */

//-PY changed it for 3D_fracture_coupling_with_multiphase
    CHR c1name1[300];         /*2 c1name is name of the input file     */
    CHR c1name2[300];         /*2 c1name is name of the input file     */

    struct YD_struct yd;     /*2 yd is Y database                     */
    YDK ydk = &yd.ydk;       /*2 ydk is Y constants database          */
    YDC ydc = &yd.ydc;       /*2 ydc is Y control database            */
    YDE yde = &yd.yde;       /*2 yde is Y element database            */
    YDI ydi = &yd.ydi;       /*2 ydi is Y interaction database        */
    YDN ydn = &yd.ydn;       /*2 ydn is Y node database               */
    YDO ydo = &yd.ydo;       /*2 ydo is Y output database             */
    YDP ydp = &yd.ydp;       /*2 ydp is Y property database           */
    YDB ydb = &yd.ydb;       /*2 ydb is Y boundary condition database */

//-PY changed it for 3D_fracture_coupling_with_multiphase
    YPAR ypar = &yd.ypar;      
    YSP ysp = &yd.ysp;       /*2 ydb is Y boundary condition database */

    YDX ydx = &yd.ydx;       /*2 ydx is Y surface boundary database */
    YDXN ydxn = &yd.ydxn;    /*2 ydxn is Y surface node database */

    DBL ts[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    TIMER t1, t2;

	DBL *d1denf;
  DBL *d1absorp;
  DBL *f1dt,*f1dt1,*f1dt2;
  DBL *d1hy_CV;

  INT index,j;


  index=1;


    /*4 get name of the problem */
    CHRcpy(c1name, argv[1]);


	CHRcpy(c1name1, argv[2]);
	CHRcpy(c1name2, argv[3]);


    /*Z initialisation (NOP: automatic variables are already zeroed) */
    ydc->finp = FILENULL; 
    ydc->fcheck = FILENULL;

    /*Z record start time */
    t = TIME;

    /*Z initilalise constants */
    Yinit(ydk); 
//-PY changed it for 3D_fracture_coupling_with_multiphase
    CHRw(stdout, "TEST OUTPUT");
    CHRwcr(stdout);

    /*4 Process while any input */
//-PY changed it for 3D_fracture_coupling_with_multiphase
    while(Yrd(c1name,c1name1,c1name2, &yd) > 0)
    { 
//		Yrd(c1name,c1name1,c1name2, &yd);

		  d1denf=TalDBL1(ydn->mnopo);
  d1absorp=TalDBL1(ydn->mnopo);
  f1dt=TalDBL1(ydn->mnopo);
  f1dt1=TalDBL1(ydn->mnopo);
  f1dt2=TalDBL1(ydn->mnopo);

	for(j=0;j<ydn->mnopo-1;j++)
	{
	
	d1denf[j]=R0;
	d1absorp[j]=R0;
	f1dt[j]=R0;
	f1dt1[j]=R0;
	f1dt2[j]=R0;

	}

      /*Z store current simulation starting timestamp */
//      DBL dctime = ydc->dctime; 

      /*Z configure per-simulation constants */
      Yconfig(ydk, ydc, ydb);

      CHRw(stdout, "NEW INPUT");
      CHRwcr(stdout);
      
      /* added by LGuo for coupling */
//      Ysurface(yde, ydx, ydn);
	 Yring(ydx,ydxn, ydn,ysp);
      /* added by LGuo */
      Ymd(yde, ydi, ydn, ydp);

      /* added by LGuo for coupling */
      YSFtoJOINT(yde,ydn,ydx,ydxn,ydp);

//	 Yring(ydx, ydn,ysp);

      /*Z run current simulation */
      for(ydc->ncstep=0; ydc->ncstep<ydc->mcstep; ydc->ncstep++)
      { 
                                                    t1 = TIME;
        /*LG mesh elements                    */
        //Ymd(yde, ydi, ydn, ydp);                    t2 = TIME; ts[0] += (t2-t1); t1 = t2;
        //CHRw(stdout, "Ymd PASSED");
        //CHRwcr(stdout);

        /*4 calculate nodal forces            */
        Yfd(yde, ydn, ydp, ydk,ydb,ydi, ydc,ydx);                    t2 = TIME; ts[1] += (t2-t1); t1 = t2;
        //CHRw(stdout, "Yfd PASSED");
        //CHRwcr(stdout);

        /*4 calculate contact detection       */
        Ycd(ydc, yde, ydi, ydn, ydp);               t2 = TIME; ts[2] += (t2-t1); t1 = t2;        
        //CHRw(stdout, "Ycd PASSED");
        //CHRwcr(stdout);


//-PY changed it for 3D_fracture_coupling_with_multiphase
        /*4 detect contact interaction        */
        Yid(ydc, yde, ydi, ydn, ydp);               t2 = TIME; ts[3] += (t2-t1); t1 = t2;             
        //CHRw(stdout, "Yid PASSED");
        //CHRwcr(stdout);

        /*4  solve equations                  */
        Ysd(ydc,yde,ydn,ydo,ydb,ydk,ydp,d1absorp,d1denf,index,f1dt,f1dt1,f1dt2,ypar,ysp);
		t2 = TIME; ts[4] += (t2-t1); t1 = t2;
        //CHRw(stdout, "Ysd PASSED");
        //CHRwcr(stdout);

        /*4  output results                   */
        Yod(c1name, &yd);                           t2 = TIME; ts[5] += (t2-t1); t1 = t2;                    
        //CHRw(stdout, "Yod PASSED");
        //CHRwcr(stdout);

        /*4  update time (Z - more accurate)  */
        //Ywrs(c1name,&yd);

	/* added by LGuo for coupling */
	YSFmesh(ydc,ydx,ydxn,ydn,ydp,yde);
		/*4  update time                   */

        ydc->dctime = ydc->dctime + (ydc->dcstec*(ydc->ncstep + 1));
      }
    }

    CHRw(stderr, "   ***** Y3D HAS ORDERLY FINISHED *****\n"); 

    /*Z display elapsed time */
    TIMERw(stderr, t, 60, "minutes\n");

    /*Z display average time per module */
    SHOW_STATS(ts)
  }
  else
  { 
    CHRw(stdout, "Double click the data file to run the program.\n");
  }
}


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/


