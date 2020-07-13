/* Copyright (C) 2000, Dr. Antonio Munjiza
 *
 * This code is provided as part of the book entitled "The Combined
 * Finite Discrete Element Method". It is distributed WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other
 * commercial or research or other purpose code is not granted without author's
 * written explicit permission. When results using whole or any part of this code
 * are published, acknowledgement to the author should be made.
 */
/* File  Yproto.h */
#ifndef YDINCL
#include "Y3Dd.h"
#define YDINCL 
#endif 

void Y3dsini(        /*   contact detection                               */
#if NeedFunctionPrototypes
 YSD ysd         /* = discrete element database                       */
 #endif
);

void Y3dscd(        /*   contact detection                               */
#if NeedFunctionPrototypes
 YSDC ysdc,       /* = control database                                */
 YSD ysd,         /* = discrete element database                       */
 YSDP ysdp,       /* = property database                               */
 YSCH ysch        /* = control database                                */
 #endif
);


void Y3dsid(        /*   procreate                                       */
#if NeedFunctionPrototypes
 YSDC ysdc,       /* = control database                                */
 YDOM ydom,       /* =  domain                                         */
 YSD ysd,         /* = discrete element database                       */
 YSDP ysdp,       /* = property database                               */
 YSCH ysch        /* = control database                                */
#endif
);

void Y3dsmove(        /*  solve equations                                  */
#if NeedFunctionPrototypes
 YSD  ysd,
 YSDP ysdp,       /* = property database                               */
 YSDC ysdc
#endif
);

INT Y3dsrd(         /*  read input                                       */
#if NeedFunctionPrototypes
 CHR *namep,      /* - name of the problem i.e. input file             */
 Y3DSD  y3dsd           /* = database                                        */
#endif
);
 
void Y3dsod(         /* output results in space-save format               */
#if NeedFunctionPrototypes
 CHR *namep,      /* - name of the problem i.e. input file             */
 Y3DSD  y3dsd           /* = database                                        */
#endif
);

