/*! \file Ytypes.h
 *  \brief Y types of objects
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


#ifndef YTYPES_H
#define YTYPES_H


/**********************************************************************/
/**********************************************************************/


/*! \name Properties */
/*! \{ */
#define YIPROPMAX  1000 /*!< maximum possible number of propertiest */ 
/*! \} */


/**********************************************************************/
/**********************************************************************/


/*! \name NODES */
/*! \{ */
#define YTN2MEC    -1  /*!< 2D mechanical x+y d.o.f. node   */
#define YTN2RIG    -2  /*!< 2D mechanical x+y d.o.f. node   */
#define YTN3MEC   -11  /*!< 3D mechanical x+y+z d.o.f. node */
/*! \} */


/**********************************************************************/
/**********************************************************************/


/*! \name ELEMENTS */
/*! \{ */
#define YTE2TRIELS   1 /*!< plain stress triangle */
#define YTE2TRIRIG   2 /*!< plain stress triangle */
#define YTE2JOINTS   3 /*!< joint                 */
#define YTE2TRISOF   4 /*!< plain stress triangle */
#define YTEG2RAD5   21 /*!< 2D 0.5 radius grain   */

//#define YTE3TETELS  10 /*!<    linear elastic tetrahedra */
#define YTE3TET10ELS 11 /*!< nonlinear elastic tetrahedra */
#define YTE3TET4ELS 13 /*!< elastic 4-noded tetrahedra   */
#define YTE3TET4SOF 14 /*!< softening 4-node tetrahedra */
#define YTE3TET4JOINT 15 /*!< joint element for 4-node tetrahedra */


/*! \} */


/**********************************************************************/
/**********************************************************************/


/*! \name FIELDS */
/*! \{ */
#define YFLDDEF 0 /*!< default = nothing */
#define YFLDSXX 1 /*!< stress sigma_xx   */
#define YFLDSXY 2 /*!< stress sigma_xy   */
#define YFLDSYY 3 /*!< stress sigma_yy   */
#define YFLDSZZ 4 /*!< stress sigma_yy   */
#define YFLDSZX 5 /*!< stress sigma_yy   */
#define YFLDSZY 6 /*!< stress sigma_yy   */
#define YFLDVEL 7 /*!< velocity          */
#define YFLDVEX 8 /*!< velocity x        */
#define YFLDVEY 9 /*!< velocity y        */
#define YFLDVEZ 10 /*!< velocity z       */
#define YFLDAEL 11 /*!< isotropic elastic damage       */
#define YFLBORM 12 /*!< borhole pressure               */
#define YFLBORP 13 /*!< borhole mass                   */
#define YFLBORV 14 /*!< borhole spec. volume           */
#define YFLEK   15 /*!< total kinetic energy           */
#define YFLG2PR 16 /*!< G2 pressure                    */
/*! \} */


/**********************************************************************/
/**********************************************************************/


#endif /* YTYPES_H */ 


/**********************************************************************/
/* EOF                                                                */
/**********************************************************************/

