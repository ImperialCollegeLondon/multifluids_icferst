/* file Ytypes.h  Y  types of objects */
/* properties */
#define YIPROPMAX  1000 /* maximum possible number of propertiest */ 
/* NODES */
#define YTN2MEC    -1  /* 2D mechanical x+y d.o.f. node   */
#define YTN2RIG    -2  /* 2D mechanical x+y d.o.f. node   */
#define YTN3MEC   -11  /* 3D mechanical x+y+z d.o.f. node */

/* ELEMENTS */
#define YTE2TRIELS  1 /* plain stress triangle */
#define YTE2TRIRIG  2 /* plain stress triangle */
#define YTE2JOINTS  3 /* joint                 */
#define YTE2TRISOF  4 /* plain stress triangle */

#define YTE3TETELS 11 /* elastic tetrahedra    */
#define YTE3TET4ELS 13 /*!< elastic 4-noded tetrahedra   */
#define YTE3TET4JOINT 15 /*!< joint element for 4-node tetrahedra */

# define YTEG2RAD5     21 /* 2D 0.5 radius grain */

/* FIELDS */
#define YFLDDEF 0 /* default = nothing */
#define YFLDSXX 1 /* stress sigma_xx   */
#define YFLDSXY 2 /* stress sigma_xy   */
#define YFLDSYY 3 /* stress sigma_yy   */
#define YFLDSZZ 4 /* stress sigma_yy   */
#define YFLDSZX 5 /* stress sigma_yy   */
#define YFLDSZY 6 /* stress sigma_yy   */
#define YFLDVEL 7 /* velocity          */
#define YFLDVEX 8 /* velocity x        */
#define YFLDVEY 9 /* velocity y        */
#define YFLDVEZ 10 /* velocity z       */
#define YFLDAEL 11 /* isotropic elastic damage       */
#define YFLBORM 12 /* borhole pressure               */
#define YFLBORP 13 /* borhole mass                   */
#define YFLBORV 14 /* borhole spec. volume           */
#define YFLEK   15 /* total kinetic energy           */
#define YFLG2PR 16 /* G2 pressure                    */
