*realformat "%20.10e"
*intformat "%7i"

*Set Cond Volume-Constraints *nodes *NoRepeat
*Set Var NConstraints0=CondNumEntities
*Set Cond Surface-Constraints *nodes *NoRepeat
*Add Cond Point-Constraints *nodes *NoRepeat
*Set Var NConstraints1=CondNumEntities
*set Cond Line-Constraints *nodes *NoRepeat
*Set Var NConstraints2=CondNumEntities
*Set Var NConstraints=NConstraints0+NConstraints1+NConstraints2

*Set Cond FACENODES *nodes *NoCanRepeat
*Set Var NFACENODES=CondNumEntities
*if(IsQuadratic(int)==0)
*Set Var nnumber=5
*Set Var eletyp=13
*endif
*if(IsQuadratic(int)==1)
*Set Var nnumber=11
*Set Var eletyp=11
*Set Var ntri=6
*endif
*Set Var mnumber=nnumber+1

*if((strcmp(GenData(15),"Linear")==0))
*Set Var nnumber=5
*Set Var eletyp=13
*Set Var ntri=3
*Set Var nnsurf=3
*endif
*if((strcmp(GenData(15),"Quadratic")==0))
*Set Var nnumber=11
*Set Var eletyp=11
*Set Var ntri=6
*Set Var nnsurf=6
*endif
*Set Var mnumber=nnumber+1

/YD/YDC/MCSTEP *GenData(1,int)
/YD/YDC/NCSTEP 0 
/YD/YDC/ISAVE *GenData(14,int) 
/YD/YDC/DCGRAX *GenData(2)
/YD/YDC/DCGRAY *GenData(3)
/YD/YDC/DCGRAZ *GenData(4)
/YD/YDC/DCSIZC *GenData(8)
/YD/YDC/DCSIZF *GenData(9)
/YD/YDC/DCSIZS *GenData(11)
/YD/YDC/DCSIZV *GenData(10)
/YD/YDC/DCSTEC *GenData(5)
/YD/YDC/DCTIME 0.0
/YD/YDC/DCRMPT *GenData(17) 
/YD/YDC/DCURELX 1.0
/YD/YDC/INITER 2
/YD/YDC/ICOUTF *GenData(6)
/YD/YDC/ICOUTI 0
*if((strcmp(GenData(7),"64bit")==0))
/YD/YDC/ICOUTP 6
*endif
*if((strcmp(GenData(7),"32bit")==0))
/YD/YDC/ICOUTP 4
*endif
/YD/YDC/IWFAST 1   


*Set elems(tetrahedra)
*Set Var melem=nelem+1

/YD/YDE/MELEM  *melem  
/YD/YDE/NELEM  *nelem 
/YD/YDE/MELST   2  
/YD/YDE/NELST 2 
/YD/YDE/MELNO   *mnumber  
/YD/YDE/NELNO *nnumber 
/YD/YDE/D2ELST 21 0 0
/YD/YDE/D1EMCT 0
/YD/YDE/I1ELCF 0
/YD/YDE/I1ELBE 0
/YD/YDE/I1ELPR *nelem 
*loop elems
*ElemsMat 
*end elems

/YD/YDE/I2ELTO 21 *nnumber *nelem 
*loop elems
*ElemsConec *ElemsLayerNum
*end elems

/YD/YDE/D3TCS 21 0 0

/YD/YDI/MICOUP   *GenData(12)  
/YD/YDI/NICOUP 0 
/YD/YDI/IIECFF  -2
/YD/YDI/DIEDI   200.0
/YD/YDI/DIEZON  *GenData(13)
/YD/YDI/D1IESL 0
/YD/YDI/I1IECN 0
/YD/YDI/I1IECT 0
/YD/YDI/D2NV    21     3     0 
/YD/YDI/D2T1V    21     3     0 
/YD/YDI/D2T2V    21     3     0 
/YD/YDI/D1DELTAT1     0 
/YD/YDI/D1DELTAT2     0 
/YD/YDI/D1DELTAN    0 

*Set Var mnopo=npoin+1

/YD/YDN/MNODIM 4 
/YD/YDN/NNODIM 3
/YD/YDN/MNOPO *mnopo 
/YD/YDN/NNOPO  *npoin  
/YD/YDN/D2NCC 21 3 *npoin  
*loop nodes
*NodesCoord 
*end nodes 
/YD/YDN/D2NCI 21 3 *npoin  
*loop nodes
*NodesCoord 
*end nodes 
/YD/YDN/D2NFC 21  0 0 
/YD/YDN/D2NFT 21  0 0 
/YD/YDN/D1NMCT  0
/YD/YDN/D2NVC  21 0 0
/YD/YDN/I1NOBF *npoin  
*Set Cond FACENODES *nodes *NoCanRepeat
*loop nodes
*cond(1,int)   
*end nodes

/YD/YDN/I1NOPR  *npoin 0 
/YD/YDN/D1NTI 0

/YD/YDP/MPROP  200
/YD/YDP/NPROP  200
/YD/YDP/D1PEKS  0
/YD/YDP/D1PELA  0
/YD/YDP/D1PEMU  0
/YD/YDP/D1PEPE  0
/YD/YDP/D1PERO  0
/YD/YDP/D1PEGFN  0
/YD/YDP/D1PEGFS  0
/YD/YDP/D1PEFT  0
/YD/YDP/D1PCOH  0
/YD/YDP/D1PICF  0
/YD/YDP/D1PEFR  0
/YD/YDP/D1PESF  0
/YD/YDP/D1PEPSF  0
/YD/YDP/D1PEVF  0
/YD/YDP/D1PEPF  0
/YD/YDP/I1PTYP  0
/YD/YDP/I1PEJP  0
/YD/YDP/I1PEMN  0
/YD/YDP/I1PSDE  0
/YD/YDP/D1TCON  0
/YD/YDP/D1CAPA  0
/YD/YDP/D1CTEX  0 
/YD/YDB/MBCON 100 
/YD/YDB/NBCON 0 
/YD/YDB/D1BNAX  0 
/YD/YDB/D1BNAY  0 
/YD/YDB/D1BNAZ  0  
/YD/YDB/D1BNFX  0
/YD/YDB/D1BNFY  0
/YD/YDB/D1BNFZ  0 
/YD/YDB/D1BNVX  0
/YD/YDB/D1BNVY  0
/YD/YDB/D1BNVZ  0 
/YD/YDB/I1BNVX  0
/YD/YDB/I1BNVY  0
/YD/YDB/I1BNVZ  0
/YD/YDB/I1BNTP  0
/YD/YDB/D1BNTP  0
/YD/YDB/I1BNHF  0
/YD/YDB/D1BNHF  0
/YD/YDB/D1BCVT  0
/YD/YDB/D1BCVC  0
/YD/YDB/I1BCVT  0
Gen_out_only *GenData(16)

MATERIAL_LIST *nmats 
*loop materials
*MatProp(1) *MatProp(2) *MatProp(3) *MatProp(4) *MatProp(5) *MatProp(6) *MatProp(7) *MatProp(8) *MatProp(9) *MatProp(10) *MatProp(11) *MatProp(12,int) *MatProp(13) *MatProp(14) *MatProp(15,int) *MatProp(16) *MatProp(17) *MatProp(18) *MatProp(19) *MatProp(20) *eletyp 
*end materials

*set Cond Element  *elems *CanRepeat
*Set Var NElement=CondNumEntities
Element_List *NElement
*loop elems *OnlyInCond
*ElemsNum *cond(1) *cond(2,int)
*end elems

*set Cond Volume-Initial-Data  *nodes *NoRepeat
*Add Cond Surface-Initial-Data  *nodes *NoRepeat
*Add Cond Line-Initial-Data  *nodes *NoRepeat
*Add Cond Point-Initial-Data *nodes *NoRepeat
*Set Var NInitial=CondNumEntities

Initial_Data *NInitial
*loop nodes *OnlyInCond
*NodesNum *cond(1,int) *cond(2,real) *cond(3,real) *cond(4,real) *cond(5,int) *cond(6,real) *cond(7,real) *cond(8,real) *cond(9,int) *cond(10,real)
*end nodes

*Set Cond Volume-Constraints *nodes *NoRepeat
*Add Cond Surface-Constraints *nodes *NoRepeat
*Add Cond Line-Constraints *nodes *NoRepeat
*Add Cond Point-Constraints *nodes *NoRepeat
*Set Var NConstraints=CondNumEntities

Constraints_List *NConstraints
*loop nodes *OnlyInCond

*NodesNum *cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6) *cond(7) *cond(8) *cond(9) *cond(10) *cond(11) *cond(12) *cond(13) *cond(14) *cond(15) *cond(16)  *cond(17)  *cond(18)  *cond(19)  *cond(20)  *cond(21)   *cond(22)
*end nodes





/YD/YDO/MOHYS 1 /YD/YDO/NOHYS 1
/YD/YDO/DOHYP 0.0005
/YD/YDO/D1OHYC 1
1.0
/YD/YDO/D1OHYF 1
1.0 
/YD/YDO/D1OHYS 1
0
/YD/YDO/D1OHYT 1
0
/YD/YDO/D1OHYX 1
0
/YD/YDO/D1OHYY 1
0
/YD/YDO/I1OHYT 1
15

*Set elems(triangle)
/YD/YDX/NELEM *nelem
/YD/YDX/NELNO *nnsurf
/YD/YDX/DTHICK *GenData(18) 
/YD/YDX/I2ELTO 21 *ntri *nelem 
*loop elems
*ElemsConec
*end elems


*Set elems(all)
*Set Cond Surface-Load *elems *CanRepeat
*Set Var NConstraints=CondNumEntities
LOAD_LIST *NConstraints
*loop elems *OnlyInCond
*ElemsNum *globalnodes  *cond(1) 
*end elems
$YDOIT
$YSTOP
