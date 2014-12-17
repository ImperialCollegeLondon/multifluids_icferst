#include "confdefs.h"

subroutine test_shape_functions_p2_3d

  use fields
  use field_derivatives
  use state_module
  use spud
  use populate_state_module
  use vtk_interfaces
  use cv_advection, only: SCVDETNX_new
  use unittest_tools
  use shape_functions_prototype

  implicit none

  type(state_type), pointer, dimension(:) :: state
  type(mesh_type), pointer :: mesh

  integer :: indx
  integer, parameter :: CV_NGI_3d=10, SCVNGI=60, CV_SNLOC=6
  integer, parameter :: u_nloc=4, CV_NLOC=10
!  integer, parameter :: CV_NGI_3d=4, SCVNGI=18, CV_SNLOC=3
!  integer, parameter :: u_nloc=4, CV_NLOC=4
  integer, parameter :: nface=4, sbcvngi=6, u_snloc=3
  integer, parameter :: cv_ele_type=8

  logical:: CV_ON_FACE( CV_NLOC, SCVNGI )
  logical:: CVFEM_ON_FACE( CV_NLOC, SCVNGI )
  logical:: U_ON_FACE( U_NLOC, SCVNGI )
  logical:: UFEM_ON_FACE( U_NLOC, SCVNGI )
  integer, pointer :: NCOLGPTS
  integer, dimension(:), pointer ::  FINDGPTS, COLGPTS
  integer, dimension(:,:), pointer ::  CV_NEILOC, CV_SLOCLIST, U_SLOCLIST
  real, pointer, dimension( :, :, : ) :: SCVFENX_ALL
  real, pointer :: VOLUME, SVOLUME
  real, dimension(:), pointer :: detwei, ra, scvdetwei, scvra
  real, dimension(:,:,:), pointer :: inv_jac, scvinv_jac

  logical :: fail

  REAL, DIMENSION( : , : ), pointer :: CVN, CVN_SHORT, CVFEN, CVFEN_SHORT,&
      SCVFEN, SBCVFEN,  SBCVN, SBUFEN, SBUFENSLX, SBUFENSLY, SCVFENSLX,&
      SCVFENSLY, SUFEN, SUFENSLX, SUFENSLY, ufen, SBCVFENSLX,  SBCVFENSLY
  REAL, DIMENSION( : , :, : ), pointer :: CVFENLX_ALL, CVFENLX_SHORT_ALL,&
       SUFENLX_ALL, SBCVFENLX_ALL, SBUFENLX_ALL,&
       SCVFENLX_ALL, UFENLX_ALL, CVFENX_ALL
  real, dimension(:), pointer :: CVWEIGHT, CVWEIGHT_SHORT, SBCVFEWEIGH,&
       SCVFEWEIGH,  SELE_OVERLAP_SCALE

  real, parameter :: B=1875/9000.0d0, A=1.0d0-2*B
  real, dimension(cv_nloc,CV_NGI_3d) :: CVN_test

  type(vector_field), pointer :: positions

  real, dimension(3,cv_ngi_3d) :: x_all
 
  real, dimension(scvngi) :: l1,l2,l3,l4, normx,normy,normz, sarea, area
  real, dimension(cv_nloc) :: x_loc,Y_loc,Z_loc
  
  integer, dimension(cv_ngi_3d) :: dglno
  real, dimension(scvngi,3) :: normal

  integer :: i,j, par
  
  call load_options("data/shape_test_3d.mpml")
  call populate_state(state)

  positions=>extract_vector_field(state,"Coordinate")



  if (cv_nloc==10) then

  X_all(:,1)=positions%val(:,1)
  X_all(:,3)=positions%val(:,2)
  X_all(:,6)=positions%val(:,3)
  X_all(:,10)=positions%val(:,4)

  X_all(:,2)=0.5*(positions%val(:,1)+positions%val(:,2))
  X_all(:,4)=0.5*(positions%val(:,1)+positions%val(:,3))
  X_all(:,7)=0.5*(positions%val(:,1)+positions%val(:,4))
  X_all(:,5)=0.5*(positions%val(:,3)+positions%val(:,2))
  X_all(:,8)=0.5*(positions%val(:,4)+positions%val(:,2))
  X_all(:,9)=0.5*(positions%val(:,3)+positions%val(:,4))
  dglno=[1,2,3,4,5,6,7,8,9,10]
  else

     X_all(:,1:4)=positions%val(:,1:4)
     
     dglno(1:4)=[1,2,3,4]
  end if

  indx=0

  cvN_Test=0.0
 
  CALL cv_fem_shape_funs_plus_storage( &
       ! Volume shape functions...
       3, CV_ELE_TYPE,  &
       CV_NGI_3d, CV_NGI_3d, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
       CVWEIGHT, CVFEN, CVFENLX_ALL, &
       CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT_ALL, &
       UFEN, UFENLX_ALL, &
       ! Surface of each CV shape functions...
       SCVNGI, CV_NEILOC, CV_ON_FACE, CVFEM_ON_FACE, &
       SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
       SCVFENLX_ALL,  &
       SUFEN, SUFENSLX, SUFENSLY,  &
       SUFENLX_ALL,  &
       ! Surface element shape funcs...
       U_ON_FACE, UFEM_ON_FACE, NFACE, &
       SBCVNGI, SBCVN, SBCVFEN,SBCVFENSLX,&
       SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX_ALL, &
       SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX_ALL, &
       CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
       ! Define the gauss points that lie on the surface of the CV...
       FINDGPTS, COLGPTS, NCOLGPTS, &
       SELE_OVERLAP_SCALE, .false.,&
       STATE,'Test',indx )


  par=0
  CALL DETNLXR_INVJAC_PLUS_STORAGE(1, X_ALL,dglno,1,cv_ngi_3d, &
       CV_NLOC, SCVNGI, &
       SCVFEN, SCVFENLX_ALL, SCVFEWEIGH, SCVDETWEI, SCVRA, SVOLUME, .false. , &
       SCVFENX_ALL, &
       3, scVINV_JAC,state, "Surface_INVJAC", par)

  par=0
 CALL DETNLXR_INVJAC_PLUS_STORAGE(1, X_ALL,dglno,1,cv_ngi_3d, &
       CV_NLOC, CV_NGI_3d, &
       CVFEN, CVFENLX_ALL, CVWEIGHT, DETWEI, RA, VOLUME, .false. , &
       CVFENX_ALL, &
       3, INV_JAC,state, "INVJAC", par)

 print*, volume, sum(detwei)

area=area_from_position(positions%val)
call james_quadrature_quad_tet(l1, l2, l3,l4,  normx, normy, normz, sarea, &
                                        X_LOC, Y_LOC, Z_LOC, CV_NEILOC, cv_nloc, scvngi)

print*, '     Node #', '     quad pt', '              volume', '        total node volume'
  do i=1, size(CVN,1)
     do j=1, size(CVN,2)
        if (cvn(i,j)>0.0) &
             print*, i,j, detwei(j), sum(CVN(i,:)*detwei)
     end do
  end do

  do j=1, size(CVN,2)
             print*, 'x_quad',j,matmul(X_ALL(:,:),CVFEN(:,j))
  end do

print*, '     Node #', '   reference weight', '     ref weight/area'

  do j=1, size(scvdetwei)
     print*, j, scvdetwei(j), area(j), scvdetwei(j)/area(j)
  end do


do j=1,scvngi
  call SCVDETNX_new( 1,j,        &
       cv_nloc,     ScVNGI, 1,3,  &
       dglno,cv_nloc,&
       SCVDETWEI,normal,SCVFEN, SCVFENSLX,    &
       SCVFENSLY,  SCVFEWEIGH, X_ALL(:,1),     &
       X_ALL,        &
       .false.,       .true.,       .false. )
end do
  
print*, '     Node #', '         l1', '         l2', '           l3' , '      scvdetwei', '          area', '  scvdetwei/area', '     ref weight/area'

 do j=1, size(scvdetwei)
     print*, j, l1(j),l2(j),l3(j), scvdetwei(j), area(j), scvdetwei(j)/area(j)
  end do

  do j=1, cv_ngi_3d

     print*, j, sum(CVFENLX_ALL(1,:,j)*X_all(1,:))

  end do

do j=1, size(scvdetwei)
     print*, cv_neiloc(:,j)
  end do
  
  call report_test("[P1 Dual shape_functions]", fail, .false., "P1 shape function incorrect")

  call deallocate(state)

  call report_test_no_references()

contains

  function area_from_position(X) result(area)
    real, dimension(3,4) :: X
    real, dimension(3,10) :: p
    real, dimension(scvngi) :: area
    
    area=0.0

    p(:,[1,3,6,10])=X(:,:)
    p(:,2)=0.5*(X(:,1)+X(:,2))
    p(:,4)=0.5*(X(:,1)+X(:,3))
    p(:,5)=0.5*(X(:,2)+X(:,3))
    p(:,7)=0.5*(X(:,1)+X(:,4))
    p(:,8)=0.5*(X(:,2)+X(:,4))
    p(:,9)=0.5*(X(:,3)+X(:,4))

    !! face 1,3,6

    area(1)=face_area(p,[1,2,4])
    area(2)=face_area(p,[2,4,1])+face_area(p,[2,4,5])+face_area(p,[2,3,5])
    area(3)=face_area(p,[3,5,2])
    area(4)=face_area(p,[4,1,2])+face_area(p,[4,2,5])+face_area(p,[4,5,6])
    area(5)=face_area(p,[5,4,6])+face_area(p,[5,4,2])+face_area(p,[5,2,3])
    area(6)=face_area(p,[6,4,5])

    !! face 1,3,10

    area(7)=face_area(p,[1,2,7])
    area(8)=face_area(p,[2,7,1])+face_area(p,[2,7,8])+face_area(p,[2,3,8])
    area(9)=face_area(p,[3,8,2])
    area(10)=face_area(p,[7,1,2])+face_area(p,[7,2,8])+face_area(p,[7,8,10])
    area(11)=face_area(p,[8,7,10])+face_area(p,[8,7,2])+face_area(p,[8,2,3])
    area(12)=face_area(p,[10,7,8])


!! face 1,6,10

    area(13)=face_area(p,[1,4,7])
    area(14)=face_area(p,[4,7,1])+face_area(p,[4,7,9])+face_area(p,[4,6,9])
    area(15)=face_area(p,[6,9,4])
    area(16)=face_area(p,[7,1,4])+face_area(p,[7,4,9])+face_area(p,[7,9,10])
    area(17)=face_area(p,[9,7,10])+face_area(p,[9,7,4])+face_area(p,[9,4,6])
    area(18)=face_area(p,[10,7,9])

 !! face 3,6,10

    area(19)=face_area(p,[3,5,8])
    area(20)=face_area(p,[5,8,3])+face_area(p,[5,8,9])+face_area(p,[5,6,9])
    area(21)=face_area(p,[6,9,5])
    area(22)=face_area(p,[8,3,5])+face_area(p,[8,5,9])+face_area(p,[8,9,10])
    area(23)=face_area(p,[9,8,10])+face_area(p,[9,8,5])+face_area(p,[9,5,6])
    area(24)=face_area(p,[10,8,9])


  !! interior faces sub-tet 1,2,4,7
    
    area(25)=internal_area_tet(p,[1,2,4,7])
    area(26)=internal_area_tet(p,[1,4,2,7])
    area(27)=internal_area_tet(p,[1,7,2,4])
    area(28)=internal_area_tet(p,[2,4,7,1])
    area(29)=internal_area_tet(p,[2,7,1,4])
    area(30)=internal_area_tet(p,[4,7,1,2])

 !! interior faces sub-tet 2,3,5,8
    
    area(31)=internal_area_tet(p,[2,3,5,8])
    area(32)=internal_area_tet(p,[2,5,3,8])
    area(33)=internal_area_tet(p,[2,8,3,5])
    area(34)=internal_area_tet(p,[3,5,2,8])
    area(35)=internal_area_tet(p,[3,8,2,5])
    area(36)=internal_area_tet(p,[5,8,2,3])

  !! interior faces sub-tet 4,5,6,9
    
    area(37)=internal_area_tet(p,[4,5,6,9])
    area(38)=internal_area_tet(p,[4,6,5,9])
    area(39)=internal_area_tet(p,[4,9,5,6])
    area(40)=internal_area_tet(p,[5,6,4,9])
    area(41)=internal_area_tet(p,[5,9,4,6])
    area(42)=internal_area_tet(p,[6,9,4,5])

 !! interior faces sub-tet 7,8,9,10
    
    area(43)=internal_area_tet(p,[7,8,9,10])
    area(44)=internal_area_tet(p,[7,9,8,10])
    area(45)=internal_area_tet(p,[7,10,8,9])
    area(46)=internal_area_tet(p,[8,9,7,10])
    area(47)=internal_area_tet(p,[8,10,7,9])
    area(48)=internal_area_tet(p,[9,10,7,8])


  !! interior faces : octogon

    area(49)=internal_area_oct(p,[2,4,5,8,7,9])
    area(50)=internal_area_oct(p,[2,5,8,7,4,9])
    area(51)=internal_area_oct(p,[2,8,7,4,5,9])
    area(52)=internal_area_oct(p,[2,7,4,5,8,9])
    area(53)=internal_area_oct(p,[4,2,7,9,5,8])
    area(54)=internal_area_oct(p,[4,2,5,9,7,8])
    area(55)=internal_area_oct(p,[4,5,9,7,2,8])
    area(56)=internal_area_oct(p,[5,2,8,9,4,7])
    area(57)=internal_area_oct(p,[5,4,9,8,2,7])
    area(58)=internal_area_oct(p,[7,2,8,9,4,5])
    area(59)=internal_area_oct(p,[7,4,9,8,2,5])
    area(60)=internal_area_oct(p,[8,5,9,7,2,4])

    


  end function area_from_position
  
  function internal_area_oct(p,ids) result(area)
    real, dimension(:,:) :: p
    integer, dimension(:) :: ids
    real :: area

    area=quad_area(0.5*(p(:,ids(1))+p(:,ids(3))),&
         1.0/3.0*(p(:,ids(1))+p(:,ids(2))+p(:,ids(3))),&
         1.0/6.0*(p(:,ids(1))+p(:,ids(2))+p(:,ids(3))+p(:,ids(4))+p(:,ids(5))+p(:,ids(6))),&
         1.0/3.0*(p(:,ids(1))+p(:,ids(3))+p(:,ids(4))))
  end function internal_area_oct

  function internal_area_tet(p,ids) result(area)
    real, dimension(:,:) :: p
    integer, dimension(:) :: ids
    real :: area

    area=quad_area(0.5*(p(:,ids(1))+p(:,ids(2))),&
         1.0/3.0*(p(:,ids(1))+p(:,ids(2))+p(:,ids(3))),&
         1.0/4.0*(p(:,ids(1))+p(:,ids(2))+p(:,ids(3))+p(:,ids(4))),&
         1.0/3.0*(p(:,ids(1))+p(:,ids(2))+p(:,ids(4))))
  end function internal_area_tet

  function face_area(p,ids) result(area)
    real, dimension(:,:) :: p
    integer, dimension(:) :: ids
    real :: area
    
    area=quad_area(p(:,ids(1)),0.5*(p(:,ids(1))+p(:,ids(2))),&
         1.0/3.0*(p(:,ids(1))+p(:,ids(2))+p(:,ids(3))),&
         0.5*(p(:,ids(1))+p(:,ids(3))))

  end function face_area

  function cross4(v1,v2) result(vo)
          real, dimension(:), intent(in) :: v1,v2
          real, dimension(3) :: vo

          vo(1)=v1(2)*v2(3)-v2(2)*v1(3)
          vo(2)=v1(3)*v2(1)-v2(3)*v1(1)
          vo(3)=v1(1)*v2(2)-v2(1)*v1(2)          
        end function cross4

        function triangle_area(v1,v2) result(area)
          real, dimension(:), intent(in) :: v1,v2
          real :: area

          area=sqrt(sum(cross4(v1,v2)**2))/2.0
        end function triangle_area

        function quad_area(p1,p2,p3,p4) result(area)
          real, dimension(3), intent(in) :: p1,p2,p3,p4
          real :: area

          area=triangle_area(p3-p1,p2-p1)+&
               triangle_area(p3-p1,p4-p1)
        end function quad_area

end subroutine test_shape_functions_p2_3d
