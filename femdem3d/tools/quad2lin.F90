
!this program does the task of extracting the corner node connectivity and coordinates 
!from a quadratic tetrahedral mesh, and outputting it together with an array of pointers which 
!indicate which node in the linear mesh corresponds to the node in the quadratic mesh
program quad2lin
  implicit none
  integer :: nnodes, nelements,linearnodes
  integer :: dat_unit,inod,ielem,i,j,iloc,node,node1
  integer, allocatable, dimension(:) :: ndglno,iread_buffer,quad2lin1
  real, allocatable :: read_buffer(:),coordinates(:,:),lcoordinates(:,:)
  character(len=4096) buffer,filename
  
  write(*,*) 'Please type input filename (without extension)'
  write(*,*) '(make sure extension of filename is .dat)'
  read(*,*) filename
  ! Open node file
  dat_unit=123
  open(unit=dat_unit, file=trim(filename)//'.dat', status='old')  
  read(dat_unit, *) buffer
  ! Read the number of nodes, and elements
  read(dat_unit, *) nnodes, nelements  
  write(*,*) nnodes, nelements
  read(dat_unit, *) buffer ! Skip line
  allocate(read_buffer(4))
  allocate(coordinates(nnodes,4))
  allocate(lcoordinates(nnodes,4))
  do inod=1,nnodes
     read(dat_unit,*) read_buffer
     forall (j=1:4)
        coordinates(inod,j)=read_buffer(j)
     end forall
  end do
  deallocate(read_buffer)

  allocate(iread_buffer(5))
  allocate(ndglno(4*nelements))
  read(dat_unit, *) buffer  ! Skip line
  !Read corner nodes
  do ielem=1,nelements  !
     read(dat_unit,*) iread_buffer
     forall (j=1:4)
        ndglno((ielem-1)*4+j)=iread_buffer(j+1)
     end forall
  end do
  deallocate(iread_buffer)
  close(dat_unit)

  allocate(quad2lin1(nnodes))
  quad2lin1=0
  linearnodes=0
  nodeloop: do inod=1,nnodes
     elementloop: do ielem=1,nelements
        do iloc=1,4
           node=ndglno(4*(ielem-1)+iloc)
           if (node==inod) then
              linearnodes=linearnodes+1
              do j=1,4
                 lcoordinates(linearnodes,j)=coordinates(inod,j)                 
              end do
              quad2lin1(inod)=linearnodes
              exit elementloop
           end if
        end do
     end do elementloop
  end do nodeloop

  do inod=1,linearnodes
     node1=lcoordinates(inod,1) !node to be replaced
     do ielem=1,nelements
        do iloc=1,4
           node=ndglno(4*(ielem-1)+iloc)
           node1=lcoordinates(inod,1)
           if(node==node1) then
              ndglno(4*(ielem-1)+iloc)=inod !replaced node
           end if
        end do
     end do
  end do

  open(unit=dat_unit, file=trim(filename)//'-2.dat', status='replace')  
  write(dat_unit,'(A)') 'Nodes   Elems'
  write(dat_unit,*) linearnodes,nelements,' 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
  write(dat_unit,'(A)') 'node    x      y     z'
  do inod=1,linearnodes
     write(dat_unit,'(I0,3F12.6)') inod,lcoordinates(inod,2),lcoordinates(inod,3),lcoordinates(inod,4)
  end do
  write(dat_unit,'(A)') 'Connectivity'
  do ielem=1,nelements 
     write(dat_unit,'(A3,I0,A1,I0,A1,I0,A1,I0)')' 4 ',ndglno((ielem-1)*4+1),' ',ndglno((ielem-1)*4+2),' ',&
          ndglno((ielem-1)*4+3),' ',ndglno((ielem-1)*4+4)
  end do
  do i=1,20
     write(dat_unit,'(A)')'Face Element nodes for 01'
  end do
  close(dat_unit)
  open(unit=dat_unit, file='quad2lin.dat', status='replace')  
  write(dat_unit,'(A)') 'Quadratic Number of Nodes'
  write(dat_unit,'(I0)') nnodes
  write(dat_unit,'(A)') 'linear to quadratic array'
  do inod=1,linearnodes
     write(dat_unit,'(I0)') int(lcoordinates(inod,1))-1
  end do
  write(dat_unit,'(A)') 'quadratic to linear array'
  do inod=1,nnodes
     write(dat_unit,'(I0)') quad2lin1(inod)-1
  end do
  close(dat_unit)

  
end program quad2lin


