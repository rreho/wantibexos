subroutine r_sort(arrin)
     !
     ! Sort real(dt) values from arrin into array
     ! arrout and give permutations in indx, indx_m1.
     ! Content of indx is destroyed.
     ! indx_m1: j went to   position indx_m1(j)
     ! indx   : i came from position indx(i)
     !
     double precision ::  arrin(:)
!     double precision, optional::  arrout(:)
!     integer,  optional::  indx(:)
!     integer,  optional::  indx_m1(:)
!     double precision, optional :: r_zero
     !
     ! local variables
     !
     integer  :: j, i,n, ir, l, indxt
     double precision :: q,r_zero_
     integer, allocatable:: l_indx(:)
     double precision,allocatable:: l_arrout(:)
     !
     r_zero_=0.0
!     if(present(r_zero)) r_zero_=r_zero
     !
     n=size(arrin)
     allocate(l_indx(n),l_arrout(n))
     !
     if(n.eq.1) then
       l_arrout(1) = arrin(1)
       l_indx(1) = 1
!       if (present(arrout)) arrout=l_arrout
!       if (.not.present(arrout)) 
       arrin=l_arrout
!       if (present(indx)) indx=l_indx
       deallocate(l_indx,l_arrout)
       return
     endif
     do j=1,n
       l_indx(j)=j
     enddo
     l=n/2+1
     ir=n
  1  continue
     if (l.gt.1)then
       l=l-1
       indxt=l_indx(l)
       q=arrin(indxt)
     else
       indxt=l_indx(ir)
       q=arrin(indxt)
       l_indx(ir)=l_indx(1)
       ir=ir-1
       if (ir.eq.1)then
         l_indx(1)=indxt
         go to 3
       endif
     endif
     i=l
     j=l+l
  2  if (j.le.ir)then
      if (j.lt.ir) then
        if (arrin(l_indx(j))<arrin(l_indx(j+1))-r_zero_) j=j+1
      endif
      if (q<arrin(l_indx(j))-r_zero_) then
        l_indx(i)=l_indx(j)
        i=j
        j=j+j
      else
        j=ir+1
      endif
      go to 2
     endif
     l_indx(i)=indxt
     go to 1
  3  continue
     do i=1,n
       l_arrout(i) = arrin(l_indx(i))
     enddo
!     if (present(arrout)) arrout=l_arrout
!     if (.not.present(arrout)) 
     arrin=l_arrout
!     if (present(indx)) indx=l_indx
!     if (present(indx_m1)) forall( i=1:n) indx_m1(l_indx(i))=i
     deallocate(l_indx,l_arrout)

end subroutine r_sort

subroutine rkmesh(rk,rlat,ngrid)

	implicit none
	
	integer,dimension(3) :: ngrid
	double precision :: rk
	double precision,dimension(3,3) :: rlat,blat
	double precision,dimension(3) :: vsize
	double precision,parameter :: pi=acos(-1.0)

	call recvec(rlat(1,:),rlat(2,:),rlat(3,:),blat(1,:),blat(2,:),blat(3,:))

	call vecsize(blat(1,:),vsize(1))
	call vecsize(blat(2,:),vsize(2))
	call vecsize(blat(3,:),vsize(3))

	vsize = vsize/(2*pi)

	ngrid(1) = int(max(1.0,(rk*vsize(1))+0.5))
	ngrid(2) = int(max(1.0,(rk*vsize(2))+0.5))
	ngrid(3) = int(max(1.0,(rk*vsize(3))+0.5))

end subroutine rkmesh

subroutine rkmesh2D(rk,rlat,ngrid)

	implicit none
	
	integer,dimension(3) :: ngrid
	double precision :: rk
	double precision,dimension(3,3) :: rlat,blat
	double precision,dimension(3) :: vsize
	double precision,parameter :: pi=acos(-1.0)

	call recvec(rlat(1,:),rlat(2,:),rlat(3,:),blat(1,:),blat(2,:),blat(3,:))

	call vecsize(blat(1,:),vsize(1))
	call vecsize(blat(2,:),vsize(2))
	call vecsize(blat(3,:),vsize(3))

	vsize = vsize/(2*pi)

	ngrid(1) = int(max(1.0,(rk*vsize(1))+0.5))
	ngrid(2) = int(max(1.0,(rk*vsize(2))+0.5))
	ngrid(3) = 1

end subroutine rkmesh2D

subroutine kpath2(outputfolder,rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	character(len=70) :: outputfolder    !pasta saida
	integer :: i,j,erro
	integer :: nks,npts
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension((nks/2)*npts,4) :: kpt

	double precision,dimension(nks,3) :: ks
	double precision,dimension(nks,3) :: ksaux

	double precision,dimension(nks/2) :: kdis !distancia entre os pontos k do caminho
	double precision :: kdistot !soma do caminho todo
	double precision :: kdisaux

	OPEN(UNIT=1733, FILE=trim(outputfolder)//'KLABELS.dat',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening KLABELS output file"

	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	do i=1,nks
	
		ksaux(i,1) = ks(i,1)*blat1(1)+ks(i,2)*blat2(1)+ks(i,3)*blat3(1)
		ksaux(i,2) = ks(i,1)*blat1(2)+ks(i,2)*blat2(2)+ks(i,3)*blat3(2)
		ksaux(i,3) = ks(i,1)*blat1(3)+ks(i,2)*blat2(3)+ks(i,3)*blat3(3)
	
	end do

	kdistot= 0.0
	do i=1,nks/2

		call distvec(ksaux(2*i,:),ksaux(2*i-1,:),kdis(i))

		!kdis(i)=sqrt((ksaux(2*i,1)-ksaux(2*i-1,1))**2+(ksaux(2*i,2)-ksaux(2*i-1,2))**2+(ksaux(2*i,3)-ksaux(2*i-1,3))**2)
		kdistot=kdistot+kdis(i)

	end do

	kdis=kdis/kdistot

	kdisaux=0.0
	write(1733,*) real(kdisaux)

	 do j=1,nks/2
		
	  do i=1,npts

		!kpt((j-1)*npts+i) = (kdisaux)+(kdis(j))*dble((i)/(npts))
		kpt((j-1)*npts+i,1) = (kdisaux)+(kdis(j))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,2) = ksaux(2*j-1,1)+ (ksaux(2*j,1)-ksaux(2*j-1,1))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,3) = ksaux(2*j-1,2)+ (ksaux(2*j,2)-ksaux(2*j-1,2))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,4) = ksaux(2*j-1,3)+ (ksaux(2*j,3)-ksaux(2*j-1,3))*(i-1.)/(npts-1.)	
		
		!write(1733,*) kpt((j-1)*npts+i,2),kpt((j-1)*npts+i,3),kpt((j-1)*npts+i,4)
		

	  end do

		kdisaux=kdisaux+kdis(j)
		write(1733,*) real(kdisaux)
	 end do


	close(1733)


end subroutine kpath2

subroutine kpathbse(outputfolder,rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	character(len=70) :: outputfolder    !pasta saida
	integer :: i,j,erro
	integer :: nks,npts
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension((nks/2)*npts,4) :: kpt

	double precision,dimension(nks,3) :: ks
	double precision,dimension(nks,3) :: ksaux

	double precision,dimension(nks/2) :: kdis !distancia entre os pontos k do caminho
	double precision :: kdistot !soma do caminho todo
	double precision :: kdisaux


	
	OPEN(UNIT=1733, FILE=trim(outputfolder)//'KLABELS-BSE.dat',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening KLABELS-BSE output file"


	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	do i=1,nks
	
		ksaux(i,1) = ks(i,1)*blat1(1)+ks(i,2)*blat2(1)+ks(i,3)*blat3(1)
		ksaux(i,2) = ks(i,1)*blat1(2)+ks(i,2)*blat2(2)+ks(i,3)*blat3(2)
		ksaux(i,3) = ks(i,1)*blat1(3)+ks(i,2)*blat2(3)+ks(i,3)*blat3(3)
	
	end do



	kdistot= 0.0
	do i=1,nks/2

		call distvec(ksaux(2*i,:),ksaux(2*i-1,:),kdis(i))

		!kdis(i)=sqrt((ksaux(2*i,1)-ksaux(2*i-1,1))**2+(ksaux(2*i,2)-ksaux(2*i-1,2))**2+(ksaux(2*i,3)-ksaux(2*i-1,3))**2)
		kdistot=kdistot+kdis(i)

	end do

	kdis=kdis/kdistot


	kdisaux=0.0
	write(1733,*) real(kdisaux)

	 do j=1,nks/2
		
	  do i=1,npts

		!kpt((j-1)*npts+i) = (kdisaux)+(kdis(j))*dble((i)/(npts))
		kpt((j-1)*npts+i,1) = (kdisaux)+(kdis(j))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,2) = ksaux(2*j-1,1)+ (ksaux(2*j,1)-ksaux(2*j-1,1))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,3) = ksaux(2*j-1,2)+ (ksaux(2*j,2)-ksaux(2*j-1,2))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,4) = ksaux(2*j-1,3)+ (ksaux(2*j,3)-ksaux(2*j-1,3))*(i-1.)/(npts-1.)	
		
		!write(1733,*) kpt((j-1)*npts+i,2),kpt((j-1)*npts+i,3),kpt((j-1)*npts+i,4)
		

	  end do

		kdisaux=kdisaux+kdis(j)
		write(1733,*) real(kdisaux)
	 end do


	close(1733)


end subroutine kpathbse

subroutine kpath(outputfolder,rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	character(len=70) :: outputfolder    !pasta saida
	integer :: i,j,erro
	integer :: nks,npts
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension((nks-1)*npts,4) :: kpt

	double precision,dimension(nks,3) :: ks
	double precision,dimension(nks,3) :: ksaux

	double precision,dimension(nks-1) :: kdis !distancia entre os pontos k do caminho
	double precision :: kdistot !soma do caminho todo
	double precision :: kdisaux

	OPEN(UNIT=1733, FILE=trim(outputfolder)//'KLABELS-BSE.dat',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening KLABELS-BSE output file"

	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	do i=1,nks
	
		ksaux(i,1) = ks(i,1)*blat1(1)+ks(i,2)*blat2(1)+ks(i,3)*blat3(1)
		ksaux(i,2) = ks(i,1)*blat1(2)+ks(i,2)*blat2(2)+ks(i,3)*blat3(2)
		ksaux(i,3) = ks(i,1)*blat1(3)+ks(i,2)*blat2(3)+ks(i,3)*blat3(3)
	
	end do


	kdistot= 0.0
	do i=1,nks-1

		call distvec(ksaux(i,:),ksaux(i+1,:),kdis(i))
		kdistot=kdistot+kdis(i)
	end do

	kdis=kdis/kdistot


	do i=1,npts

		kpt(i,1) = (kdis(1))*(i-1.)/(npts-1.)
		kpt(i,2) = ksaux(1,1)+ (ksaux(2,1)-ksaux(1,1))*(i-1.)/(npts-1.)
		kpt(i,3) = ksaux(1,2)+ (ksaux(2,2)-ksaux(1,2))*(i-1.)/(npts-1.)
		kpt(i,4) = ksaux(1,3)+ (ksaux(2,3)-ksaux(1,3))*(i-1.)/(npts-1.)

	end do

	if (nks .gt. 2) then

	kdisaux=0.0
	write(1733,*) real(kdisaux)
	 do j=2,nks-1
		kdisaux=kdisaux+kdis(j-1)
		write(1733,*) real(kdisaux)
	  do i=1,npts

		kpt((j-1)*npts+i,1) = (kdisaux)+(kdis(j))*(i)/(npts)
		kpt((j-1)*npts+i,2) = ksaux(j,1)+ (ksaux(j+1,1)-ksaux(j,1))*(i)/(npts)
		kpt((j-1)*npts+i,3) = ksaux(j,2)+ (ksaux(j+1,2)-ksaux(j,2))*(i)/(npts)
		kpt((j-1)*npts+i,4) = ksaux(j,3)+ (ksaux(j+1,3)-ksaux(j,3))*(i)/(npts)
		
		!write(1733,*) kpt((j-1)*npts+i,2),kpt((j-1)*npts+i,3),kpt((j-1)*npts+i,4)	

	  end do
	 end do

	else 

	write(1733,*) real(0.0000)
	write(1733,*) real(1.0000)
	 continue

	end if


end subroutine kpath

subroutine distvec(vec1,vec2,distvecout) !subroutina para calcular a distância entre dois vetores

	implicit none

	double precision,dimension(3) :: vec1,vec2
	double precision :: distvecout

	distvecout=sqrt((vec1(1)-vec2(1))**2+(vec1(2)-vec2(2))**2+(vec1(3)-vec2(3))**2)


end subroutine distvec


subroutine monhkhorst_packq(q,n1,n2,n3,shift,rlat1,rlat2,rlat3,qpt)

	implicit none
	
	integer :: n1,n2,n3
	double precision,dimension(3) :: shift,kshift
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(n1*n2*n3,3) :: qpt

	double precision,dimension(4) :: q

	integer :: i,j,k, counter
	double precision,dimension(3) :: blat1,blat2,blat3

	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	kshift(1) = blat1(1)*shift(1)+blat2(1)*shift(2)+blat3(1)*shift(3) 
	kshift(2) = blat1(2)*shift(1)+blat2(2)*shift(2)+blat3(2)*shift(3)
	kshift(3) = blat1(3)*shift(1)+blat2(3)*shift(2)+blat3(3)*shift(3)

	counter = 1

	do i=0,n1-1
	 do j=0,n2-1
	  do k=0,n3-1

	    !qpt(counter,1) = q(2)+(blat1(1)/dble(n1))*(dble(i)+shift(1))+(blat2(1)/dble(n2))*(dble(j)+shift(2))&
		!	     +(blat3(1)/dble(n3))*(dble(k)+shift(3))

	    !qpt(counter,2) = q(3)+(blat1(2)/dble(n1))*(dble(i)+shift(1))+(blat2(2)/dble(n2))*(dble(j)+shift(2))&
	!		     +(blat3(2)/dble(n3))*(dble(k)+shift(3))

	    !qpt(counter,3) =q(4)+(blat1(3)/dble(n1))*(dble(i)+shift(1))+(blat2(3)/dble(n2))*(dble(j)+shift(2))&
	!		     +(blat3(3)/dble(n3))*(dble(k)+shift(3))	


	    qpt(counter,1) = q(2)+(blat1(1)/dble(n1))*(dble(i))+(blat2(1)/dble(n2))*(dble(j))&
			     +(blat3(1)/dble(n3))*(dble(k))+kshift(1)

	    qpt(counter,2) = q(3)+(blat1(2)/dble(n1))*(dble(i))+(blat2(2)/dble(n2))*(dble(j))&
			     +(blat3(2)/dble(n3))*(dble(k))+kshift(2)

	    qpt(counter,3) =q(4)+(blat1(3)/dble(n1))*(dble(i))+(blat2(3)/dble(n2))*(dble(j))&
			     +(blat3(3)/dble(n3))*(dble(k))+kshift(3)	
   
	    counter = counter+1

	  end do
	 end do
	end do

end subroutine monhkhorst_packq


subroutine monhkhorst_pack(n1,n2,n3,shift,rlat1,rlat2,rlat3,kpt)

	implicit none
	
	integer :: n1,n2,n3
	double precision,dimension(3) :: shift,kshift
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(n1*n2*n3,3) :: kpt

	integer :: i,j,k, counter
	double precision,dimension(3) :: blat1,blat2,blat3

	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	counter = 1

	kshift(1) = blat1(1)*shift(1)+blat2(1)*shift(2)+blat3(1)*shift(3) 
	kshift(2) = blat1(2)*shift(1)+blat2(2)*shift(2)+blat3(2)*shift(3)
	kshift(3) = blat1(3)*shift(1)+blat2(3)*shift(2)+blat3(3)*shift(3)

	do i=0,n1-1
	 do j=0,n2-1
	  do k=0,n3-1

	    !kpt(counter,1) = (blat1(1)/dble(n1))*(dble(i)+shift(1))+(blat2(1)/dble(n2))*(dble(j)+shift(2))&
		!	     +(blat3(1)/dble(n3))*(dble(k)+shift(3))

	    !kpt(counter,2) = (blat1(2)/dble(n1))*(dble(i)+shift(1))+(blat2(2)/dble(n2))*(dble(j)+shift(2))&
		!	     +(blat3(2)/dble(n3))*(dble(k)+shift(3))

	    !kpt(counter,3) =(blat1(3)/dble(n1))*(dble(i)+shift(1))+(blat2(3)/dble(n2))*(dble(j)+shift(2))&
		!	     +(blat3(3)/dble(n3))*(dble(k)+shift(3))	

	    kpt(counter,1) = (blat1(1)/dble(n1))*(dble(i))+(blat2(1)/dble(n2))*(dble(j))&
			     +(blat3(1)/dble(n3))*(dble(k))+kshift(1)

	    kpt(counter,2) = (blat1(2)/dble(n1))*(dble(i))+(blat2(2)/dble(n2))*(dble(j))&
			     +(blat3(2)/dble(n3))*(dble(k))+kshift(2)

	    kpt(counter,3) =(blat1(3)/dble(n1))*(dble(i))+(blat2(3)/dble(n2))*(dble(j))&
			     +(blat3(3)/dble(n3))*(dble(k))+kshift(3)	

	    counter = counter+1

	  end do
	 end do
	end do

end subroutine monhkhorst_pack


subroutine get_ijk(kpoint,ijk,nkpt,rlat)
  implicit none

  integer, intent(in) :: nkpt(3)
  integer, intent(out) :: ijk(nkpt(1)*nkpt(2)*nkpt(3),3)
  double precision, intent(in) :: kpoint(nkpt(1)*nkpt(2)*nkpt(3),3), rlat(3,3)
  
  integer :: id,n

  do id = 1, 3
    do n = 1, nkpt(1)*nkpt(2)*nkpt(3)
      ijk(n,id) = int(nkpt(id)*dot_product(rlat(id,1:3),kpoint(n,1:3)))
    enddo
  enddo

end subroutine get_ijk

subroutine find_neighbour(nkpt,dlat,kpt,neighbours,dk_red,dk_car)
  implicit none
  integer, intent(in) :: nkpt(3)
  double precision, intent(in) :: dlat(3,3), kpt(nkpt(1)*nkpt(2)*nkpt(3),3)
  double precision :: x(nkpt(1),nkpt(2)), y(nkpt(1),nkpt(2)), z(nkpt(1),nkpt(2),nkpt(3)), tmp_x, tmp_y, tmp_z
  double precision :: r_h(3), q(3),qip1_j(3), qip1_jp1(3), qi_jp1(3), phi_n, phase, Nkpoints(3), &
              nx,ny,nz, sum_phases,delta, deltakkp, k_red(nkpt(1)*nkpt(2)*nkpt(3),3), dk_red(3), dk_car(3)

  integer :: i, j, k, l, n, counter, m, di, dj, dk, ik, nmax, id, i_q, i_qip1_j, i_qip1_jp1, i_qjp1    
  double precision :: blat(3,3)
  integer, intent(inout) :: neighbours(nkpt(1)*nkpt(2)*nkpt(3),6,2)
  
  call recvec(dlat(1,:),dlat(2,:),dlat(3,:),blat(1,:),blat(2,:),blat(3,:))   

  neighbours(:,:,2) = 1

  ! 1st step is to convert all k-points to reduced coordinates
  nmax = nkpt(1)*nkpt(2)*nkpt(3)
  do n = 1, nmax
    k_red(n,1:3) = matmul(dlat(1:3,1:3),kpt(n,1:3))
  enddo

  ! Shift all points with coordinates that are outside the cube ]-0.5,0.5]^3 to inside the cube by adding or subtracting 1
  do ik = 1, nmax
    do id  = 1, 3
      if ((k_red(ik,id) .le. -0.5*(1.-0.0005))) then    ! check flag for kpoint tolerance 
        k_red(ik,id) = k_red(ik,id) + 1.
      elseif ((k_red(ik,id) .gt. 0.5*(1.+0.0005))) then  ! same as above
        k_red(ik,id) = k_red(ik,id) - 1.
      endif
    enddo
  enddo

  ! call r_sort to sort the points along the 1st dimension  
!  call r_sort(k_red(:,1))
  ! but they must also be sorted along other dimensions 
  ! I stick with the 2D case for now, as I am unsure if everything works
  ! this should be equivalent to having k=0 in the listing of k-points
  do i = 1, nkpt(1)
    do j = 1, nkpt(2)
      x(i,j) = k_red(1 + nkpt(3)*j+nkpt(3)*nkpt(2)*(i-1),1)
      y(i,j) = k_red(1 + nkpt(3)*j+nkpt(3)*nkpt(2)*(i-1),2)
    enddo
  enddo

  do j = 1, nkpt(2)
    do i = 1, nkpt(1)
      do l = 1, nkpt(2)
        do k = 1, nkpt(1)
          if (y(i,j) < y(k,l)) then
            tmp_y = y(i,j)
            y(i,j) = y(k,l)
            y(k,l) = tmp_y
           end if
         end do
      end do
    end do
  end do

  dk_red = (/1./nkpt(1),1./nkpt(2),1./nkpt(3)/)
      
  dk_car = matmul(blat,dk_red)

! Following Riccardo's implementation, this should mean that we can do the loop using reduced coordinates in the following way
!!(i,j+1) <----------(i+1,j+1)
!!        |         |
!!        |         |
!!        |         |
!!  (i,j) ---------->(i+1,j)

  ! Calculate the flux through each plaquette
! total_flux =0.00
! all_phases =0.00
  do k = 1, nkpt(3) 
    do j = 1, nkpt(2) 
      do i = 1, nkpt(1) 
        q        = (/x(i,j), y(i,j),dble(0.0)/)
        qip1_j   = (/x(i+1,j),y(i+1,j),dble(0.0)/)
        qip1_jp1 = (/x(i+1,j+1), y(i+1,j+1),dble(0.0)/)
        qi_jp1   = (/x(i,j+1), y(i,j+1), dble(0.0)/)
        i_q        = 1 + nkpt(3)*j + nkpt(3)*nkpt(2)*i
        i_qip1_j   = 1 + nkpt(3)*j + nkpt(3)*nkpt(2)*(i+1)
        i_qip1_jp1 = 1 + nkpt(3)*(j+1) + nkpt(3)*nkpt(2)*(i+1)
        i_qjp1     = 1 + nkpt(3)*(j+1) + nkpt(3)*nkpt(2)*i

        neighbours(n,1,1) = 1 + k     + nkpt(3)*j      + nkpt(3)*nkpt(2)*(i+di) ! neighbour in +dx
        neighbours(n,2,1) = 1 + k     + nkpt(3)*j      + nkpt(3)*nkpt(2)*(i-di) ! neighbour in -dx
        neighbours(n,3,1) = 1 + k     + nkpt(3)*(j+dj) + nkpt(3)*nkpt(2)*i      ! neighbour in +dy
        neighbours(n,4,1) = 1 + k     + nkpt(3)*(j-dj) + nkpt(3)*nkpt(2)*i      ! neighbour in -dx
        neighbours(n,5,1) = 1 + k + 1 + nkpt(3)*j      + nkpt(3)*nkpt(2)*i      ! neighbour in +dz
        neighbours(n,6,1) = 1 + k - 1 + nkpt(3)*j      + nkpt(3)*nkpt(2)*i      ! neighbour in -dx
!      overlap = S_lambda_qqp(i1,i_qip1_j,i_q)* &
!               S_lambda_qqp(i1,i_qip1_jp1,i_qip1_j)* &
!               S_lambda_qqp(i1,i_q,i_qjp1)
!      write (*,*) 'overlap', overlap
!      flux(i,j) = (x(i+1,j) - x(i,j))*(y(i,j+1) - y(i,j)) * overlap - &
!                  (x(i,j+1) - x(i,j))*(y(i+1,j) - y(i,j)) * overlap
!      phi_n = aimag(log(overlap))!atan2(aimag(overlap),real(overlap))
!      write (*,*) 'phi_n', phi_n
!      all_phases(i,j) = phi_n
!      total_flux = total_flux + flux(i,j)
!      sum_phases = sum_phases + all_phases(i,j)/(Pi)
      enddo
    end do
  end do



!   i = int(n1*(rlat1(1)*kpt(n,1) + rlat2(1)*kpt(n,2) + rlat3(1)*kpt(n,3)))
!   j = int(n2*(rlat1(2)*kpt(n,1) + rlat2(2)*kpt(n,2) + rlat3(2)*kpt(n,3)))
!   k = int(n3*(rlat1(3)*kpt(n,1) + rlat2(3)*kpt(n,2) + rlat3(3)*kpt(n,3)))

!   di = int(rlat1(1)*blat1(1)/n1 + rlat2(1)*blat2(1)/n2 + rlat3(1)*blat3(1)/n3)
!   dj = int(rlat1(2)*blat1(1)/n1 + rlat2(2)*blat2(1)/n2 + rlat3(2)*blat3(1)/n3)
!   dk = int(rlat1(3)*blat1(1)/n1 + rlat2(3)*blat2(1)/n2 + rlat3(3)*blat3(1)/n3)

!   ! no, these are the neighbours along b1, b2, and b3, the reciprocal lattice vectors, not x,y,z in the reciprocal lattice
!   neighbours(n,1,1) = 1 + k    + n3*j      + n3*n2*(i+di) ! neighbour in +dx
!   neighbours(n,2,1) = 1 + k    + n3*j      + n3*n2*(i-di) ! neighbour in -dx
!   neighbours(n,3,1) = 1 + k    + n3*(j+dj) + n3*n2*i      ! neighbour in +dy
!   neighbours(n,4,1) = 1 + k    + n3*(j-dj) + n3*n2*i      ! neighbour in -dx
!   neighbours(n,5,1) = 1 + k+dk + n3*j      + n3*n2*i      ! neighbour in +dz
!   neighbours(n,6,1) = 1 + k-dk + n3*j      + n3*n2*i      ! neighbour in -dx

!   do m = 1, 6
!     if (neighbours(n,m,1) > n1*n2*n3 .or. neighbours(n,m,1) < 0) neighbours(n,m,2) = 0
!   enddo


end subroutine find_neighbour


subroutine recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3) !calcula os vetores da rede recíproca, a partir dos vetores da rede real

	implicit none
	
	double precision,parameter :: pi=acos(-1.)
	
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension(3) :: v23,v31,v12
	double precision :: vol

	call prodvec(rlat2,rlat3,v23)
	call prodvec(rlat3,rlat1,v31)
	call prodvec(rlat1,rlat2,v12)

	vol= abs((rlat1(1)*v23(1))+(rlat1(2)*v23(2))+(rlat1(3)*v23(3)))

	blat1 = ((2.0*pi)/vol)*v23
	blat2 = ((2.0*pi)/vol)*v31
	blat3 = ((2.0*pi)/vol)*v12


end subroutine recvec



!subrotinas para gerar grid hexagonal

subroutine gridgenhexcount(ngrid,a,counter)


	implicit none

	double precision,parameter :: pi=acos(-1.)

	integer :: ngrid,i,j

	integer :: counter

	double precision :: a!constante da rede

	double precision :: a0

        double precision :: xmax, ymax

	double precision :: kx,ky


	a0 = a/sqrt(3.)

        xmax=2*pi/(3.d0*a0)
        ymax=4*pi/(3.d0*a)


	counter = 0


	do i=1,ngrid
                kx=xmax*(i-1)/(ngrid-1)
		do j=1,ngrid
                         ky=ymax*(j-1)/(ngrid-1)
                         if (ky.lt.ymax/2.d0) then
                            counter=counter+1
                            if (i.ne.1.and.j.ne.1) counter=counter+1 
                            if (i.ne.1) counter=counter+1
	                    if (j.ne.1) counter=counter+1
                          
                         else if (ky.lt.(ymax-kx*tan(pi/6.d0))) then
                            counter=counter+1
                            if (i.ne.1.and.j.ne.1) counter=counter+1 
                            if (i.ne.1) counter=counter+1
	                    if (j.ne.1) counter=counter+1
                         end if
		
		end do
	end do 


end subroutine gridgenhexcount


subroutine gridgenhex(ngrid,a,ncounter,kg)


	implicit none

	double precision,parameter :: pi=acos(-1.)

	integer :: ngrid,i,j

	integer :: counter,ncounter

	double precision :: a!constante da rede

	double precision :: a0

        double precision :: xmax, ymax

	double precision :: kx,ky

	double precision,dimension(ncounter,2) :: kg


	a0 = a/sqrt(3.)

        xmax=2*pi/(3.d0*a0)
        ymax=4*pi/(3.d0*a)

	counter = 0

	do i=1,ngrid
                kx=xmax*(i-1)/(ngrid-1)
		do j=1,ngrid
                         ky=ymax*(j-1)/(ngrid-1)
                         if (ky.lt.ymax/2.d0) then
                            counter=counter+1 
                            kg(counter, 2)=kx
                            kg(counter, 1)=ky
                            if (i.ne.1.and.j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=-ky
                            end if
                            if (i.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=ky
			    end if
	                    if (j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=kx
                                kg(counter, 1)=-ky
                            end if
                         else if (ky.lt.(ymax-kx*tan(pi/6.d0))) then
                            counter=counter+1 
                            kg(counter, 2)=kx
                            kg(counter, 1)=ky
                            if (i.ne.1.and.j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=-ky
                            end if
                            if (i.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=ky
			    end if
	                    if (j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=kx
                                kg(counter, 1)=-ky
                            end if
                         end if
                
		
		end do
        end do

	


end subroutine gridgenhex
