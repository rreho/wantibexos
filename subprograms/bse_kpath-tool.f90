!gfortran -mcmodel=large bse_kpath-tool.f90 -o bse_kpath-tool.x -llapack95 -lopenblas -fopenmp 

subroutine bsebnds(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
         ebse0,ebsef,numbse,sme,ktol,params,params_mmn,kpaths,kpathsbse,orbw,ediel, &
         exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,berryexc,excwf0,excwff)

  use omp_lib
  use hamiltonian_input_variables
  use bse_types

  implicit none

  double precision,parameter:: pi=acos(-1.)
  integer :: dimbse, nbands,nntot
  double precision,allocatable,dimension(:,:) :: exk !(dimbse,nqpts*npath)
  double complex, allocatable, dimension(:,:,:) :: wfk

  integer :: ncaux,nvaux, index

  integer,allocatable,dimension(:,:) :: stt !(ngrid*ngrid*nc*nv,4)
  double precision,allocatable,dimension(:,:) :: kpt,qpt, kpt_red, qpt_red !pontos k do grid (ngrid*ngrid,2)

  double precision,dimension(4) :: q

  double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
  double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

  double precision,allocatable,dimension(:,:) :: energy !variavel que guarda os autovalores para cada ponto k
  double complex,allocatable,dimension(:,:,:) :: vector !variavel que guarda os autovetores para cada ponto k

  double precision,allocatable,dimension(:,:) :: energyq !variavel que guarda os autovalores para cada ponto k
  double complex,allocatable,dimension(:,:,:) :: vectorq !variavel que guarda os autovetores para cada ponto k

  double precision,allocatable,dimension(:,:) :: qauxv !(nqpts*npath,4) 

  integer:: counter,c,v,i,j,f,l,h,erro,i2,k,ip,it

  double complex:: matrizelbsekq

  double precision :: a,r0,ed
  real, parameter :: tol = 1.0e-6 
  integer :: ngkpt

  !variaveis relacionadas a marcacao do tempo

  double precision:: t0,tf
  integer,dimension(8) :: values,values2

  !variaveis kpoints

  integer :: nks
  integer :: nkpts, nkpts_mmn
  double precision,allocatable,dimension(:,:) :: ks

  integer,allocatable,dimension(:) :: nocpk,nocpq


  !definicoes diagonalizacao 
  INTEGER   ::       ifail
  double precision,parameter :: ABSTOL=1.0e-6
  INTEGER          INFO
  double precision,allocatable,dimension(:) :: W,RWORK
  COMPLEX*16,allocatable,dimension(:,:) :: hbse

  INTEGER ::         LWMAX
  INTEGER ::         LWORK
  INTEGER ::         LIWORK, LRWORK
  INTEGER,allocatable,dimension(:) :: IWORK
        double complex,allocatable,dimension (:) :: WORK

  !modificacoes versao 2.1

  integer :: nthreads
  integer,dimension(3) :: ngrid
  integer :: nc,nv
  !integer :: ncrpa,nvrpa
  !integer :: ncbz,nvbz
  double precision :: edos0,edosf,numdos
  double precision :: ebse0,ebsef,numbse
  double precision :: sme,exc,rk
  double precision,dimension(3) :: mshift
  double precision :: ktol
  character(len=70) :: params   !parametros TB
  character(len=70) :: params_mmn !parameter mmn
  character(len=70) :: orbw     !peso orbitais lcount
  character(len=70) :: kpaths    !kpath
  character(len=70) :: kpathsbse    !kpath
  !character(len=70) :: diein    !ambiente dieletrico
  character(len=70) :: outputfolder    !pasta saida
  character(len=70) :: calcparms
  character(len=70) :: meshtype
  character(len=5) :: coultype
  character(len=30) :: Format
  double precision,dimension(3) :: ediel
  double precision :: ez,w1,lc
  logical :: bsewf
  logical :: berryexc
  integer :: excwf0,excwff  
  type(bse_coeff) :: bse_coefficient !
  complex, allocatable    :: A_table(:,:,:,:,:), bse_table(:,:,:), exc_exc_overlap(:,:,:)
  double complex, allocatable, dimension(:,:,:,:) :: Mmn ! overlap matrix
  double complex, allocatable, dimension(:,:,:,:) :: o_mmn

  complex, parameter :: czero = (0.0,0.0)
  complex, parameter :: cone = (1.0,0.0)
  complex, parameter :: imun = (0.0,1.0)
  double precision, parameter :: half = 0.5

  integer :: nqpts, dir1, dir2, dir3, ic, iq, n, id1, id2, id3, iq_d1, iq_d2, iq_d3 
  integer :: ijk(ngrid(1)*ngrid(2)*ngrid(3),3), ic1, ic2,iv, iv1, iv2, ik, id, ii
  integer :: kmq_index(ngrid(1)*ngrid(2)*ngrid(3),ngrid(1)*ngrid(2)*ngrid(3)), kmq_ijk(3), i_kmq
  double precision, allocatable :: kmq_grid(:,:,:)  ! ik-iq point =  kmq_grid(ik,iq,1:3)
  integer, allocatable :: q_neigh(:,:)     ! 
  complex, allocatable :: chern_exc_exc(:,:) ! exciton-exciton contribution to the total chern number, chern_exc_exc(dimbse,3)
  complex, allocatable :: chern_exc_hole(:,:) ! exciton-hole contribution to the total chern number, chern_exc_hole(dimbse,3)
  complex, allocatable :: chern_exc_electron(:,:) ! exciton-electron contribution to the total chern number, chern_exc_electron(dimbse,3)
  integer, allocatable :: kmq_neighbours(:,:,:) ! neighbours in the k-q grid
  complex ::  rot_overlap, tmp_o(3)
  double precision :: blat(3,3), dk_red(3), dk_car(3)

  interface
  function find_kpoint_index(kpoints, nkpts, point) result(index)
      integer, intent(in) :: nkpts
      double precision, dimension(nkpts, 3), intent(in) :: kpoints
      double precision, dimension(3), intent(in) :: point
      integer :: index
  end function find_kpoint_index
  end interface  
  !call input_read
  ! INPUT 
  OPEN(UNIT=500, FILE= kpathsbse,STATUS='old', IOSTAT=erro)
  if (erro/=0) stop "Error opening bse-kpath input file"
  !OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
  !if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"

  !OUTPUT : criando arquivos de saida
  OPEN(UNIT=300, FILE=trim(outputfolder)//"log_bse_kpath.dat",STATUS='unknown', IOSTAT=erro)
  if (erro/=0) stop "Error opening log_bse_kpath output file"

  OPEN(UNIT=400, FILE=trim(outputfolder)//"bands_bse.dat",STATUS='unknown', IOSTAT=erro)
  if (erro/=0) stop "Error opening bands_bse output file"

  OPEN(UNIT=600, FILE=trim(outputfolder)//"bands_wf_bse.dat",STATUS='unknown', IOSTAT=erro)
  if (erro/=0) stop "Error opening bands_bse output file"    

  OPEN(UNIT=700, FILE=trim(outputfolder)//"bse_coeff.dat",STATUS='unknown', IOSTAT=erro)
  if (erro/=0) stop "Error opening bands_bse output file"  


  call cpu_time(t0)
  call date_and_time(VALUES=values)
  call OMP_SET_NUM_THREADS(nthreads)

  call hamiltonian_input_read(200,params)
! if (berryexc) then
!   if (.not. allocated(Mmn)) allocate(Mmn(nbands,nbands,nkpts_mmn,nntot))
!   call mmn_input_read_elements(800,params_mmn,w90basis,nkpts,nntot)
! end if
  !ediel(2) = edielh

  read(500,*) nks
  read(500,*) nkpts
  print*, 'nks=', nks
  allocate(ks(nks,3))

  do i=1,nks
    read(500,*) ks(i,1),ks(i,2),ks(i,3)
  end do

  !termino leitura parametros
  !parametros do calculo
  !termino parametros calculo 

  ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
  print*, 'ngkpt', ngkpt
  dimbse = ngkpt*nc*nv

  ! unknown for why we use half
  ! trying for full grid, since dimbse is defined in this way
! nqpts = half*nks*ngkpt
  nqpts = ngkpt
  !call alat(systype,rlat,a)
  !lc = rlat3(3)

  !if (systype .eq. "2D") then
  !  ed = (ediel(1)+ediel(3))/2.
  !  r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
  !else
  !  r0 = 1.0
  !  ed = ediel(2)
  !end if

  !definindo kpath
  !allocate(qauxv(nkpts*(nks-1),4))
  allocate(qauxv((nks/2)*nkpts,4))  

  call kpathbse(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,qauxv)

  !Informações para o arquivo de log do calculo

  write(300,*) 'threads:', nthreads
  write(300,*)
  write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
  write(300,*)
  write(300,*) 'ktol-coulomb:', ktol
  write(300,*)
  write(300,*) 'kmesh shift:',mshift(1),mshift(2),mshift(3)
  write(300,*)
  write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
  write(300,*)
  write(300,*) 'number of kpoints in the path:','   ',(nks/2)*nkpts
  write(300,*)
  write(300,*)
  write(300,*) 'begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
  write(300,*) 

  call flush(300)

  allocate(kpt(ngkpt,3))
  allocate(kpt_red(ngkpt,3))
  
  !shift= 0.0
  call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)

  allocate(eaux(w90basis),vaux(w90basis,w90basis))
  allocate(energy(ngkpt,nc+nv),vector(ngkpt,nc+nv,w90basis))

  allocate(nocpk(ngkpt))
  allocate(nocpq(ngkpt))

! $omp parallel do
  do i=1,ngkpt
    call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
          rlat,rvec,hopmatrices,ihopmatrices,efermi,eaux,vaux)
    do j=1,nc+nv
      energy(i,j)= eaux(nocpk(i)-nv+j)
    end do

    do l=1,nc+nv
      do h=1,w90basis
        vector(i,l,h)=vaux(nocpk(i)-nv+l,h)
      end do
    end do
  end do
! $omp end parallel do

  deallocate(eaux,vaux)


  !definindo os numeros quanticos dos estados
  !write(*,*) "quantum numbers for exciton basis set finished"

  write(300,*) "quantum numbers for exciton basis set finished"
  call flush(300)  

  counter=counter-1 !numero total de estados para equação bse

  allocate(qpt(ngkpt,3))
  allocate(exk(dimbse,nkpts*(nks-1)))
  allocate(wfk(dimbse,dimbse,nkpts*(nks-1)))

  print*,"Dimensions:", w90basis, nc,dimbse,(nks/2)*nkpts
  print*,"    w90basis: ", w90basis
  print*,"    nc:       ", nc
  print*,"    dimbse:   ", dimbse
  print*,"    k-points: ", (nks/2)*nkpts
  print*,"    q-points: ", nqpts
  allocate(A_table(excwff-excwf0+1,w90basis,w90basis,ngkpt,(nks/2)*nkpts))
  allocate(bse_table(dimbse,dimbse,nqpts),exc_exc_overlap(dimbse,nqpts,nqpts))
  print*, 'Shape of exc_exc_overlap: ', shape(exc_exc_overlap)
  print*, 'Shape of bse_table: ', shape(bse_table)
! if (berryexc) then
!   allocate(Mmn(nbands,nbands,nkpts_mmn,nntot))  
!   call mmn_input_read_elements(800,params_mmn,w90basis,nkpts,nntot)
! end if
  A_table = czero;  bse_table = czero; exc_exc_overlap = czero

  write(700,'(2a18,6a8)') 'Re(A^λ_cvkq)','Im(A^λ_cvkq)', 'i2 dim_bse', 'c', 'v', 'k', 'Q idx', 'E(eV)'
  Format = "(2F15.8,5I,F15.4)"

  do i=1,nks/2*nkpts
    !definindo os pontos q
    
    q(1)= qauxv(i,1)
    q(2)= qauxv(i,2) 
    q(3)= qauxv(i,3)
    q(4)= qauxv(i,4)

    !gerando o grid k+q , here I effectively already have the k+q point but I am not sure I know the inedex of this q
    
    call monhkhorst_packq(q,ngrid(1),ngrid(2),ngrid(3),mshift,&
              rlat(1,:),rlat(2,:),rlat(3,:),qpt)


    allocate(eaux(w90basis),vaux(w90basis,w90basis))

    allocate(energyq(ngkpt,nc+nv),vectorq(ngkpt,nc+nv,w90basis))

    allocate (stt(ngkpt*nc*nv,4))

  ! $omp parallel do
    do i2=1,ngkpt
      !write(*,*) i2
      call eigsys(nthreads,scs,exc,nocpq(i2),ffactor,qpt(i2,1),qpt(i2,2),qpt(i2,3),w90basis,nvec,&
                   rlat,rvec,hopmatrices,ihopmatrices,efermi,eaux,vaux)
      do j=1,nc+nv
        energyq(i2,j)= eaux(nocpq(i2)-nv+j) 
      end do

      do l=1,nc+nv
        do h=1,w90basis
          vectorq(i2,l,h)=vaux(nocpq(i2)-nv+l,h)
        end do
      end do
    end do
  ! $omp end parallel do
    call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpq,stt)

    deallocate(eaux,vaux)

    allocate(hbse(dimbse,dimbse),W(dimbse))

    hbse=0.

  !$omp parallel do 
        !collapse(2)

    do i2=1,dimbse
      do j=i2,dimbse
        hbse(i2,j)= matrizelbsekq(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,q,rlat,stt(i2,:),energyq(stt(i2,4),stt(i2,3))&
          ,energy(stt(i2,4),stt(i2,2)),vectorq(stt(i2,4)&
          ,stt(i2,3),:) ,vector(stt(i2,4),stt(i2,2),:),kpt(stt(i2,4),:),stt(j,:)&
          ,energyq(stt(j,4),stt(j,3)),energy(stt(j,4),stt(j,2))&
          ,vectorq(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:))
      end do
    end do

  !$omp end parallel do

        !call LA_HEEVR( hbse, W, JOBZ='N', UPLO='U', ABSTOL=ABSTOL, INFO=INFO ) 
        !LWORK = 2*dimbse+dimbse**2
        !LIWORK = 3 + 5*dimbse
        !LRWORK = 1 + 5*dimbse + 2*dimbse**2

  !allocate(WORK(LWORK))
  !allocate (IWORK(LIWORK))
  !allocate (RWORK(LRWORK))
        !CALL ZHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
        !LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        !LRWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        !LIWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

        !CALL ZHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,&
  !        RWORK, LRWORK, IWORK, LIWORK,INFO )
        !IF( INFO.GT. 0 ) THEN
        !WRITE(*,*)'The algorithm failed to compute eigenvalues.'
        ! STOP
        !END IF   

    allocate (RWORK(3*dimbse-2))
    LWMAX = 2*dimbse-1
    allocate(WORK(LWMAX))
  
    if (bsewf) then
      LWORK = -1
      CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      IF( INFO.GT. 0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF 
    else
      LWORK = -1
      CALL ZHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      CALL ZHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
        STOP
      END IF 
    end if  

    !write(*,*) W(1)

    do i2=1,dimbse
      exk(i2,i)=W(i2)
      wfk(i2,:,i)=hbse(i2,:) ! should it be hbse(i2,:)?
    !  write(400,*) qauxv(i,1),W(i2)
    !  call flush(400)
    end do

    bse_table(1:dimbse,1:dimbse,i) = hbse(1:dimbse,1:dimbse)

    !initialize bse coefficiente
    ! bse_coefficient%nc = nc
    ! bse_coefficient%nv = nv
    ! bse_coefficient%nkpts = ngkpt
    ! bse_coefficient%nT = excwff-excwf0
    ! bse_coefficient%lmbd = dimbse

    if (bsewf) then
      do i2=excwf0,excwff
        !write(*,*) "starting cycle i2", i2
        call excwfi2(outputfolder,ngkpt,kpt,q,nc,nv,nocpk,stt,W(i2),i2,i,hbse(:,i2)) !i is qptnum stt(ip,4) gives the k-point index
        if(berryexc) then
!         print*, '1st berryexc true, '
          do ip=1,ngkpt*nc*nv
            !write(*,*) "starting cycle", ip
            !write(*,*) "shapes", shape(stt), 'ta', shape(hbse), 'stt(ip)', stt(ip,4),i
            A_table(i2,nocpk(stt(ip,4))+stt(ip,3)-nv,nocpk(stt(ip,4))-nv+stt(ip,2),stt(ip,4),i) = hbse(ip,i2)

            ! write(*,*) 'bse indices bef', i2,nocpk(stt(ip,4))+stt(i,3)-nv,nocpk(stt(ip,4))-nv+stt(ip,2),stt(ip,4),i
            ! write(*,*) 'bse_cofficient', A_table(i2,nocpk(stt(ip,4))+stt(i,3)-nv,nocpk(stt(ip,4))-nv+stt(ip,2),stt(ip,4),i)
            ! write(*,*) 'bse indices aft', i2,nocpk(stt(ip,4))+stt(i,3)-nv,nocpk(stt(ip,4))-nv+stt(ip,2),stt(ip,4),i

            write(700,*) real(hbse(ip,i2)),aimag(hbse(ip,i2)),i2,nocpk(stt(ip,4))+stt(ip,3)-nv,&
&            nocpk(stt(ip,4))-nv+stt(ip,2),stt(ip,4),i, W(ip)
          end do
        endif
      end do
    else
     continue  
    end if

    deallocate(hbse,W)
    deallocate(energyq,vectorq)
    deallocate(stt)
    deallocate(RWORK,WORK)
    !deallocate (IWORK)

    qpt=0.0
    nocpq=0

    write(300,*) 'progress:',i,'/',(nks/2)*nkpts
    call flush(300)

  end do !end do i =qpoint

  if (berryexc) then
    print*, 'allocating  neighbours'

    allocate(q_neigh(ngrid(1)*ngrid(2)*ngrid(3),6))
    allocate(chern_exc_exc(dimbse,3))
    allocate(chern_exc_hole(dimbse,3))
    allocate(chern_exc_electron(dimbse,3))
    allocate(kmq_grid(ngrid(1)*ngrid(2)*ngrid(3),ngrid(1)*ngrid(2)*ngrid(3),3)) ! ik-iq point =  kmq_grid(ik,iq,1:3)
    allocate(kmq_neighbours(ngrid(1)*ngrid(2)*ngrid(3),ngrid(1)*ngrid(2)*ngrid(3),6)) ! neighbours in the k-q and k+q grids
    
    chern_exc_exc = czero; chern_exc_hole = czero; chern_exc_electron = czero
    blat = 0.
    
    do n = 1, nqpts
      kpt_red(n,1:3) = matmul(rlat(1:3,1:3),kpt(n,1:3))/(2*pi)
    enddo
    do i=1, ngkpt
      print*, kpt_red(i,:) 
    enddo      
    print*, 'mat', rlat(1:3,1:3)
    print*, 'finding neighbours for each k-point'
    ! get the list of neighbours in reduced coordinate space
    call find_neighbour(ngrid,kpt_red,q_neigh)

    print*, 'calling neighbour indices'
    ! get the ijk indices of each kpoint
    call get_ijk(kpt_red,ijk,ngrid,rlat)

    ! get the reciprocal lattice vector matrix
    call recvec(rlat(1,:),rlat(2,:),rlat(3,:),blat(1,:),blat(2,:),blat(3,:))
 
    ! get the reduced lattice differential dk_red
    dk_red(1:3) = 1.0/ngrid(1:3)
print*, 'delta-k in reduced coordinates: ', dk_red 
    ! get the cartesian lattice differential dk_car
    dk_car = matmul(blat,dk_red)
print*, 'delta-k in cartesian coordinates: ', dk_red
!   A_table(:,:,:,:,:) = 1.0
    ! calculate all overlaps between exciton states at different q-points
    call calculate_exc_exc_overlap(bse_table,exc_exc_overlap,dimbse,nqpts)
    print*, 'passed calculate_exc_exc_overlap'
    
    do ic = 1,3
      print*, 'bse ic = ', ic
      do n  = 1, dimbse 
        do iq = 1, ngrid(1)*ngrid(2)*ngrid(3)
          if (kpt(iq,ic) == 0.0) then
            ! what are the i,j,k indices of this point?
            i = ijk(iq,1); j = ijk(iq,2); k = ijk(iq,3)
            ! tmp_o is not needed anymore but I keep it because the name is shorter
            tmp_o(1) = exc_exc_overlap(n,iq,q_neigh(iq,1)) ! overlap <(i,j,k)|(i+1,j,k)>
            tmp_o(2) = exc_exc_overlap(n,iq,q_neigh(iq,2)) ! overlap <(i,j,k)|(i,j+1,k)>
            tmp_o(3) = exc_exc_overlap(n,iq,q_neigh(iq,3)) ! overlap <(i,j,k)|(i,j,k+1)>
            
            rot_overlap = czero
            do id1 = 1, 3
              do id2 = 1, 3 
                ! rotation from reduced coordinates to cartesian coordinates
                rot_overlap = rot_overlap + rlat(id2,id1)*dk_car(id1)*tmp_o(id2)/dk_red(id2)               
              enddo
            enddo
            ! boundary check
            if (iq_d1 .ge. ngrid(dir1) .or. iq_d2 .ge. ngrid(2) .or. iq_d3 .ge. ngrid(3)) then
              chern_exc_exc(n,ic) = chern_exc_exc(n,ic) + rot_overlap
            ! if inside the BZ add 2*Im[\nabla A^* x \nabla A]
            else
              chern_exc_exc(n,ic) = chern_exc_exc(n,ic) + 2*aimag(rot_overlap)
            endif
          endif 
        enddo !iq
      enddo !n
    enddo !ic 

chern_exc_exc = 1.0/ngkpt * chern_exc_exc ! volume factor
print*, 'exc-exc cycle ended'
    print*, 'Reading file: ', params_mmn
    call mmn_input_read_dimensions(800,params_mmn)
    if (.not. allocated(o_mmn)) allocate(o_mmn(w90basis,w90basis,nqpts,nqpts))
    call mmn_input_read_elements(800,params_mmn,w90basis,nqpts,nqpts,o_mmn)
print*, 'starting exc-hole cycle'
print*, 'finding indices of k-q grid'

    ! reinitialize arrays
    kmq_grid = 0.0
    kmq_ijk = 0.0
    kmq_index = 0.0

    ! build 1st the k-q grid, which is the same for all exction states, and find the neighours of all k-q points
    do ik = 1, nqpts ! hole k-point
      do iq = 1, nqpts ! exciton q-point
        ! compute k-q in reduced coordinates
        kmq_grid(ik,iq,1:3) = kpt_red(ik,1:3) - kpt_red(iq,1:3)
        ! if the coordinate is larger/smaller than 0 or larger than 1 we need to bring it back into the {0,...,1 - 1/NGX/Y/Z} set
        ! this is equivalent to applying a translation ablong b1, b2, and or b3
        do ic = 1, 3
          if (kmq_grid(ik,iq,ic) .lt. 0.0-tol) then
            kmq_grid(ik,iq,ic) = kmq_grid(ik,iq,ic) + 1.
          elseif (kmq_grid(ik,iq,ic) .ge. 1.0-tol) then 
            kmq_grid(ik,iq,ic) = kmq_grid(ik,iq,ic) - 1.
          endif          
        enddo !ic
        kmq_index(ik,iq) = find_kpoint_index(kpt_red, nqpts, kmq_grid(ik,iq,1:3))
      enddo !iq
! PM: there's no need to call the routine again. Since find_kpoint_index gives the index, we already have the neighbours for each
! index from the exc-exc section
!    call find_neighbour(ngrid,kmq_grid(ik,:,:),kmq_neighbours(ik,:,:))
    enddo !ik

    do id = 1,3
      do n  = 1, dimbse
        do iq = 1, ngrid(1)*ngrid(2)*ngrid(3)
          if (kpt(iq,id) == 0.0) then
             do ic = 1, nc
              do iv1 = 1, nv
                do iv2 = 1, nv
                  do ik = 1, ngrid(1)*ngrid(2)*ngrid(3)
                    i_kmq = kmq_index(ik,iq)
                    do dir1 = 1, 3
                      do dir2 = 1,3
                        chern_exc_hole(n,id) = chern_exc_hole(n,id) + &
& A_table(n,ic,iv1,ik,iq)*conjg(A_table(n,ic,iv2,ik,iq))*o_mmn(iv1,iv2,i_kmq,q_neigh(i_kmq,dir2))*rlat(dir2,dir1)*dk_car(dir1)/dk_red(dir2)- &
& A_table(n,ic,iv1,ik,iq)*conjg(A_table(n,ic,iv2,ik,iq))*o_mmn(iv1,iv2,i_kmq,i_kmq)*              rlat(dir2,dir1)*dk_car(dir1)/dk_red(dir2)
                      enddo
                    enddo
                  enddo !ik 
                enddo !nv
              enddo !nv
            enddo !nc
          endif
        enddo !iq
      enddo !n
    enddo !ic    
print*, 'overlaps computed'

    chern_exc_hole = 1.0/ngkpt * chern_exc_hole
print*, 'exc-hole cycle ended', chern_exc_hole

print*, 'starting exc-electron cycle'
print*, 'building k+q grid and finding its neighbours'
    kmq_grid = 0.0
    kmq_ijk = 0.0
    kmq_index = 0.0

    ! build 1st the k-q grid, which is the same for all exction states, and find the neighours of all k+q points
    do ik = 1, nqpts ! hole k-point
      do iq = 1, nqpts ! exciton q-point
        do ic = 1,3
         ! compute k+q
         ! we can use kmq to store k+q after k-q is no longer needed
         kmq_grid(ik,iq,ic) = kpt_red(ik,ic) + kpt_red(iq,ic)
         ! if the coordinate is larger/smaller than 0 or larger than 1 we need to bring it back into the {0,...,1 - 1/NGX/Y/Z} set
         ! this is equivalent to applying a translation along b1, b2, and or b3
         if (kmq_grid(ik,iq,ic) .lt. 0.0-tol) then 
           kmq_grid(ik,iq,ic) = kmq_grid(ik,iq,ic) + 1.
         elseif (kmq_grid(ik,iq,ic) .ge. 1.0-tol) then
           kmq_grid(ik,iq,ic) = kmq_grid(ik,iq,ic) - 1.
         endif          
       enddo !ic
       kmq_index(ik,iq) = find_kpoint_index(kpt_red, nqpts, kmq_grid(ik,iq,1:3))
     enddo !iq
   enddo !ik
print*, 'found all k+q grid neighbours'

print*, 'computing overlaps'
    do id = 1,3
      do n  = 1, dimbse
        do iq = 1, ngrid(1)*ngrid(2)*ngrid(3)
          if (kpt(iq,id) == 0.0) then
            do ic1 = 1, nc
              do ic2 = 1, nc
                do iv = 1, nv
                  do ik = 1, ngrid(1)*ngrid(2)*ngrid(3)
                    i_kmq = kmq_index(ik,iq)
                    do dir1 = 1, 3
                      do dir2 = 1, 3
                        chern_exc_electron(n,id) = chern_exc_electron(n,id) + &
& A_table(n,ic1,iv,ik,iq)*conjg(A_table(n,ic2,iv,ik,iq))*o_mmn(ic1,ic2,i_kmq,q_neigh(i_kmq,dir2))*rlat(dir2,dir1)*dk_car(dir1)/dk_red(dir2) - &
& A_table(n,ic1,iv,ik,iq)*conjg(A_table(n,ic2,iv,ik,iq))*o_mmn(ic1,ic2,i_kmq,i_kmq)*              rlat(dir2,dir1)*dk_car(dir1)/dk_red(dir2)
                      enddo
                    enddo
                  enddo
                enddo !iv
              enddo !ic2
            enddo !ic1
          endif
        enddo !iq
      enddo !n
    enddo !ic

    chern_exc_electron = 1.0/ngkpt * chern_exc_electron
print*, 'exc-electron cycle ended'
print*, 'passed berryexc'
  endif !berryexc

  do i2=1,dimbse
    do i=1,(nks/2)*nkpts
      write(400,*) real(qauxv(i,1)),exk(i2,i)
!      do it = 1,dimbse
      write(600,*) real(qauxv(i,1)), 'it', it, sum(wfk(i2,:,i))
      !end do
      call flush(400)
      call flush(600)
    end do
    write(400,*)
    write(600,*)
  end do
  deallocate(energy,vector)
  deallocate(kpt)
  deallocate(kpt_red)
  deallocate(qpt)
  deallocate(exk,wfk)
  deallocate(nocpq,nocpk)
  deallocate(rvec,hopmatrices)
  deallocate(ihopmatrices,ffactor)
  deallocate(A_table, bse_table)
  if (allocated(qauxv)) deallocate(qauxv)
  if (allocated(ks)) deallocate(ks)
  if (allocated(o_mmn)) deallocate(o_mmn)
  if (allocated(chern_exc_exc)) deallocate(chern_exc_exc)
  if (allocated(exc_exc_overlap)) deallocate(exc_exc_overlap)
  call cpu_time(tf)
  call date_and_time(VALUES=values2)
  ! if (allocated(kmq_grid)) deallocate(kmq_grid)
  ! if (allocated(kmq_neighbours)) deallocate(kmq_neighbours)
  if (allocated(q_neigh)) deallocate(q_neigh)
  if (allocated(chern_exc_hole)) deallocate(chern_exc_hole)
  if (allocated(chern_exc_electron)) deallocate(chern_exc_electron)

  write(300,*)
  write(300,*) 'end','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
  write(300,*)
  close(200)
  close(300)
  close(400)
  close(500)
  close(700)
  close(600)
  close(800)
end subroutine bsebnds

! PM: there is no reason to have this as a subroutine since it only computes the product of 4 variables 
! and then returns the imaginary part of its logarithm
subroutine exc_overlap(Mmn,A_table,nbands,nkpoints,nntot,n,ik,iq,iqp,c1,c2,v1,v2, excwfftot)

  double complex  :: phase
  double complex  :: overlap
  double complex  :: tmp_A1,tmp_A2

  tmp_A1 = A_table(n,c1,v1,ik,iq)
  tmp_A2 = A_table(n,c2,v2,ik,iqp)
  overlap = conjg(tmp_A1)*tmp_A2*Mmn(c1,c2,ik+iq,ik+iqp)*Mmn(v2,v1,ik,ik)
  phase = aimag(log(overlap))
end subroutine exc_overlap

! computes the exciton-exciton term of the Chern number
! following \sum_{cv\kk}A_{cv\kk}^{\gl,*}(\qq_i)\[A_{cv\kk}^{\gl}(\qq_i + \gD\qq_i)-A_{cv\kk}^\gl(\qq_i)\]
! the best strategy is to agregate all {cv\kk} triplets and create an array that depends only on \gl, \qq_i, 
! right now it is easier to just compute the overlap between all q-points
subroutine calculate_exc_exc_overlap(bse_mat,overlap,matdim,nq)
  double complex, intent(in) :: bse_mat(1:matdim,1:matdim,1:nq), overlap(1:matdim,1:nq,1:nq)
  integer, intent(in) :: matdim, nq
  integer :: il
  print*, 'shape of arrays inside calculate_exc_exc_overlap:'
  print*, 'shape of bse_mat: ', shape(bse_mat)
  print*, 'shape of overlap: ', shape(overlap)
  print*, 'matdim:           ', matdim
  print*, 'nq:               ', nq
  do il = 1, matdim
    call zgemm('n','n', nq, nq, matdim, cone, conjg(bse_mat(il,1:matdim,1:nq)), matdim, bse_mat(il,1:matdim,1:nq), matdim, czero,overlap(il,1:nq,1:nq), nq)
  enddo
  print*, 'passed zgemm call for ', nq
end subroutine calculate_exc_exc_overlap

! computes the weighted elements of the M_nmkpb matrix created by pw2wannier90.x and the BSE coefficients
subroutine calculate_hole_hole_overlap(a_mat,m_mat,matdim,nc,nv,nk,nq,overlap,ntot)

  double complex :: a_mat(matdim,nc,nv,nk,nq), m_mat(nc+nv,nc+nv,nk,nk), overlap(matdim,nq)
  integer, intent(in) :: matdim, nc, nv, nk, nq, ntot
  integer :: il, iq, ic, iv1, iv2, ik1, ik2
  double complex :: tmp

!NEEDS TO BE REWORKED 
! WE CAN IDENTIFY WHICH PAIRS OF THE M_nmkpb matrix are actually points within the 1st BZ, but I have to rearrange the loop

  do il = 1, matdim
    do iq = 1, nq
      do iv1 = 1, nv
        do iv2 = 1, nv
          do ik1 = 1, nk
            do ik2 = 1, ntot
              do ic = 1, nc
! this is wrong but I cannot see now how this is handled in the code
! there has to be a way to identify which index is k - q - \Delta q_i
! and I do not think it is as it is done in exc_overlap, since there is no
! guarantee that you do not fall outside the size of the array
                overlap(il,iq) = overlap(il,iq) + conjg(a_mat(il,ic,iv1,ik,iq))*a_mat(il,ic,iv2,ik,iq)*m_mat(iv1,iv2,ik1-iq,ik2-iq)
              enddo ! nc
            enddo ! nk
          enddo ! nk
        enddo ! nv
      enddo ! nv
    enddo ! nq
  enddo  ! matdim

end subroutine calculate_hole_hole_overlap

