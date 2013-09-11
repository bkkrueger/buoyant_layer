subroutine init(ndump,nhist,nspec,mype,npes)
  use hydro_parameters
  use variables
  implicit none
  integer :: ndump,nhist,nspec
  integer :: mype,npes

! read input parameters (contains the grid size,...)
  if (mype==0) write (*,*) 'Reading input data...'
  call read_input
#if WITHMPI==1
  if (mype > 0) verbose=.false.
#endif

!calculate state vector sizes
  iu1=-2 ;  iu2=nx+3    ! iu1 = 1 - 3, iu2 = nx + 3  -->  3 guard cells
  ju1= 1 ;  ju2=1
  ku1= 1 ;  ku2=1
#if NDIM > 1
  ju1 = -2 ; ju2=ny+3  
#endif
#if NDIM==3
  ku1 = -2 ; ku2=nz+3  
#endif
! calculate flux vector sizes
  if1=1; if2=nx+1
  jf1=1; jf2=1
  kf1=1; kf2=1
#if NDIM > 1
  jf1=1; jf2=ny+1
#endif
#if NDIM==3
  kf1=1; kf2=nz+1
#endif

  if (verbose) write (*,*) 'Allocating arrays...'
  allocate(uin         (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3       )) ; uin=0.d0
  allocate(old_uin     (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3       )) ; old_uin=0.d0
  allocate(gravin      (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim         )) ; gravin=0.d0
  allocate(flux        (1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar  ,1:ndim)) ; flux=0.d0
  allocate(flux_pre_tot(1:nvector,if1:if2,jf1:jf2,kf1:kf2         ,1:ndim)) ; flux_pre_tot=0.d0
  allocate(emfx        (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2                )) ; emfx=0.d0
  allocate(emfy        (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2                )) ; emfy=0.d0
  allocate(emfz        (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2                )) ; emfz=0.d0

  allocate(x(iu1:iu2))
  allocate(y(ju1:ju2))
  allocate(z(ku1:ku2))

  time=0.d0
  dt=0.d0

  call init_grid(mype,npes)

  if (verbose) write (*,*) 'Initialising problem...'
  call condinit(mype)
  call bval(uin,ngrid,time,dt,dx,dy,dz,bval_type,mype)

  ndump = 0
  nhist = 1
  nspec = 1

  if (restart>0) then
      ndump=restart
      call restart_run(uin,x,y,z,time,dt,dx,dy,dz,ndump,nhist,nspec,mype)
  endif
  
  call init_geom

  !call user_init
  call user_init(mype, npes)

  return
end subroutine init
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine read_input
  use hydro_parameters
  use variables
#if WITHMPI==1
  use mpi_var
#endif
  implicit none

  integer :: nx_glob,ny_glob,nz_glob
  real(dp), dimension(ndim) :: bval_in,bval_out

  namelist /start_params/restart,tlim,verbose
  namelist /scheme_params/nx_glob,ny_glob,nz_glob,riemann,riemann2d,slope_type,&
                         & bval_in,bval_out,courant,fargo
  namelist /model_params/geom,xmin,xmax,ymin,ymax,zmin,zmax,gamma,Omega0,ciso,nu,eta
  namelist /output_params/dthist,dtdump,dtspec,use_collective_io,store_ghost,&
                         & store_sides,inline_io

  call default_params(nx_glob,ny_glob,nz_glob,bval_in,bval_out)

  open(unit=1,file='input',status='old')
  read(1,start_params)
  read(1,scheme_params)
  read(1,model_params)
  read(1,output_params)

#if WITHMPI==1
  call mpi_param
  nx=nx_glob/nxslice ; ny=ny_glob/nyslice ; nz=nz_glob/nzslice
#else
  nx=nx_glob ; ny=ny_glob ; nz=nz_glob
#endif

  close(1)

  bval_type(1:ndim,1)=bval_in (1:ndim)
  bval_type(1:ndim,2)=bval_out(1:ndim)

  cartesian=.false. ; cylindrical=.false. ; spherical=.false.
  select case (geom)
  case (1)
     cartesian=.true.
  case (2)
     cylindrical=.true.
  case (3)
     spherical=.true.
  end select

  call init_solver

  return
end subroutine read_input
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine default_params(nx_glob,ny_glob,nz_glob,bval_in,bval_out)
  use hydro_parameters
  use variables
  implicit none

  integer :: nx_glob,ny_glob,nz_glob
  real(dp), dimension(ndim) :: bval_in,bval_out

  !'start_params' namelist parameters
  restart=0
  tlim   =1.d0
  verbose=.false.

  !scheme_params namelist parameters
  nx_glob=1 ; ny_glob=1 ; nz_glob=1
  riemann  ='hlld'
  riemann2d='hlld'
  slope_type=1
  bval_in(1:ndim)=1
  bval_out(1:ndim)=1
  courant=0.8
  fargo=.false.

  !models_params namelist parameters
  geom=1
  xmin=0.d0 ; xmax=1.d0
  ymin=0.d0 ; ymax=1.d0
  zmin=0.d0 ; zmax=1.d0
  gamma=1.001d0
  Omega0=0.d0
  ciso  =0.d0
  nu    =0.d0
  eta   =0.d0

  !output_params namelist parameters (negative value imply no output)
  dthist=-1.d0
  dtdump=-1.d0
  dtspec=-1.d0
  use_collective_io = .false.
  inline_io   = .false.
  store_ghost = .false.
  store_sides = 'none'

  return
end subroutine default_params
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_grid(mype,npes)
#if WITHMPI==1
  use mpi_var
#endif
  use variables
  use hydro_parameters
  implicit none

  integer :: mype,npes

  integer :: i,j,k

  if (verbose) write (*,*) 'Initialising grid...'

#if WITHMPI==1
  call grid_struct(mype,npes)
#endif

#if WITHMPI==1
  dx=(xmax-xmin)/(nx*nxslice)
  dy=(ymax-ymin)/(ny*nyslice)
  dz=(zmax-zmin)/(nz*nzslice)
#else
  dx=(xmax-xmin)/nx
  dy=(ymax-ymin)/ny
  dz=(zmax-zmin)/nz
#endif

  if (mype == 0) then
     write(*,*) ""
     write(*,*) 'Mesh size:'
  end if
  if (mype==0) write(*,*) '  dx',dx
#if WITHMPI==1
  x=xmin + dx*nx*xposition + dx/2 + (/ (i-1,i=iu1,iu2) /)*dx
#else
  x=xmin +                   dx/2 + (/ (i-1,i=iu1,iu2) /)*dx
#endif
#if NDIM > 1
  if (mype==0) write(*,*) '  dy',dy
#if WITHMPI==1
  y=ymin + dy*ny*yposition + dy/2 + (/ (j-1,j=ju1,ju2) /)*dy
#else
  y=ymin +                   dy/2 + (/ (j-1,j=ju1,ju2) /)*dy
#endif
#endif
#if NDIM==3
  if (mype==0) write(*,*) '  dz',dz
#if WITHMPI==1
  z=zmin + dz*nz*zposition + dz/2 + (/ (k-1,k=ku1,ku2) /)*dz
#else
  z=zmin +                   dz/2 + (/ (k-1,k=ku1,ku2) /)*dz
#endif
#endif
   if (mype==0) write(*,*) ""

  return
end subroutine init_grid
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_geom
  use hydro_parameters
  use const
  use variables
  implicit none

  integer :: ilo,ihi,jlo,jhi,klo,khi
  integer :: i,j,k,idim

  allocate(dv(iu1:iu2,ju1:ju2,ku1:ku2))        ; dv=0.d0
  allocate(ds(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)) ; ds=0.d0

  ilo=min(1,iu1+1); ihi=max(1,iu2-1)
  jlo=min(1,ju1+1); jhi=max(1,ju2-1)
  klo=min(1,ku1+1); khi=max(1,ku2-1)

  if (cartesian) then
     do k=klo,khi
        do j=jlo,jhi
           do i=ilo,ihi
              dv(i,j,k)=dx*dy*dz
              do idim=1,ndim
                 if (idim==1) ds(i,j,k,idim)=dy*dz
                 if (idim==2) ds(i,j,k,idim)=dx*dz
                 if (idim==3) ds(i,j,k,idim)=dx*dy
              enddo
           end do
        end do
     end do
  endif
  if (cylindrical) then
     do k=klo,khi
        do j=jlo,jhi
           do i=ilo,ihi
              dv(i,j,k)=dx*x(i)*dy*dz
              do idim=1,ndim
                 if (idim==1) ds(i,j,k,idim)=half*(x(i)+x(i-1))*dy*dz
                 if (idim==2) ds(i,j,k,idim)=dx*dz
                 if (idim==3) ds(i,j,k,idim)=x(i)*dx*dy
              enddo
           end do
        end do
     end do
  endif
  if (spherical) then
     do k=klo,khi
        do j=jlo,jhi
           do i=ilo,ihi
              dv(i,j,k)=dx*x(i)*dy*x(i)*sin(y(j))*dz
              do idim=1,ndim
                 if (idim==1) ds(i,j,k,idim)=half*(x(i)+x(i-1))*dy*dz
                 if (idim==2) ds(i,j,k,idim)=half*(x(i)+x(i-1))*dx*dz
                 if (idim==3) ds(i,j,k,idim)=x(i)*sin(z(k))*dx*dy
              enddo
           end do
        end do
     end do
  endif

end subroutine init_geom
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_solver
  use hydro_parameters
  implicit none

  select case (riemann)
  case ('roe')
     iriemann=iroe
  case ('llf')
     iriemann=illf
  case ('hll')     
     iriemann=ihll
  case ('hlld')
     iriemann=ihlld
  case ('upwind')
     iriemann=iupwind
  case ('hydro')
     iriemann=iacoustic
  end select

  select case (riemann2d)
  case ('roe')
     iriemann2d=iroe
  case ('llf')
     iriemann2d=illf
  case ('hll')     
     iriemann2d=ihll
  case ('hlld')
     iriemann2d=ihlld
  case ('upwind')
     iriemann2d=iupwind
  case ('hydro')
     iriemann2d=iacoustic
  case ('hllf')
     iriemann2d=ihllf
  case ('hlla')
     iriemann2d=ihlla
  end select

  return
end subroutine init_solver
