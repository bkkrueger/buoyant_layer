#if HDF5==1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine restart_run(uin,x,y,z,time,dt,dx,dy,dz,ndump,nhist,nspec,mype)
  use hdf5
#if WITHMPI==1
  use mpi_var
#endif
  use hydro_parameters
  use variables, only : use_collective_io
  use rdwrt_h5
  use steady_state, only : ss0, SS_DENS, SS_MOMX, SS_ENER, SS_NVAR
  use perturb, only : x_advc_old, x_geom_old, &
                      dv_grav_old, dv_buoy_old, dv_drag_old
  implicit none
  real(dp)::time,dt,dx,dy,dz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(iu1:iu2)::x
  real(dp),dimension(ju1:ju2)::y
  real(dp),dimension(ku1:ku2)::z
  integer::ndump,nhist,nspec,mype

  character(LEN=80) :: filename

  ! HDF5 vars
  integer        :: ierr
  integer(hid_t) :: file_id       ! File ID

  integer, parameter :: nb_real=10,nb_int=9,nb_mpi=9

  integer , dimension(nb_int ) :: para_int
  integer , dimension(nb_mpi ) :: para_mpi
  real(dp), dimension(nb_real) :: para_real

  integer :: n1,n2,n3

  if (verbose) write (*,*) 'Entering newstart...'

#if PARAHDF5==1
     call get_dumpname(ndump,0,filename)
#else
     call get_dumpname(ndump,mype,filename)
#endif

  call H5open_f(ierr)
  call H5Fopen_f(trim(filename),H5F_ACC_RDONLY_F,file_id,ierr)

  !Dataset "para_real"; Type real; dim=1; size=nb_real
  call get_1d_array_h5(file_id,"para_real",para_real,nb_real)
  time=para_real(1)
  dt  =para_real(2)
  dx  =para_real(3)
  dy  =para_real(4)
  dz  =para_real(5)
  x_advc_old = para_real(6)
  x_geom_old = para_real(7)
  dv_grav_old = para_real(8)
  dv_buoy_old = para_real(9)
  dv_drag_old = para_real(10)

  !Dataset "para_int"; Type integer; dim=1; size=nb_int
  call get_1d_array_int_h5(file_id,"para_int",para_int,nb_int)
  ndump  =para_int(1)
  nhist  =para_int(2)
  nspec  =para_int(3)
  nx     =para_int(4)
  ny     =para_int(5)
  nz     =para_int(6)
#if WITHMPI==1
  nxslice=para_int(7)
  nyslice=para_int(8)
  nzslice=para_int(9)
#endif
  
  ! Size of the datasets
  n1=iu2-iu1+1
  n2=ju2-ju1+1
  n3=ku2-ku1+1
  !Dataset "x", "y" and "z"; Type real; dim=1, size=n1, n2 and n3
  if(use_collective_io) then
#if PARAHDF5==1
     call get_2d_array_h5(file_id, "x", x, n1, 1, mype)
     call get_2d_array_h5(file_id, "y", y, n2, 1, mype)
     call get_2d_array_h5(file_id, "z", z, n3, 1, mype)
#endif
  else
     call get_1d_array_h5(file_id,"x",x,n1)
     call get_1d_array_h5(file_id,"y",y,n2)
     call get_1d_array_h5(file_id,"z",z,n3)
  endif

  !Stable state
  allocate(ss0(iu1:iu2,SS_NVAR))
  if(use_collective_io) then
#if PARAHDF5==1
     call get_1d_array_h5(file_id,"ss_dens",ss0(:,SS_DENS),n1,mype)
     call get_1d_array_h5(file_id,"ss_momx",ss0(:,SS_MOMX),n1,mype)
     call get_1d_array_h5(file_id,"ss_ener",ss0(:,SS_ENER),n1,mype)
#endif
  else
     call get_2d_array_h5(file_id,"ss",ss0,n1,SS_NVAR)
  endif

  !Dataset "uin"; Type real; dim=4; size=(n1,n2,n3,nvar+3)
  if(use_collective_io) then
#if PARAHDF5==1
     call get_3d_array_h5(file_id,"rho",     uin(1,:,:,:,1),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"E",       uin(1,:,:,:,5),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"rho_vx",  uin(1,:,:,:,2),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"rho_vy",  uin(1,:,:,:,3),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"rho_vz",  uin(1,:,:,:,4),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"Bx",      uin(1,:,:,:,6),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"By",      uin(1,:,:,:,7),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"Bz",      uin(1,:,:,:,8),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"Bxr",     uin(1,:,:,:,9),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"Byr",     uin(1,:,:,:,10),n1,n2,n3,mype)
     call get_3d_array_h5(file_id,"Bzr",     uin(1,:,:,:,11),n1,n2,n3,mype)
#endif
  else
     call get_4d_array_h5(file_id,"uin",uin,n1,n2,n3,nvar+3)
  endif

#if WITHMPI==1
  !Dataset "para_mpi"; Type integer; dim=1; size=nb_mpi
  if(use_collective_io) then
#if PARAHDF5==1
     call get_2d_array_int_h5(file_id,"para_mpi",para_mpi,nb_mpi &
          & ,1 , mype)
     xleft    =para_mpi(1)
     xright   =para_mpi(2)
     yleft    =para_mpi(3)
     yright   =para_mpi(4)
     zleft    =para_mpi(5)
     zright   =para_mpi(6)
     xposition=para_mpi(7)
     yposition=para_mpi(8)
     zposition=para_mpi(9)
#endif
  else
     call get_1d_array_int_h5(file_id,"para_mpi",para_mpi,nb_mpi)
     xleft    =para_mpi(1)
     xright   =para_mpi(2)
     yleft    =para_mpi(3)
     yright   =para_mpi(4)
     zleft    =para_mpi(5)
     zright   =para_mpi(6)
     xposition=para_mpi(7)
     yposition=para_mpi(8)
     zposition=para_mpi(9)
  endif
#endif
  
  ! Terminate access to HDF5
  call H5Fclose_f(file_id,ierr)
  call H5close_f(ierr)

  return
end subroutine restart_run
#elif PNETCDF==1
!###########################################################
!###########################################################
!###########################################################
subroutine restart_run(uin,x,y,z,time,dt,dx,dy,dz,ndump,nhist,nspec,mype)
   stop
  use hydro_parameters
  use pnetcdf
#if WITHMPI==1
  use mpi_var
#endif
  use variables, only : inline_io

  implicit none

#if WITHMPI==1
  include 'mpif.h'
#endif

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) ::uin
  real(dp),dimension(iu1:iu2) ::x
  real(dp),dimension(ju1:ju2) ::y
  real(dp),dimension(ku1:ku2) ::z
  real(dp) ::time, dt, dx, dy, dz
  integer  :: ndump, nhist, nspec, mype, npes
  integer, parameter :: nb_real=5, nb_int=9, nb_mpi=9
  integer(8) :: nbreal, nbint, nbmpi
  integer , dimension(nb_int)  :: para_int
  real(dp), dimension(nb_real) :: para_real
  real(dp), dimension(nb_mpi)  :: para_mpi
  real(dp), dimension(6) :: boxSize

  character(LEN=80) :: filename

  ! PnetCDF variables
  integer :: nout, ncid
  integer(kind=MPI_OFFSET_KIND), dimension(3) :: dims, start, count
  integer :: ndims, nvars, natts, nulms
  character(Len=20) :: name
  integer :: type, ndimvar, dimid(3), id, diminline, i

  real(dp) :: d0,beta_init,beta,B0_init,B0
  namelist /restart_params/d0,beta_init,beta

  if(verbose) write(*,*) 'Entering restart (with PnetCDF)!'

  call get_dumpname(ndump,0,filename)

  ! Open netCDF file
  nout = nfmpi_open(MPI_COMM_WORLD, filename, ior(NF_NOWRITE, NF_64BIT_OFFSET) &
       & ,MPI_INFO_NULL, ncid)

  ! Get basic informations (number of dimensions, variables...)
  nout = nfmpi_inq(ncid, ndims, nvars, natts, nulms)

  ! Loop over variables to retrieve metadata
  do i = 1, nvars
     nout = nfmpi_inq_var(ncid, i, name, type, ndimvar, dimid, natts)

     if(ndimvar.eq.1) then
        if(trim(name).eq.'para_int')then
           nbint = nb_int
           nout    = nfmpi_get_vara_int_all(ncid, i, (/1_8/), (/nbint/) &
                , para_int)
           ndump   = para_int(1)
           nhist   = para_int(2)
           nspec   = para_int(3)
           nx      = para_int(4)
           ny      = para_int(5)
           nz      = para_int(6)
#if WITHMPI==1
           nxslice = para_int(7)
           nyslice = para_int(8)
           nzslice = para_int(9)
#endif
        endif
        if(trim(name).eq.'para_real')then
           nbreal = nb_real
           nout   = nfmpi_get_vara_double_all(ncid, i, (/1_8/), (/nbreal/) &
                , para_real)
           time = para_real(1)
           dt   = para_real(2)
           dx   = para_real(3) 
           dy   = para_real(4) 
           dz   = para_real(5) 
        end if
        if(trim(name).eq.'boxSize')then
           nout = nfmpi_get_vara_double_all(ncid, i, (/1_8/), (/6_8/), boxSize)
           xmin = boxSize(1)
           xmax = boxSize(2)
           ymin = boxSize(3)
           ymax = boxSize(4)
           zmin = boxSize(5)
           zmax = boxSize(6)
        endif
#if WITHMPI==1
        if(trim(name).eq.'para_mpi')then
           nbmpi = nb_mpi
           nout = nfmpi_get_vara_double_all(ncid, i, (/1_8/), (/nbmpi/) &
                , para_real)
           xleft     = para_mpi(1)
           xright    = para_mpi(2)
           yleft     = para_mpi(3)
           yright    = para_mpi(4)
           zleft     = para_mpi(5)
           zright    = para_mpi(6)
           xposition = para_mpi(7)
           yposition = para_mpi(8)
           zposition = para_mpi(9)
        endif
#endif
     endif
  enddo

  ! Synchronize all threads to ensure metadata reading
  call MPI_Barrier(MPI_COMM_WORLD, nout)

  dims = (/ nx+2*nghost, ny+2*nghost, nz+2*nghost /)

  ! Loop over variables to retrieve data
  do i = 1, nvars
     nout = nfmpi_inq_var(ncid, i, name, type, ndimvar, dimid, natts)
     
     if(trim(name).eq.'rho') id = 1
     if(trim(name).eq.'E') id = 5
     if(trim(name).eq.'rho_vx') id = 2
     if(trim(name).eq.'rho_vy') id = 3
     if(trim(name).eq.'rho_vz') id = 4
     if(trim(name).eq.'Bx') id = 5
     if(trim(name).eq.'By') id = 6
     if(trim(name).eq.'Bz') id = 7
     
     if(inline_io) then
        diminline = mype*dims(3)+1
        start = (/ 1, 1, diminline /)
     else
        start = (/ xposition, yposition, zposition /)*dims + 1
     endif
     count = dims

     if(ndimvar.eq.3) then
        nout = nfmpi_get_vara_double_all(ncid, i, start, count, uin(1,:,:,:,id))
     endif
  enddo

  ! Close the file
  nout = nfmpi_close(ncid)

  return
end subroutine restart_run
#else
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine restart_run(uin,x,y,z,time,dt,dx,dy,dz,ndump,nhist,nspec,mype)
   stop
#if WITHMPI==1
  use mpi_var
#endif
  use hydro_parameters
  implicit none
  real(dp)::time,dt,dx,dy,dz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(iu1:iu2)::x
  real(dp),dimension(ju1:ju2)::y
  real(dp),dimension(ku1:ku2)::z
  integer::ndump,nhist,nspec,mype

  integer :: a,b,c
  character(LEN=6) :: snumfile,n_mype
  character(LEN=80) :: filename,filedir

  if (verbose) write (*,*) 'Entering newstart...'

  call get_dumpname(ndump,mype,filename)

  open(unit=2,file=trim(filename),status='unknown',form='unformatted')
  read (2) time,dt,dx,dy,dz
  read (2) ndump,nhist,nspec
  read (2) nx,ny,nz
#if WITHMPI==1
  read(2) nxslice,nyslice,nzslice
#else
  read(2) a,b,c
#endif
  read (2) x,y,z
  read (2) uin

#if WITHMPI==1
  read(2) xleft,xright,yleft,yright,zleft,zright
  read(2) xposition,yposition,zposition
#endif

  close(2)

  return
end subroutine restart_run
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_dumpname(ndump,mype,filename)
  implicit none

  integer::ndump,mype
  character(LEN=*) :: filename

  character(LEN=6) :: snumfile,n_mype
  character(LEN=80) :: filedir

  call convtoasc(ndump,snumfile)
  call convtoasc(mype,n_mype  )
  filedir='output_'//trim(snumfile)//'/'
  filename=trim(filedir)//'slices.'//trim(n_mype)

end subroutine get_dumpname
