!> \file output_collective.f90

!###########################################################
!###########################################################
!###########################################################
! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)
#if HDF5==1
subroutine output(uin,x,y,z,time,dt,dx,dy,dz,nhist,ndump,nspec,mype,npes)
  use hdf5
  use hydro_parameters
  use rdwrt_h5
#if WITHMPI==1
  use mpi_var
#endif
  use variables, only : use_collective_io, inline_io
  use steady_state, only : ss0, SS_DENS, SS_MOMX, SS_ENER, SS_NVAR
  use perturb, only : x_advc_old, x_geom_old, &
                      dv_grav_old, dv_buoy_old, dv_drag_old

  implicit none

#if WITHMPI==1
  include 'mpif.h'
#endif

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) ::uin
  real(dp),dimension(iu1:iu2) ::x
  real(dp),dimension(ju1:ju2) ::y
  real(dp),dimension(ku1:ku2) ::z
  real(dp) ::time,dt,dx,dy,dz
  integer :: ndump,nhist,nspec,mype,npes

  character(LEN=80) :: filename

  ! HDF5 vars
  integer        :: ierr
  integer(hid_t) :: file_id                   ! File ID
  integer(hid_t) :: fapl_id                   ! File property id (used in parallel version only)
  integer(hid_t) :: driver_id  ! low-level file driver identifier

  integer, parameter :: nb_real=10,nb_int=9,nb_mpi=9

  integer , dimension(nb_int ) :: para_int
  integer , dimension(nb_mpi ) :: para_mpi
  real(dp), dimension(nb_real) :: para_real
  real(dp), dimension(6) :: boxSize

  integer :: n1,n2,n3
  real(dp) :: tbegin, tend

  if (verbose) write (*,*) 'Entering output (hdf5)...'
  ! call cpu_time(tbegin)
#if WITHMPI==1
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tbegin = MPI_Wtime()
#endif

#if PARAHDF5==1
  if(use_collective_io) then
     ! get HDF5 filename (use mype=0 since there will be only 1 single file)
     call get_filename(ndump,0,'slices',filename)
  else
     call get_filename(ndump,mype,'slices',filename)
  endif
#else
  call get_filename(ndump,mype,'slices',filename)
#endif

  ! initialize hdf5 fortran interface
  call H5open_f(ierr)

  ! When using MPI, parameter use_collective_io allows to choose between
  ! collective IO, i.e. one file single or the standard way (one file per MPI process)
#if PARAHDF5==1
  if(use_collective_io) then
     ! create hdf5 property id for parallel file access
     call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
     call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
     call H5Pget_driver_f(fapl_id, driver_id, ierr)
     
     if( driver_id /= H5FD_MPIO_F) then
        write(*,*) 'Wrong driver information returned'
     endif
     
     call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, access_prp=fapl_id)
  else
     ! standard way (one file per MPI process)
     call H5Fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr)
  endif
#elif WITHMPI==1
  ! standard way (one file per MPI process)
  call H5Fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr)
#else
  ! sequential version
  call H5Fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr)
#endif

  !Dataset "para_real"; Type real; dim=1; size=nb_real
  !para_real=(/ time,dt,dx,dy,dz /)
  para_real=(/ time, dt, dx, dy, dz, x_advc_old, x_geom_old, &
               dv_grav_old, dv_buoy_old, dv_drag_old /)
  call dump_1d_array_h5(file_id,"para_real",para_real,nb_real)

  !Dataset "para_int"; Type integer; dim=1; size=nb_int
#if WITHMPI==1
  para_int=(/ ndump,nhist,nspec,nx,ny,nz,nxslice,nyslice,nzslice/)
#else
  para_int=(/ ndump,nhist,nspec,nx,ny,nz,1,1,1/)
#endif
  call dump_1d_array_int_h5(file_id,"para_int",para_int,nb_int)

  !Dataset "boxSize"; Type real; dim=1; size=6 - ADD from Sebastien's version
  boxSize=(/ xmin,xmax,ymin,ymax,zmin,zmax /)
  call dump_1d_array_h5(file_id,"boxSize",boxSize,6)

  ! Size of the datasets
  n1=iu2-iu1+1
  n2=ju2-ju1+1
  n3=ku2-ku1+1
  !Dataset "x", "y" and "z"; Type real; dim=1, size=n1, n2 and n3
  if(use_collective_io) then
#if PARAHDF5==1
     call dump_2d_array_h5(file_id, "x", x, n1, 1, mype)
     call dump_2d_array_h5(file_id, "y", y, n2, 1, mype)
     call dump_2d_array_h5(file_id, "z", z, n3, 1, mype)
#endif
  else
     call dump_1d_array_h5(file_id,"x",x,n1)
     call dump_1d_array_h5(file_id,"y",y,n2)
     call dump_1d_array_h5(file_id,"z",z,n3)
  endif

  !Stable state
  if(use_collective_io) then
#if PARAHDF5==1
     call dump_1d_array_h5(file_id,"ss_dens",ss0(:,SS_DENS),n1,mype)
     call dump_1d_array_h5(file_id,"ss_momx",ss0(:,SS_MOMX),n1,mype)
     call dump_1d_array_h5(file_id,"ss_ener",ss0(:,SS_ENER),n1,mype)
#endif
  else
     call dump_2d_array_h5(file_id,"ss",ss0,n1,SS_NVAR)
  endif

  !Dataset "uin"; Type real; dim=4; size=(n1,n2,n3,nvar)
  if(use_collective_io) then
#if PARAHDF5==1
     call dump_3d_array_h5(file_id,"rho",    uin(1,:,:,:,1),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"E",      uin(1,:,:,:,5),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"rho_vx", uin(1,:,:,:,2),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"rho_vy", uin(1,:,:,:,3),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"rho_vz", uin(1,:,:,:,4),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"Bx",     uin(1,:,:,:,6),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"By",     uin(1,:,:,:,7),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"Bz",     uin(1,:,:,:,8),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"Bxr",    uin(1,:,:,:,9),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"Byr",    uin(1,:,:,:,10),n1,n2,n3,mype)
     call dump_3d_array_h5(file_id,"Bzr",    uin(1,:,:,:,11),n1,n2,n3,mype)
#endif
  else
     call dump_4d_array_h5(file_id,"uin",uin,n1,n2,n3,nvar+3)
  endif

  ! Dataset "para_mpi"; Type integer; dim=1; size=nb_mpi
#if PARAHDF5==1
  if(use_collective_io) then
     para_mpi=(/ xleft,xright,yleft,yright,zleft,zright,&
          & xposition,yposition,zposition/)
     call dump_2d_array_int_h5(file_id, "para_mpi", para_mpi, nb_mpi, 1, mype)
  else
     para_mpi=(/ xleft,xright,yleft,yright,zleft,zright,&
          & xposition,yposition,zposition/)
     call dump_1d_array_int_h5(file_id,"para_mpi",para_mpi,nb_mpi)
  endif
#elif WITHMPI==1
  para_mpi=(/ xleft,xright,yleft,yright,zleft,zright,&
       & xposition,yposition,zposition/)
  call dump_1d_array_int_h5(file_id,"para_mpi",para_mpi,nb_mpi)
#endif
  
  ! Terminate access to HDF5
  call H5Fclose_f(file_id,ierr)
#if PARAHDF5==1  
  if(use_collective_io) call H5Pclose_f(fapl_id,ierr)
#endif
  call H5close_f(ierr)

#if WITHMPI==1
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! call cpu_time(tend)
  tend = MPI_Wtime()
  if(mype==0) write(*,*) 'WRITING OUTPUT | time elapsed: ', tend-tbegin, ' s'
#else
  if(mype==0) write(*,*) 'WRITING OUTPUT'
#endif

  return

end subroutine output
#elif PNETCDF==1
!###########################################################
!###########################################################
!###########################################################
subroutine output(uin,x,y,z,time,dt,dx,dy,dz,nhist,ndump,nspec,mype,npes)
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
  integer  :: ierr
  integer, parameter :: nb_real=10, nb_int=9, nb_mpi=9
  integer(dp) :: nbreal, nbint, nbmpi
  integer , dimension(nb_int)  :: para_int
  real(dp), dimension(nb_real) :: para_real
  real(dp), dimension(nb_mpi)  :: para_mpi
  real(dp), dimension(6) :: boxSize

  character(LEN=80) :: filename

  ! PnetCDF variables
  integer(kind=MPI_OFFSET_KIND) :: nxglob, nyglob, nzglob
  integer :: nout, ndid, ncid, xdimid, ydimid, zdimid
  integer :: niid, nrid, bsid, nmid
  integer, dimension(3) :: sdimid
  integer :: prid, piid, nsid, pmid
  integer :: rhoid, Eid, vxid, vyid, vzid, Bxid, Byid, Bzid
  integer(kind=MPI_OFFSET_KIND), dimension(3) :: dims, start, count
  integer :: diminline

  real(dp) :: tbegin, tend

  if(verbose) write(*,*) 'Entering output, using ParallelNetCDF'
  ! call cpu_time(tbegin)
  ! call MPI_Comm_Set_Attr(MPI_COMM_WORLD, MPI_WTIME_IS_GLOBAL, .true.)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tbegin = MPI_Wtime()

#if WITHMPI==1
  if(inline_io) then
     nxglob = (nx+2*nghost)
     nyglob = (ny+2*nghost)
     nzglob = (nz+2*nghost)*nxslice*nyslice*nzslice
  else
     nxglob = (nx+2*nghost)*nxslice
     nyglob = (ny+2*nghost)*nyslice
     nzglob = (nz+2*nghost)*nzslice
  endif
#else
  nxglob = nx+2*nghost
  nyglob = ny+2*nghost
  nzglob = nz+2*nghost
#endif

  dims = (/ nx+2*nghost, ny+2*nghost, nz+2*nghost /)

  call get_filename(ndump,0,'slices',filename)

  if(mype==0) then
     ! Write meta-data
     ! Create file
     nout = nfmpi_create(MPI_COMM_SELF, filename &
          , ior(NF_CLOBBER,NF_64BIT_OFFSET), MPI_INFO_NULL, ndid)
     call pcheck(nout, 'create file')

     ! Define dimensions
     nbint = nb_int; nbreal = nb_real
     nout = nfmpi_def_dim(ndid, "nb_int", nbint, niid)
     nout = nfmpi_def_dim(ndid, "nb_real", nbreal, nrid)
     nout = nfmpi_def_dim(ndid, "nb_box", 6_8, nsid)
     nout = nfmpi_def_dim(ndid, "nb_mpi", nbmpi, nmid)
     
     ! Create variables
     nout = nfmpi_def_var(ndid, "para_int", NF_INT, 1, (/niid/), piid)
     nout = nfmpi_def_var(ndid, "para_real", NF_DOUBLE, 1, (/nrid/), prid)
     nout = nfmpi_def_var(ndid, "boxSize", NF_DOUBLE, 1, (/nsid/), bsid)
     nout = nfmpi_def_var(ndid, "para_mpi", NF_DOUBLE, 1, (/nmid/), pmid)
     
     ! End of definitions
     nout = nfmpi_enddef(ndid)
     
     ! Write metadata
#if WITHMPI==1
     para_int = (/ ndump,nhist,nspec,nx,ny,nz,nxslice,nyslice,nzslice /)
     para_mpi = (/ xleft,xright,yleft,yright,zleft,zright,&
          & xposition,yposition,zposition /)
#else
     para_int = (/ ndump,nhist,nspec,nx,ny,nz,1,1,1 /)
#endif
     !para_real = (/ time,dt,dx,dy,dz /)
     para_real = (/ time, dt, dx, dy, dz, x_advc_old, x_geom_old, &
                    dv_grav_old, dv_buoy_old, dv_drag_old /)
     boxSize   = (/ xmin,xmax,ymin,ymax,zmin,zmax /)
     
     nout = nfmpi_put_vara_int_all(ndid, piid, (/1_8/), (/nbint/), para_int)
     nout = nfmpi_put_vara_double_all(ndid, prid, (/1_8/), (/nbreal/) &
          &, para_real)
     nout = nfmpi_put_vara_double_all(ndid, bsid, (/1_8/), (/6_8/), boxSize)
#if WITHMPI==1
     nout = nfmpi_put_vara_double_all(ndid, pmid, (/1_8/), (/nbmpi/), para_mpi)
#endif
     
     ! Close file
     nout = nfmpi_close(ndid)
  endif
  
  ! Synchronize the threads to ensure the metadata creation
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! Write data
  ! Open file
  nout = nfmpi_open(MPI_COMM_WORLD, filename, ior(NF_WRITE, NF_64BIT_OFFSET) &
       & ,MPI_INFO_NULL, ncid)
  nout = nfmpi_redef(ncid)

  ! Define dimensions
  nout = nfmpi_def_dim(ncid, "nx", nxglob, xdimid)
  nout = nfmpi_def_dim(ncid, "ny", nyglob, ydimid)
  nout = nfmpi_def_dim(ncid, "nz", nzglob, zdimid)
  sdimid = (/ xdimid, ydimid, zdimid /)

  ! Create variables
  nout = nfmpi_def_var(ncid, "rho", NF_DOUBLE, 3, sdimid, rhoid)
  nout = nfmpi_def_var(ncid, "E", NF_DOUBLE, 3, sdimid, Eid)
  nout = nfmpi_def_var(ncid, "rho_vx", NF_DOUBLE, 3, sdimid, vxid)
  nout = nfmpi_def_var(ncid, "rho_vy", NF_DOUBLE, 3, sdimid, vyid)
  nout = nfmpi_def_var(ncid, "rho_vz", NF_DOUBLE, 3, sdimid, vzid)
  nout = nfmpi_def_var(ncid, "Bx", NF_DOUBLE, 3, sdimid, Bxid)
  nout = nfmpi_def_var(ncid, "By", NF_DOUBLE, 3, sdimid, Byid)
  nout = nfmpi_def_var(ncid, "Bz", NF_DOUBLE, 3, sdimid, Bzid)

  ! End of definitions
  nout = nfmpi_enddef(ncid)

  ! Position of the first element and number of elements in each dimension
  if(inline_io) then
     diminline = mype*dims(3)+1 
     start = (/ 1, 1, diminline /)
  else
     start = (/ xposition, yposition, zposition /)*dims + 1
  endif
  count = dims

  ! Write data
  nout = nfmpi_put_vara_double_all(ncid, rhoid, start, count, uin(1,:,:,:,1))
  nout = nfmpi_put_vara_double_all(ncid, Eid, start, count, uin(1,:,:,:,5))
  nout = nfmpi_put_vara_double_all(ncid, vxid, start, count, uin(1,:,:,:,2))
  nout = nfmpi_put_vara_double_all(ncid, vyid, start, count, uin(1,:,:,:,3))
  nout = nfmpi_put_vara_double_all(ncid, vzid, start, count, uin(1,:,:,:,4))
  nout = nfmpi_put_vara_double_all(ncid, Bxid, start, count, uin(1,:,:,:,6))
  nout = nfmpi_put_vara_double_all(ncid, Byid, start, count, uin(1,:,:,:,7))
  nout = nfmpi_put_vara_double_all(ncid, Bzid, start, count, uin(1,:,:,:,8))

  ! Close file
  nout = nfmpi_close(ncid)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  tend = MPI_Wtime()
  ! call cpu_time(tend)
  if(mype==0) write(*,*) 'WRITING OUTPUT | time elapsed: ', tend-tbegin, ' s'

  return

end subroutine output
!###########################################################
  subroutine pcheck(errcode, indication)
    ! Check subroutine for ParallelNetCDF files
    use pnetcdf
    implicit none
    integer, intent(in) :: errcode
    character(LEN=*) :: indication
    
    if(errcode /= nf_noerr) then
       print *, indication, ' failed!'
       print *, 'Error: ', trim(nfmpi_strerror(errcode))
       stop 99
    endif
  end subroutine pcheck

#else
!###########################################################
!###########################################################
!###########################################################
subroutine output(uin,x,y,z,time,dt,dx,dy,dz,nhist,ndump,nspec,mype,npes)
   stop
#if WITHMPI==1
  use mpi_var
#endif
  use hydro_parameters
  implicit none

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) ::uin
  real(dp),dimension(iu1:iu2) ::x
  real(dp),dimension(ju1:ju2) ::y
  real(dp),dimension(ku1:ku2) ::z
  real(dp) ::time,dt,dx,dy,dz
  integer :: ndump,nhist,nspec,mype,npes

  character(LEN=80) :: filename

  if (verbose) write (*,*) 'Entering output...'

  call get_filename(ndump,mype,'slices',filename)

  open(unit=10,file=filename,status='unknown',form='unformatted')
  write(10) time,dt,dx,dy,dz
  write(10) ndump,nhist,nspec
  write(10) nx,ny,nz
#if WITHMPI==1
  write(10) nxslice,nyslice,nzslice
#else
  write(10) 1,1,1
#endif
  write(10) x,y,z
  write(10) uin
#if WITHMPI==1
  write(10) xleft,xright,yleft,yright,zleft,zright
  write(10) xposition,yposition,zposition
#endif
  close(10)

  return
end subroutine output
#endif
!###########################################################
!###########################################################
!###########################################################
subroutine get_filename(ndump,mype,prefix,filename)
  use hydro_parameters
  implicit none
#if WITHMPI==1
#include "mpif.h"
#endif

  integer :: ndump,mype

  integer :: info
  character(LEN=6) :: snumfile,n_mype
  character(LEN=*) :: prefix
  character(LEN=80) :: filename,filedir,filecmd
#if WITHMPI==1
  integer :: ierr
#endif

  call convtoasc(ndump,snumfile)
  call convtoasc(mype,n_mype  )
  filedir='output_'//trim(snumfile)//'/'
  filecmd='mkdir -p '//trim(filedir)
#ifdef NOSYSTEM
  if (mype==0) call PXFMKDIR(trim(filedir),len(trim(filedir)),O'755',info)
#elif BLUEGENE==1
  if (mype==0) call mkdir(trim(filedir)//'\0', %val(511)) ! dec(511) = oct(777)
#else
  if (mype==0) call system(filecmd)
#endif
#if WITHMPI==1
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
  filename=trim(filedir)//trim(prefix)//'.'//trim(n_mype)
  
  return
end subroutine get_filename
!###########################################################
!###########################################################
!###########################################################
subroutine convtoasc(number,sstring)
!=======================================================================
! To convert an integer smaller than 999999 to a 6 characters string
!=======================================================================
  implicit none
  integer :: number, istring, num, nums10, i
  character(LEN=6) :: sstring
  character(LEN=10),parameter :: nstring="0123456789"
     
  num=1000000
  nums10=num/10
  do i=1,6
      istring=1+mod(number,num)/nums10
      sstring(i:i)=nstring(istring:istring)
      num=num/10
      nums10=nums10/10
  enddo
  
end subroutine convtoasc
!###########################################################
!###########################################################
!###########################################################
subroutine check(string,error)
  character(len=*) :: string
  integer :: error
  if (error .lt. 0) then
     write(*,*) '######################## ', string, " FAILED"
  end if
  return
end subroutine check
