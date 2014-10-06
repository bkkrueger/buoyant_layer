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
  use perturb, only : x_mass_old, x_advc_old, time_old, &
                      dv_grav_old, dv_buoy_old, dv_drag_old, osc_omega
  use, intrinsic :: iso_fortran_env, only : output_unit
  implicit none
#if WITHMPI==1
  include 'mpif.h'
#endif
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

  integer, parameter :: nb_real=12,nb_int=9,nb_mpi=9

  integer , dimension(nb_int ) :: para_int
  integer , dimension(nb_mpi ) :: para_mpi
  real(dp), dimension(nb_real) :: para_real

  integer :: n1,n2,n3

  integer :: commx, commyz, root
  integer(hid_t) :: fapl_id
  real(dp), dimension(:), allocatable :: mpi_buffer
  integer :: i, v, offset, idx

  character (len=*), parameter :: f_str = "('mype = ',i3,', pos = (',i3,',',i3,',',i3,'), root = ',i3)"

  if (verbose) write (*,*) 'Entering newstart...'

#if PARAHDF5==1
     call get_dumpname(ndump,0,filename)
#else
     call get_dumpname(ndump,mype,filename)
#endif

  call H5open_f(ierr)
  call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
  call H5Pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
  call H5Fopen_f(trim(filename),H5F_ACC_RDONLY_F,file_id,ierr, &
                 access_prp=fapl_id)

  !Dataset "para_real"; Type real; dim=1; size=nb_real
  call get_1d_array_h5(file_id,"para_real",para_real,nb_real)
  time=para_real(1)
  dt  =para_real(2)
  dx  =para_real(3)
  dy  =para_real(4)
  dz  =para_real(5)
  x_mass_old  = para_real(6)
  x_advc_old  = para_real(7)
  time_old    = para_real(8)
  dv_grav_old = para_real(9)
  dv_buoy_old = para_real(10)
  dv_drag_old = para_real(11)
  osc_omega   = para_real(12)

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
  
  ! Steady state (non-collective) - - - - - - - - - - - - -
  if (.not. use_collective_io) then
     call get_4d_array_h5(file_id, "ss", ss0, n1, SS_NVAR, 1, 1)
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Close the file
  call H5Fclose_f(file_id,ierr)
  call H5Pclose_f(fapl_id, ierr)

  ! Steady state (collective) - - - - - - - - - - - - - - -
  if (use_collective_io) then
#if PARAHDF5==1
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     ! Create new communicator: columns along x (:,y,z lines)
     call create_comm_line(commx, 1, mype)
     ! Only a single line will perform this read
     if (yposition == 0 .and. zposition == 0) then
        ! Create the property list
        call H5Pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
        ! Attach line communicator to property list
        call H5Pset_fapl_mpio_f(fapl_id, commx, MPI_INFO_NULL, ierr)
        ! Open the file
        call H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr, &
                       access_prp=fapl_id)
        ! Read the data
        call get_1d_array_h5_special(file_id, "ss_dens", ss0(:,1), n1, mype)
        call get_1d_array_h5_special(file_id, "ss_momx", ss0(:,2), n1, mype)
        call get_1d_array_h5_special(file_id, "ss_ener", ss0(:,3), n1, mype)
        ! Close the file
        call H5Fclose_f(file_id, ierr)
        ! Close the property list
        call H5Pclose_f(fapl_id, ierr)
     end if
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     ! Release the line communicator
     call MPI_Comm_Free(commx, ierr)

     ! Create new communicator: planes of constant x (x,:,: planes)
     call create_comm_plane(commyz, root, 1, mype)
     ! Pack the send buffer
     allocate(mpi_buffer((iu2-iu1+1)*SS_NVAR))
     if (yposition == 0 .and. zposition == 0) then
        do i = iu1, iu2
           offset = (i - iu1) * SS_NVAR
           do v = 1, SS_NVAR
              idx = offset + v
              mpi_buffer(idx) = ss0(i,v)
           end do
        end do
     endif
     ! Broadcast across planes
     call MPI_Bcast(mpi_buffer, (iu2-iu1+1)*SS_NVAR, &
                    MPI_DOUBLE_PRECISION, root, commyz, ierr)
     ! Unpack the receive buffer
     if (yposition /= 0 .or. zposition /= 0) then
        do i = iu1, iu2
           offset = (i - iu1) * SS_NVAR
           do v = 1, SS_NVAR
              idx = offset + v
              ss0(i, v) = mpi_buffer(idx)
           end do
        end do
     end if
     deallocate(mpi_buffer)
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     ! Release the plane communicator
     call MPI_Comm_Free(commyz, ierr)
#endif
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Terminate access to HDF5
  call H5close_f(ierr)

  return
end subroutine restart_run
#elif PNETCDF==1
!###########################################################
!###########################################################
!###########################################################
subroutine restart_run(uin,x,y,z,time,dt,dx,dy,dz,ndump,nhist,nspec,mype)
  use hydro_parameters
  use pnetcdf
#if WITHMPI==1
  use mpi_var
#endif
  use variables, only : inline_io
  use steady_state, only : ss0, SS_NVAR
  use perturb, only : x_mass_old, x_advc_old, time_old, &
                      dv_grav_old, dv_buoy_old, dv_drag_old, &
                      osc_omega

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
  integer, parameter :: nb_real=12, nb_int=9, nb_mpi=9
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

  ! For steady state
  integer :: i_ss, offset, idx, v
  integer :: commx, commyz, root, ierr
  real(dp), dimension(:), allocatable :: mpi_buffer
  integer(kind=MPI_OFFSET_KIND) :: nxglob, temp
  integer(kind=MPI_OFFSET_KIND), dimension(2) :: startss, countss

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
           x_mass_old  = para_real(6)
           x_advc_old  = para_real(7)
           time_old    = para_real(8)
           dv_grav_old = para_real(9)
           dv_buoy_old = para_real(10)
           dv_drag_old = para_real(11)
           osc_omega   = para_real(12)
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
! #if WITHMPI==1
!         if(trim(name).eq.'para_mpi')then
!            nbmpi = nb_mpi
!            nout = nfmpi_get_vara_double_all(ncid, i, (/1_8/), (/nbmpi/) &
!                 , para_real)
!            xleft     = para_mpi(1)
!            xright    = para_mpi(2)
!            yleft     = para_mpi(3)
!            yright    = para_mpi(4)
!            zleft     = para_mpi(5)
!            zright    = para_mpi(6)
!            xposition = para_mpi(7)
!            yposition = para_mpi(8)
!            zposition = para_mpi(9)
!         endif
! #endif
     endif
  enddo

  ! Synchronize all threads to ensure metadata reading
  call MPI_Barrier(MPI_COMM_WORLD, nout)

  dims = (/ nx+2*nghost, ny+2*nghost, nz+2*nghost /)

  ! Loop over variables to retrieve data
  do i = 1, nvars
     nout = nfmpi_inq_var(ncid, i, name, type, ndimvar, dimid, natts)
     
     if (trim(name).eq.'rho') then
        id = 1
     else if (trim(name).eq.'E') then
        id = 5
     else if (trim(name).eq.'rho_vx') then
        id = 2
     else if (trim(name).eq.'rho_vy') then
        id = 3
     else if (trim(name).eq.'rho_vz') then
        id = 4
     else if (trim(name).eq.'Bx') then
        id = 6
     else if (trim(name).eq.'By') then
        id = 7
     else if (trim(name).eq.'Bz') then
        id = 8
     else if (trim(name).eq.'Bxr') then
        id = 9
     else if (trim(name).eq.'Byr') then
        id = 10
     else if (trim(name).eq.'Bzr') then
        id = 11
     else if (trim(name).eq.'ss0') then
        i_ss = i
     else
        continue
     end if
     
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

  ! Open the file and extract the steady state along the (x,0,0) line
  nxglob = (nx+2*nghost)
  call create_comm_line(commx, 1, mype)
  if (yposition == 0 .and. zposition == 0) then
     ! Open netCDF file
     nout = nfmpi_open(commx, filename, ior(NF_NOWRITE, NF_64BIT_OFFSET) &
          & ,MPI_INFO_NULL, ncid)
     ! Define start and count
     temp = SS_NVAR
     countss = (/ nxglob, temp /)
     startss = xposition * countss + 1
     ! Grab the data
     nout = nfmpi_get_vara_double_all(ncid, i_ss, startss, countss, ss0(:,:))
     ! Close the file
     nout = nfmpi_close(ncid)
  end if
  call MPI_Barrier(commx, ierr)
  call MPI_Comm_Free(commx,ierr)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! Broadcast along yz slices: all processors with same xposition have same ss0
  call create_comm_plane(commyz, root, 1, mype)
  ! Pack the send buffer
  allocate(mpi_buffer((iu2-iu1+1)*SS_NVAR))
  if (yposition == 0 .and. zposition == 0) then
     do i = iu1, iu2
        offset = (i - iu1) * SS_NVAR
        do v = 1, SS_NVAR
           idx = offset + v
           mpi_buffer(idx) = ss0(i,v)
        end do
     end do
  endif
  ! Broadcast
  call MPI_Bcast(mpi_buffer, (iu2-iu1+1)*SS_NVAR, &
                 MPI_DOUBLE_PRECISION, root, commyz, ierr)
  ! Unpack the receive buffer
  if (yposition /= 0 .or. zposition /= 0) then
     do i = iu1, iu2
        offset = (i - iu1) * SS_NVAR
        do v = 1, SS_NVAR
           idx = offset + v
           ss0(i, v) = mpi_buffer(idx)
        end do
     end do
  end if
  deallocate(mpi_buffer)
  call MPI_Barrier(commyz, ierr)
  call MPI_Comm_Free(commyz,ierr)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

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
