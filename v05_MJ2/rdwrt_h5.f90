module rdwrt_h5
#if HDF5==1
  use hdf5

contains

  subroutine dump_1d_array_int_h5(loc_id,dsetname,array,nx)
    use hdf5
    implicit none

    integer(hid_t)         :: loc_id
    character(LEN=*)       :: dsetname
    integer                :: nx
    integer, dimension(nx) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace
    integer(hid_t) :: h5_dset
    
    integer(hsize_t), dimension(1) :: dims
    integer :: rank=1
    
    ! Dump the array in dataset dsetname
    dims = (/ nx /)
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_INTEGER, &
         & h5_dspace, h5_dset, ierr)
    call H5Dwrite_f(h5_dset, H5T_NATIVE_INTEGER, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
    
    return
  end subroutine dump_1d_array_int_h5
  !====================================================================
  !====================================================================
  !====================================================================
  subroutine dump_1d_array_h5(loc_id,dsetname,array,nx)
    use hdf5
    implicit none
  
    integer(hid_t)        :: loc_id
    character(LEN=*)      :: dsetname
    integer               :: nx
    real*8, dimension(nx) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace
    integer(hid_t) :: h5_dset
    
    integer(hsize_t), dimension(1) :: dims
    integer :: rank=1
    
    ! Dump the array in dataset dsetname
    dims = (/ nx /)
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
         & h5_dspace, h5_dset, ierr)
    call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
    
    return
  end subroutine dump_1d_array_h5
  !====================================================================
  !====================================================================
  !====================================================================
  subroutine dump_1d_array_h5_special(loc_id,dsetname,array,nx,mype)
    use hdf5
    use hydro_parameters, only : nghost
#if WITHMPI==1
    use mpi_var, only : nxslice, xposition, yposition, zposition
#endif
    use variables, only : use_collective_io
    implicit none
    
    integer(hid_t)        :: loc_id
    character(LEN=*)      :: dsetname
    integer               :: nx,mype
    real*8, dimension(nx) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace, h5_dspace_file
    integer(hid_t) :: h5_dset
    integer(hid_t) :: plist_id
    
    integer(hsize_t), dimension(1) :: dims, dims_file
    integer(hssize_t), dimension(1) :: offset
    integer :: rank=1

    integer :: nbody, ilo, ihi
    
    ! Dump the array in dataset dsetname
    if(use_collective_io) then
#if PARAHDF5==1
       nbody = nx - 2*nghost
       if (xposition == 0) then
          dims = (/ nbody + nghost /)
          ilo = 1
          ihi = nx - nghost
       else if (xposition == nxslice - 1) then
          dims = (/ nbody + nghost /)
          ilo = 1 + nghost
          ihi = nx
       else
          dims = (/ nbody /)
          ilo = 1 + nghost
          ihi = nx - nghost
       end if
       dims_file = (/ nbody*nxslice+2*nghost /)
       if (xposition == 0) then
          offset(1) = 0
       else
          offset(1) = xposition*nbody+nghost
       endif
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dspace_file, ierr)
       call h5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
                      & h5_dspace_file, h5_dset, ierr)
       call H5Dget_space_f(h5_dset, h5_dspace_file, ierr)
       call H5Sselect_hyperslab_f(h5_dspace_file, H5S_SELECT_SET_F, offset,&
            & dims, ierr)
       call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array(ilo:ihi), dims, ierr &
            & , file_space_id = h5_dspace_file, mem_space_id = h5_dspace &
            & , xfer_prp = plist_id)
       call H5Sclose_f(h5_dspace_file, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Pclose_f(plist_id, ierr)
#endif
    else
       dims = (/ nx /)
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
            & h5_dspace, h5_dset, ierr)
       call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
    endif
    
    return
  end subroutine dump_1d_array_h5_special
  !====================================================================
  !====================================================================
  !====================================================================
  subroutine dump_2d_array_int_h5(loc_id,dsetname,array,nx,ny,mype)
    use hdf5
#if WITHMPI==1
    use mpi_var, only : nxslice,nyslice,nzslice
#endif
    use variables, only : use_collective_io
    implicit none
    
    integer(hid_t)            :: loc_id
    character(LEN=*)          :: dsetname
    integer                   :: nx,ny,mype
    integer, dimension(nx,ny) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace, h5_dspace_file
    integer(hid_t) :: h5_dset
    integer(hid_t) :: plist_id
    
    integer(hsize_t), dimension(2) :: dims, dims_file
    integer(hssize_t), dimension(2) :: offset
    integer :: rank=2
    
    ! Dump the array in dataset dsetname
    if(use_collective_io) then
#if PARAHDF5==1
       dims = (/ nx,ny /)
       dims_file = (/ nx, ny*nxslice*nyslice*nzslice /)
       offset(1) = 0
       offset(2) = mype*ny
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dspace_file, ierr)
       call h5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_INTEGER, &
                      & h5_dspace_file, h5_dset, ierr)
       call H5Dget_space_f(h5_dset, h5_dspace_file, ierr)
       call H5Sselect_hyperslab_f(h5_dspace_file, H5S_SELECT_SET_F, offset, &
            & dims, ierr)
       call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       call H5Dwrite_f(h5_dset, H5T_NATIVE_INTEGER, array, dims, ierr &
            & , file_space_id = h5_dspace_file, mem_space_id = h5_dspace &
            & , xfer_prp = plist_id)
       call H5Sclose_f(h5_dspace_file, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Pclose_f(plist_id, ierr)
#endif
    else
       dims = (/ nx,ny /)
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_INTEGER, &
            & h5_dspace, h5_dset, ierr)
       call H5Dwrite_f(h5_dset, H5T_NATIVE_INTEGER, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
    endif
    
    return
  end subroutine dump_2d_array_int_h5
  !====================================================================
  !====================================================================
  !====================================================================
  subroutine dump_2d_array_h5(loc_id,dsetname,array,nx,ny,mype)
    use hdf5
#if WITHMPI==1
    use mpi_var, only : nxslice,nyslice,nzslice
#endif
    use variables, only : use_collective_io
    implicit none
    
    integer(hid_t)           :: loc_id
    character(LEN=*)         :: dsetname
    integer                  :: nx,ny,mype
    real*8, dimension(nx,ny) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace, h5_dspace_file
    integer(hid_t) :: h5_dset
    integer(hid_t) :: plist_id
    
    integer(hsize_t), dimension(2) :: dims, dims_file
    integer(hssize_t), dimension(2) :: offset
    integer :: rank=2
    
    ! Dump the array in dataset dsetname
    if(use_collective_io) then
#if PARAHDF5==1
       dims = (/ nx,ny /)
       dims_file = (/ nx, ny*nxslice*nyslice*nzslice /)
       offset(1) = 0
       offset(2) = mype*ny
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dspace_file, ierr)
       call h5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
                      & h5_dspace_file, h5_dset, ierr)
       call H5Dget_space_f(h5_dset, h5_dspace_file, ierr)
       call H5Sselect_hyperslab_f(h5_dspace_file, H5S_SELECT_SET_F, offset, &
            & dims, ierr)
       call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr &
            & , file_space_id = h5_dspace_file, mem_space_id = h5_dspace &
            & , xfer_prp = plist_id)
       call H5Sclose_f(h5_dspace_file, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Pclose_f(plist_id, ierr)
#endif
    else
       dims = (/ nx,ny /)
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
            & h5_dspace, h5_dset, ierr)
       call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
    endif
    
    return
  end subroutine dump_2d_array_h5
  !====================================================================
  !====================================================================
  !====================================================================
  subroutine dump_3d_array_h5(loc_id,dsetname,array,n1,n2,n3,mype)
    use hdf5
    use variables,      only : use_collective_io, store_ghost, inline_io
    use amr_parameters, only : nx,ny,nz
    use hydro_parameters, only : nghost
#if WITHMPI==1
    use mpi_var,        only : nxslice,nyslice,nzslice, xposition, yposition, zposition
#endif

    implicit none
    
    integer(hid_t)              :: loc_id
    character(LEN=*)            :: dsetname
    integer                     :: n1,n2,n3,mype
    integer                     :: nx_glob,ny_glob,nz_glob
    real*8, dimension(n1,n2,n3) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace
    integer(hid_t) :: h5_dset
    
    ! the following parameters are only useful for the collective output method
    integer(hid_t)                 :: h5_dspace_file
    integer(hsize_t), dimension(3) :: start ! Offset of start of hyperslab
    integer(hsize_t), dimension(3) :: count ! Number of blocks to select in data space
    integer(hsize_t), dimension(3) :: stride
    integer(hsize_t), dimension(3) :: blockSize
    integer(hid_t)                 :: prop_data_chunk_id
    integer(hid_t)                 :: dxpl_id  ! dataset transfer property list

    ! other variables
    integer(hsize_t), dimension(3) :: dims, dims_file, dims_chunk
    integer :: rank=3
    
    ! This version contains the ghostzones
    dims = (/ n1,n2,n3 /)

#if WITHMPI==1
    if(inline_io) then
       nx_glob = (nx+2*nghost)
       ny_glob = (ny+2*nghost)
       nz_glob = (nz+2*nghost)*nxslice*nyslice*nzslice
    else if(store_ghost) then
       nx_glob = (nx+2*nghost)*nxslice
       ny_glob = (ny+2*nghost)*nyslice
       nz_glob = (nz+2*nghost)*nzslice
    else
       nx_glob = nx*nxslice+2*nghost
       ny_glob = ny*nyslice+2*nghost
       nz_glob = nz*nzslice+2*nghost
    end if
#else
    nx_glob = nx+2*nghost
    ny_glob = ny+2*nghost
    nz_glob = nz+2*nghost
#endif    

#if WITHMPI==1
    if(use_collective_io) then
#if PARAHDF5==1
       ! one single hdf5 file
       
       ! dims DOES contain ghostzones
       dims_file = (/ nx_glob, ny_glob, nz_glob /)
       
       call H5Screate_simple_f(rank, dims     , h5_dspace     , ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dspace_file, ierr)
       
       ! hyperslab for selecting data in h5_dspace WITH ghostcells  
       start     = (/ 0,0,0 /)
       stride    = (/ 1,1,1 /)
       count     = dims
       blockSize = (/ 1,1,1 /)
       call H5Sselect_hyperslab_f(h5_dspace, H5S_SELECT_SET_F, start,&
            count, ierr, stride, blockSize)
       
       ! hyperslab for selecting location in h5_dspace_file (to set the
       ! correct location in file where we want to put our piece of data)
       if(inline_io) then
          start(1) = 0
          start(2) = 0
          start(3) = dims(3)*mype
       else if(store_ghost) then
          start  = (/ xposition, yposition, zposition /)*dims
       else
          start  = (/ xposition, yposition, zposition /)*(dims-2*nghost)
       end if
       stride    = (/ 1,1,1 /)
       count     = dims
       blockSize = (/ 1,1,1 /)
       call H5Sselect_hyperslab_f(h5_dspace_file, H5S_SELECT_SET_F, start,&
            count, ierr, stride, blockSize)
       
       ! create property id for our data chunk
       dims_chunk = dims
       call H5Pcreate_f (H5P_DATASET_CREATE_F, prop_data_chunk_id, ierr)
       call H5Pset_chunk_f(prop_data_chunk_id, 3, dims_chunk, ierr)
       
       ! enable parallel collective IO
       call H5Pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
       call check("h5pcreate_f", ierr)
       call H5Pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F,ierr)
       call check("H5Pset_dxpl_mpio_f", ierr)
       
       ! create data set
       call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE,&
            & h5_dspace_file, h5_dset, ierr, prop_data_chunk_id) 
       
       ! finally write data to file
       call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims,&
            & ierr, mem_space_id=h5_dspace, file_space_id=h5_dspace_file, xfer_prp=dxpl_id)
       call check("H5Dwrite_f",ierr)
       
       ! clean open hdf5 Ids
       call H5Pclose_f(prop_data_chunk_id, ierr)
       call H5Pclose_f(dxpl_id, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Sclose_f(h5_dspace_file, ierr)
#endif
    else
       ! one HDF5 file per MPI process
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
            & h5_dspace, h5_dset, ierr)
       call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)          
    endif
#else
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
         & h5_dspace, h5_dset, ierr)
    call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
#endif
    
    return
  end subroutine dump_3d_array_h5
  !====================================================================
  !====================================================================
  !====================================================================
  subroutine dump_4d_array_h5(loc_id,dsetname,array,n1,n2,n3,n4)
    use hdf5
    implicit none
    
    integer(hid_t)                 :: loc_id
    character(LEN=*)               :: dsetname
    integer                        :: n1,n2,n3,n4
    real*8, dimension(n1,n2,n3,n4) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace
    integer(hid_t) :: h5_dset
    
    integer(hsize_t), dimension(4) :: dims
    integer :: rank=4
    
    ! Dump the array in dataset dsetname
    dims = (/ n1,n2,n3,n4 /)
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
         & h5_dspace, h5_dset, ierr)
    call H5Dwrite_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
    
    return
  end subroutine dump_4d_array_h5
  !===========================================================================
  !===========================================================================
  !===========================================================================
  subroutine get_1d_array_int_h5(loc_id,dsetname,array,nx)
    use hdf5
    implicit none
    
    integer(hid_t)               :: loc_id
    character(LEN=*)             :: dsetname
    integer                      :: nx
    integer, dimension(nx)       :: array
  
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dset
    
    integer(hsize_t), dimension(1) :: dims
    
    dims=(/ nx /)
    call H5Dopen_f(loc_id,trim(dsetname),h5_dset,ierr)
    call H5Dread_f(h5_dset,H5T_NATIVE_INTEGER,array,dims,ierr)
    call H5Dclose_f(h5_dset,ierr)
    
    return
  end subroutine get_1d_array_int_h5
  !===========================================================================
  !===========================================================================
  !===========================================================================
  subroutine get_1d_array_h5(loc_id,dsetname,array,nx)
    use hdf5
    implicit none
    
    integer(hid_t)              :: loc_id
    character(LEN=*)            :: dsetname
    integer                     :: nx
    real*8, dimension(nx)       :: array
  
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dset
    
    integer(hsize_t), dimension(1) :: dims
    
    dims=(/ nx /)
    call H5Dopen_f(loc_id,trim(dsetname),h5_dset,ierr)
    call H5Dread_f(h5_dset,H5T_NATIVE_DOUBLE,array,dims,ierr)
    call H5Dclose_f(h5_dset,ierr)
    
    return
  end subroutine get_1d_array_h5
  !====================================================================
  !====================================================================
  !====================================================================
  subroutine get_1d_array_h5_special(loc_id,dsetname,array,nx,mype)
    use hdf5
    use hydro_parameters, only : nghost
#if WITHMPI==1
    use mpi_var, only : nxslice, xposition, yposition, zposition
#endif
    use variables, only : use_collective_io
    implicit none
    
    integer(hid_t)        :: loc_id
    character(LEN=*)      :: dsetname
    integer               :: nx,mype
    real*8, dimension(nx) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace, h5_dspace_file
    integer(hid_t) :: h5_dset
    integer(hid_t) :: plist_id
    
    integer(hsize_t), dimension(1) :: dims, dims_file
    integer(hssize_t), dimension(1) :: offset
    integer :: rank=1

    integer :: nbody, ilo, ihi
    
    ! Dump the array in dataset dsetname
    if(use_collective_io) then
#if PARAHDF5==1
       nbody = nx - 2*nghost
       dims = (/ nx /)
       ilo = 1
       ihi = nx
       nbody = nx - 2*nghost
       dims_file = (/ nbody*nxslice+2*nghost /)
       offset(1) = xposition*nbody
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dspace_file, ierr)
       call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
       call H5Dget_space_f(h5_dset, h5_dspace_file, ierr)
       call H5Sselect_hyperslab_f(h5_dspace_file, H5S_SELECT_SET_F, offset, &
            & dims, ierr)
       call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array(ilo:ihi), dims, ierr &
            & , file_space_id = h5_dspace_file, mem_space_id = h5_dspace &
            & , xfer_prp = plist_id)
       call H5Sclose_f(h5_dspace_file, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Pclose_f(plist_id, ierr)
#endif
    else
       dims = (/ nx /)
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
       call H5Dread_f(h5_dset, H5T_NATIVE_INTEGER, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
    endif
    
    return
  end subroutine get_1d_array_h5_special
  !====================================================================
  !====================================================================
  !====================================================================
  subroutine get_2d_array_int_h5(loc_id,dsetname,array,nx,ny,mype)
    use hdf5
#if WITHMPI==1
    use mpi_var, only : nxslice,nyslice,nzslice
#endif
    use variables, only : use_collective_io
    implicit none
    
    integer(hid_t)            :: loc_id
    character(LEN=*)          :: dsetname
    integer                   :: nx,ny,mype
    integer, dimension(nx,ny) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace, h5_dspace_file
    integer(hid_t) :: h5_dset
    integer(hid_t) :: plist_id
    
    integer(hsize_t), dimension(2) :: dims, dims_file
    integer(hssize_t), dimension(2) :: offset
    integer :: rank=2
    
    ! Dump the array in dataset dsetname
    if(use_collective_io) then
#if PARAHDF5==1
       dims = (/ nx,ny /)
       dims_file = (/ nx, ny*nxslice*nyslice*nzslice /)
       offset(1) = 0
       offset(2) = mype*ny
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dspace_file, ierr)
       call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
       call H5Dget_space_f(h5_dset, h5_dspace_file, ierr)
       call H5Sselect_hyperslab_f(h5_dspace_file, H5S_SELECT_SET_F, offset, &
            & dims, ierr)
       call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       call H5Dread_f(h5_dset, H5T_NATIVE_INTEGER, array, dims, ierr &
            & , file_space_id = h5_dspace_file, mem_space_id = h5_dspace &
            & , xfer_prp = plist_id)
       call H5Sclose_f(h5_dspace_file, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Pclose_f(plist_id, ierr)
#endif
    else
       dims = (/ nx,ny /)
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
       call H5Dread_f(h5_dset, H5T_NATIVE_INTEGER, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
endif
    
    return
  end subroutine get_2d_array_int_h5
  !====================================================================
  !====================================================================
  !====================================================================
  subroutine get_2d_array_h5(loc_id,dsetname,array,nx,ny,mype)
    use hdf5
#if WITHMPI==1
    use mpi_var, only : nxslice,nyslice,nzslice
#endif
    use variables, only : use_collective_io
    implicit none
    
    integer(hid_t)           :: loc_id
    character(LEN=*)         :: dsetname
    integer                  :: nx,ny,mype
    real*8, dimension(nx,ny) :: array
    
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace, h5_dspace_file
    integer(hid_t) :: h5_dset
    integer(hid_t) :: plist_id
    
    integer(hsize_t), dimension(2) :: dims, dims_file
    integer(hssize_t), dimension(2) :: offset
    integer :: rank=2
    
    ! Dump the array in dataset dsetname
    if(use_collective_io) then
#if PARAHDF5==1
       dims = (/ nx,ny /)
       dims_file = (/ nx, ny*nxslice*nyslice*nzslice /)
       offset(1) = 0
       offset(2) = mype*ny
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dspace_file, ierr)
       call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
       call H5Dget_space_f(h5_dset, h5_dspace_file, ierr)
       call H5Sselect_hyperslab_f(h5_dspace_file, H5S_SELECT_SET_F, offset, &
            & dims, ierr)
       call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
       call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr &
            & , file_space_id = h5_dspace_file, mem_space_id = h5_dspace &
            & , xfer_prp = plist_id)
       call H5Sclose_f(h5_dspace_file, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Pclose_f(plist_id, ierr)
#endif
    else
       dims = (/ nx,ny /)
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
       call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
endif
    
    return
  end subroutine get_2d_array_h5
  !===========================================================================
  !===========================================================================
  !===========================================================================
  subroutine get_3d_array_h5(loc_id,dsetname,array,n1,n2,n3,mype)
    use hdf5
    use variables,      only : use_collective_io, store_ghost, inline_io
    use amr_parameters, only : nx,ny,nz
    use hydro_parameters, only : nghost
#if WITHMPI==1
    use mpi_var,        only : nxslice,nyslice,nzslice, xposition, yposition, zposition
#endif

    implicit none
    
    integer(hid_t)              :: loc_id
    character(LEN=*)            :: dsetname
    integer                     :: n1,n2,n3,mype
    integer                     :: nx_glob,ny_glob,nz_glob
    real*8, dimension(n1,n2,n3) :: array
  
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dspace
    integer(hid_t) :: h5_dset
    
    ! the following parameters are only useful for the collective output method
    integer(hid_t)                 :: h5_dspace_file
    integer(hsize_t), dimension(3) :: start ! Offset of start of hyperslab
    integer(hsize_t), dimension(3) :: count ! Number of blocks to select in data space
    integer(hsize_t), dimension(3) :: stride
    integer(hsize_t), dimension(3) :: blockSize
    integer(hid_t)                 :: dxpl_id ,dcpl  ! dataset transfer property list

    ! other variables
    integer(hsize_t), dimension(3) :: dims, dims_file, dims_chunk
    integer :: rank=3

    ! This version contains the ghostzones
    dims = (/ n1,n2,n3 /)

#if WITHMPI==1
    if(inline_io) then
       nx_glob = (nx+2*nghost)
       ny_glob = (ny+2*nghost)
       nz_glob = (nz+2*nghost)*nxslice*nyslice*nzslice
    else if(store_ghost) then
       nx_glob = (nx+2*nghost)*nxslice
       ny_glob = (ny+2*nghost)*nyslice
       nz_glob = (nz+2*nghost)*nzslice
    else
       nx_glob = nx*nxslice+2*nghost
       ny_glob = ny*nyslice+2*nghost
       nz_glob = nz*nzslice+2*nghost
    end if
#else
    nx_glob = nx+2*nghost              
    ny_glob = ny+2*nghost              
    nz_glob = nz+2*nghost              
#endif    

#if WITHMPI==1
    if(use_collective_io) then
#if PARAHDF5==1
       ! one single hdf5 file
       
       ! dims DOES contain ghostzones
       dims_file = (/ nx_glob, ny_glob, nz_glob /)
       
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Screate_simple_f(rank, dims_file, h5_dspace_file, ierr)
       
       ! create data set
       call H5Dopen_f(loc_id, trim(dsetname), h5_dset, ierr)
       call check("h5dopen_f", ierr)
       
       ! hyperslab for selecting data in h5_dspace WITH ghostcells  
       start     = (/ 0,0,0 /)
       stride    = (/ 1,1,1 /)
       count     = dims
       blockSize = (/ 1,1,1 /)
       call H5Sselect_hyperslab_f(h5_dspace, H5S_SELECT_SET_F, start,&
            count, ierr, stride, blockSize)
       
       ! hyperslab for selecting location in h5_dspace_file (to set the
       ! correct location in file where we want to put our piece of data)
       if(inline_io) then
          start(1) = 0
          start(2) = 0
          start(3) = dims(3)*mype
       else if(store_ghost) then
          start  = (/ xposition, yposition, zposition /)*dims
       else
          start  = (/ xposition, yposition, zposition /)*(dims-2*nghost)
       end if
       stride    = (/ 1,1,1 /)
       count     = dims
       blockSize = (/ 1,1,1 /)
       call H5Sselect_hyperslab_f(h5_dspace_file, H5S_SELECT_SET_F, start,&
            count, ierr, stride, blockSize)
       
       ! enable parallel collective IO
       call H5Pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
       call check("h5pcreate_f", ierr)
       call H5Pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F,ierr)
       call check("H5Pget_dxpl_mpio_f", ierr)
       
       ! finally read data to file
       call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims,&
            & ierr, mem_space_id=h5_dspace, file_space_id=h5_dspace_file,&
            & xfer_prp=dxpl_id)
       call check("H5Dread_f",ierr)
       
       ! clean open hdf5 Ids
       ! call H5Pclose_f(prop_data_chunk_id, ierr)
       call H5Pclose_f(dxpl_id, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)
       call H5Sclose_f(h5_dspace_file, ierr)
#endif
    else
       ! one HDF5 file per MPI process
       call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
       call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
            & h5_dspace, h5_dset, ierr)
       call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
       call H5Dclose_f(h5_dset, ierr)
       call H5Sclose_f(h5_dspace, ierr)          
    endif
#else
    call H5Screate_simple_f(rank, dims, h5_dspace, ierr)
    call H5Dcreate_f(loc_id, trim(dsetname), H5T_NATIVE_DOUBLE, &
         & h5_dspace, h5_dset, ierr)
    call H5Dread_f(h5_dset, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call H5Dclose_f(h5_dset, ierr)
    call H5Sclose_f(h5_dspace, ierr)
#endif
    
    return
  end subroutine get_3d_array_h5
  !===========================================================================
  !===========================================================================
  !===========================================================================
  subroutine get_4d_array_h5(loc_id,dsetname,array,n1,n2,n3,n4)
    use hdf5
    implicit none
    
    integer(hid_t)                 :: loc_id
    character(LEN=*)               :: dsetname
    integer                        :: n1,n2,n3,n4
    real*8, dimension(n1,n2,n3,n4) :: array
  
    ! HDF5 vars
    integer        :: ierr
    integer(hid_t) :: h5_dset
    
    integer(hsize_t), dimension(4) :: dims
    
    dims=(/ n1,n2,n3,n4 /)
    call H5Dopen_f(loc_id,trim(dsetname),h5_dset,ierr)
    call H5Dread_f(h5_dset,H5T_NATIVE_DOUBLE,array,dims,ierr)
    call H5Dclose_f(h5_dset,ierr)
    
    return
  end subroutine get_4d_array_h5
  !
#endif
end module rdwrt_h5
