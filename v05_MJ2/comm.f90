! Written by Marc Joos -- creates communicators either in a plane (all
! processors that share the specified coordinate) or in a line (all processors
! that vary only along the specified coordinate).  I think; check the routines
! to be sure.  --BKK

!===============================================================================
!> This routine creates communicators in process planes
!===============================================================================
subroutine create_comm_plane(comm, zero, direction, mype)
#if WITHMPI == 1
  use hydro_parameters
  use mpi_var
  implicit none
#if WITHMPI==1
  include 'mpif.h'
#endif

  integer :: direction, mype
  integer, intent(out) :: comm, zero
  integer :: tmp_comm
  integer :: base_grp, grp_comm, ierr
  integer :: N, nslice, posn, p0
  integer, allocatable, dimension(:) :: list, list2
  integer :: i, j, id

  ! Create a master group to create other groups
  call MPI_Comm_Group(MPI_COMM_WORLD, base_grp, ierr)

  if (direction == 1) then
     N      = nyslice*nzslice
     nslice = nxslice
  else if (direction == 2) then
     N      = nxslice*nzslice
     nslice = nyslice
  else
     N      = nxslice*nyslice
     nslice = nzslice
  endif

  allocate(list(N))
  allocate(list2(N))

  do j = 0, nslice-1
     id = 1
     do i = 0, nxslice*nyslice*nzslice-1
        if (direction == 1) posn = mod(i, nxslice)
        if (direction == 2) posn = mod(i/(nxslice*nzslice), nyslice)
        if (direction == 3) posn = mod(i/nxslice, nzslice)
        if (posn == j) then
           list(id) = i
           id       = id + 1
        endif
     enddo

     ! Create a group per slice
     call MPI_Group_Incl(base_grp, N, list, grp_comm, ierr)
     call MPI_Comm_Create(MPI_COMM_WORLD, grp_comm, tmp_comm, ierr)


     ! Determine the current processor information
     if (direction == 1) then
        posn = xposition
        p0 = xposition
     else if (direction == 2) then
        posn = yposition
        p0 = yposition * nxslice * nzslice
     else
        posn = zposition
        p0 = zposition * nxslice
     end if
     if (posn == j) then
        comm = tmp_comm
        call MPI_Group_translate_ranks(base_grp,N,list, grp_comm,list2, ierr)
        do i = 1, N
           if (list(i) == p0) then
              zero = list2(i)
           end if
        end do
     endif

     ! Free the new group
     call MPI_Group_Free(grp_comm, ierr)
  enddo

  call MPI_Group_Free(base_grp, ierr)
  deallocate(list)

#endif
  return
end subroutine create_comm_plane
!===============================================================================
!> This routine creates communicators in a given direction
!===============================================================================
subroutine create_comm_line(comm, direction, mype)
#if WITHMPI == 1
  use hydro_parameters
  use mpi_var
  implicit none
#if WITHMPI==1
  include 'mpif.h'
#endif

  integer :: direction, mype
  integer, intent(out) :: comm
  integer :: tmp_comm
  integer :: base_grp, grp_comm, ierr
  integer :: dim, n1slice, n2slice, position1, position2
  integer, allocatable, dimension(:) :: list
  integer :: i, j, k, id

  ! Create a master group to create other groups
  call MPI_Comm_Group(MPI_COMM_WORLD, base_grp, ierr)

  if (direction == 1) then
     dim     = nxslice
     n1slice = nyslice
     n2slice = nzslice
  else if (direction == 2) then
     dim     = nyslice
     n1slice = nxslice
     n2slice = nzslice
  else
     dim     = nzslice
     n1slice = nxslice
     n2slice = nyslice
  endif

  allocate(list(dim))

  do k = 0, n2slice-1
     do j = 0, n1slice-1
        id = 1
        do i = 0, nxslice*nyslice*nzslice-1
           if (direction == 1) then
              position1 = mod(i/(nxslice*nzslice), nyslice)
              position2 = mod(i/nxslice, nzslice)
           else if (direction == 2) then
              position1 = mod(i, nxslice)
              position2 = mod(i/nxslice, nzslice)
           else
              position1 = mod(i, nxslice)
              position2 = mod(i/(nxslice*nzslice), nyslice)
           endif
           if (position1 == j .and. position2 == k) then
              list(id) = i
              id       = id + 1
           endif
        enddo
        ! Create a group per line
        call MPI_Group_Incl(base_grp, dim, list, grp_comm, ierr)
        call MPI_Comm_Create(MPI_COMM_WORLD, grp_comm, tmp_comm, ierr)
        ! Determine the group of the current thread
        if (direction == 1) then
           position1 = mod(mype/(nxslice*nzslice), nyslice)
           position2 = mod(mype/nxslice, nzslice)
        else if (direction == 2) then
           position1 = mod(mype, nxslice)
           position2 = mod(mype/nxslice, nzslice)
        else
           position1 = mod(mype, nxslice)
           position2 = mod(mype/(nxslice*nzslice), nyslice)
        endif
        if (position1 == j .and. position2 == k) comm = tmp_comm
        call MPI_Group_Free(grp_comm, ierr)
     enddo
  enddo

  call MPI_Group_Free(base_grp, ierr)
  deallocate(list)

#endif
end subroutine create_comm_line

