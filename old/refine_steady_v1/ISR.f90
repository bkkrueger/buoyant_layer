
! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

module initial_state_refinement

   use amr_parameters, only : dp
   use steady, only : SS_nv

   implicit none

   integer :: ISR_N                 ! Number of elements
   integer :: glb_mype, glb_npes    ! Proc ID and # procs in MPI_COMM_WORLD
#  if WITHMPI==1
   integer :: ISR_comm              ! Communicator for the vertical slice
   integer :: ISR_mype, ISR_npes    ! Proc ID and # procs in vertical slice
   integer :: ISR_left, ISR_right   ! Neighbors in vertical slice
#  endif
   real(dp) :: ss_bc(SS_nv)         ! Upper BCs for ss0

   contains

   ! ==========================================================================

   subroutine ISR_setup(mype, npes)

      use amr_parameters, only : nx
      use heating_layer, only : Mup
      use hydro_parameters, only : iu1, iu2, gamma
      use steady, only : SS_nv
#     if WITHMPI==1
      use mpi_var, only : nxslice,   nyslice,   nzslice,    &
                        & xposition, yposition, zposition
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! Declare variables -----------------------------------------------------
      integer, intent(in) :: mype, npes
#     if WITHMPI==1
      integer :: mpi_color
      integer :: ierr
#     endif

      ! Set up ----------------------------------------------------------------

      ! Save upper BCs for ss0
      ss_bc(1) = 1.0d0
      ss_bc(2) = -Mup
      ss_bc(3) = 1.0d0/(gamma*(gamma-1.0d0)) + 0.5d0*Mup**2

#     if WITHMPI==1
      ! MPI VERSION - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ISR_N = nx * nxslice * SS_nv

      glb_mype = mype
      glb_npes = npes

      ! Define partition identifier: each vertical slice constitutes a single
      ! communicator
      mpi_color = yposition*nzslice + zposition

      ! Partition
      call MPI_Comm_split(MPI_COMM_WORLD, mpi_color, xposition, ISR_comm, ierr)
      if ((ierr /= MPI_SUCCESS) .or. (ISR_comm == MPI_COMM_NULL)) then
         call ISR_stop("ERROR: Failure in MPI_Comm_split")
      end if

      ! Get the size of the new communicator
      call MPI_Comm_size(ISR_comm, ISR_npes, ierr)
      if (ISR_npes /= nxslice) then
         call ISR_stop("ERROR: Unexpected size for ISR_comm")
      end if

      ! Get the rank within the new communicator
      call MPI_Comm_rank(ISR_comm, ISR_mype, ierr)
      if (ISR_mype /= xposition) then
         call ISR_stop("ERROR: Unexpected rank in ISR_comm")
      end if

      ! Define neighbors
      if (ISR_mype == 0) then
         ISR_left = ISR_npes - 1
      else
         ISR_left = ISR_mype - 1
      end if
      if (ISR_mype == ISR_npes-1) then
         ISR_right = 0
      else
         ISR_right = ISR_mype + 1
      end if

#     else
      ! NON-MPI VERSION - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ISR_N = nx * SS_nv
#     endif

      return
   end subroutine ISR_setup

   ! ==========================================================================

   subroutine ISR_cleanup ()

      implicit none

#     if WITHMPI==1
#     include "mpif.h"

      integer :: ierr
      call MPI_Comm_free(ISR_comm, ierr)
#     endif

      return
   end subroutine ISR_cleanup

   ! ==========================================================================

   subroutine ISR_stop (message)

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      character(len=*), intent(in) :: message
      character(len=75) :: f_str1, f_str2
      integer :: L0, L4
#     if WITHMPI==1
      integer :: ierr
#     endif

      L0 = len_trim(message)
      L4 = L0 + 4
      write(f_str1,"(a5,i2,a1)") "(3x,a", L4, ")"
      write(f_str2,"(a10,i2,a6)") "(3x,'* ',a", L0, ",' *')"

#     if WITHMPI==1
      if (glb_mype == 0) then
#     endif
         write(*,*) ""
         write(*,f_str1) repeat("*",L4)
         write(*,f_str2) trim(message)
         write(*,f_str1) repeat("*",L4)
#     if WITHMPI==1
      end if
#     endif
#     if WITHMPI==1
      call MPI_Finalize(ierr)
#     endif

      stop

      return
   end subroutine ISR_stop

end module initial_state_refinement

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

