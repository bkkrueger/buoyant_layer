
! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

module initial_state_refinement

   use amr_parameters, only : dp

   implicit none

   integer :: ISR_N                 ! Number of elements
   integer :: glb_mype, glb_npes    ! Proc ID and # procs in MPI_COMM_WORLD
   integer :: ISR_mype, ISR_npes    ! Proc ID and # procs in vertical slice
#  if WITHMPI==1
   integer :: ISR_comm              ! Communicator for the vertical slice
#  endif

   contains

   ! ==========================================================================

   subroutine ISR_setup(mype, npes)

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! Declare variables -----------------------------------------------------
      integer, intent(in) :: mype, npes
#     if WITHMPI==1
      integer :: ierr
#     endif

      ! Set up ----------------------------------------------------------------

      ! Size of vector segments
      ISR_N = 18

      ! Defaults for MPI stuff (in case of no MPI being used)
      ISR_npes = 1
      ISR_mype = 0

#     if WITHMPI==1
      ! MPI stuff

      glb_mype = mype
      glb_npes = npes

      ! Partition
      call MPI_Comm_split(MPI_COMM_WORLD, 0, mype, ISR_comm, ierr)
      if ((ierr /= MPI_SUCCESS) .or. (ISR_comm == MPI_COMM_NULL)) then
         call ISR_stop("ERROR: Failure in MPI_Comm_split")
      end if

      ! Get the size of the new communicator
      call MPI_Comm_size(ISR_comm, ISR_npes, ierr)
      if (ISR_npes /= npes) then
         call ISR_stop("ERROR: Unexpected size for ISR_comm")
      end if

      ! Get the rank within the new communicator
      call MPI_Comm_rank(ISR_comm, ISR_mype, ierr)
      if (ISR_mype /= mype) then
         call ISR_stop("ERROR: Unexpected rank in ISR_comm")
      end if

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

