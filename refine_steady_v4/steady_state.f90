module steady

   use amr_parameters

   implicit none

   ! Verbosity for this module: do I print a lot of status messages or keep
   ! quiet?
   logical :: ss_verbose

   ! The steady state
   ! -- Note: the steady state need only be defined where the source/sink terms
   !    are non-zero, hence the existence of ss_nx.  However, when I print the
   !    steady state I want to print all values of x that exist for the
   !    simulation, which includes large regions where the source/sink terms
   !    are zero.  In those regions, we extend by boundary conditions: the
   !    upper boundary is the prescribed inflow and the lower boundary is
   !    zero-gradient.
   real(dp), dimension(:,:), allocatable :: ss0    ! storage for the state
   integer, parameter :: ss_nv = 3                 ! number of variables
   integer :: ss_nx                                ! number of cells
   real(dp), ss_bc(ss_nv)                          ! upper boundary values

   ! The number of guard cells
   ! -- This also tells us information about the stencil.  The stencil is no
   !    larger than i-ss_ng to i+ss_ng.  It may be smaller, but this allows us
   !    to make use of information already in DUMSES to ensure that we are
   !    consistent with certain types of future changes (e.g., implementing a
   !    larger stencil, which necessitates more guard cells, which would
   !    propagate into our method by this construction).
   integer :: ss_ng

#  if WITHMPI==1
   integer :: ss_comm, ss_mype, ss_npes
   integer :: ss_neighlo, ss_neighhi
#  endif

   contains

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine ss_setup ()

      ! Read the config-file parameters _______________________________________
      ! Defaults
      ss_verbose = .true.
      ! Read from file
      open(unit=LUN, file='input', status='old')
      read(LUN, ssc_params)
      close(LUN)

      ! Determine ss_nx _______________________________________________________
      ss_ng = max(1 - iu1, iu2 - nx)
      ihi = iu1
      ilo = iu2
      do idx = iu1, iu2
         if (shape_function(x(idx)) /= 0.0d0) then
            if (idx > ihi) then
               ihi = idx
            end if
            if (idx < ilo) then
               ilo = idx
            end if
         end if
      end do
      ilo = ilo - 2*ss_ng
      ihi = ihi + 2*ss_ng
      ss_nx = ihi - ilo + 1

      ! Allocate ss0
      allocate(ss0(ilo:ihi,ss_nv))

      ! Set up MPI ____________________________________________________________
#     if WITHMPI==1

      ! Limit the number of processors so that slices are big enough to contain
      ! at least one complete stencil
      ss_npes = min(npes, ss_nx / (2*ss_ng+1))

      ! Partition: create a communicator holding only the processors to be used
      ! for the correction of the steady state
      if (mype < ss_npes) then
         mpi_color = 1
      else
         mpi_color = 2
      end if
      call MPI_Comm_split(MPI_COMM_WORLD, mpi_color, mype, ss_comm, ierr)
      if ((ierr /= MPI_SUCCESS) .or. (ss_comm == MPI_COMM_NULL)) then
         call ss_abort("ERROR: Failure in MPI_Comm_split")
      end if

      ! Verify the partition is the correct size
      call MPI_Comm_Size(ss_comm, tmp_int, ierr)
      if (tmp_int /= ss_npes) then
         call ss_abort("ERROR: Unexpected size for ss_comm")
      end if

      ! Get the rank within the new communicator
      call MPI_Comm_rank(ss_comm, ss_mype, ierr)

      ! Define neighbors
      if (ss_mype == 0) then
         ss_neigh_lo = ss_npes - 1
      else
         ss_neighlo = ss_mype - 1
      end if
      if (ss_mype == ss_npes - 1) then
         ss_neighhi = 0
      else
         ss_neighhi = ss_mype + 1
      end if

      ! Adjust verbosity: only the "master" writes
      if (ss_mype /= 0) then
         ss_verbose = .false.
      end if
#     endif

      ! Upper boundary conditions _____________________________________________
      ss_bc(1) = 1.0d0
      ss_bc(2) = -Mup
      ss_bc(3) = 1.0d0/(gamma*(gamma-1.0d0)) + 0.5d0*Mup**2

      ! Build initial guess ___________________________________________________
      call ss_initial_guess()

   end subroutine ss_setup

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine ss_cleanup ()

      ! Save steady state _____________________________________________________
      call ss_write_to_file()

      ! Load initial conditions _______________________________________________
      do idx = iu1, iu2
         if (iu1 < ilo) then
            uin(1,idx,:,:,:) = uin(1,ilo,:,:,:)
         else if (iu2 > ihi) then
            uin(1,idx,:,:,:) = 0.0d0
            uin(1,idx,:,:,1) = ss_bc(1)
            uin(1,idx,:,:,2) = ss_bc(2)
            uin(1,idx,:,:,5) = ss_bc(3)
         else
            uin(1,idx,:,:,:) = 0.0d0
            uin(1,idx,:,:,1) = ss0(idx,1)
            uin(1,idx,:,:,2) = ss0(idx,2)
            uin(1,idx,:,:,5) = ss0(idx,3)
         end if
      end do
      time = 0.0d0
      dt   = 0.0d0

      ! Deallocate ____________________________________________________________
      deallocate(ss0)

      ! Clean up MPI __________________________________________________________
      call MPI_Comm_free(ss_comm, ierr)

   end subroutine ss_cleanup

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine ss_abort ()

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
      if (ss_mype == 0) then
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
   end subroutine ss_abort

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine ss_initial_guess ()

      ! Put in Jerome's method to build the initial guess at the steady state

   end subroutine ss_initial_guess

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine ss_write_to_file ()

      ! Print the steady state to a file

   end subroutine ss_write_to_file

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

end module steady
