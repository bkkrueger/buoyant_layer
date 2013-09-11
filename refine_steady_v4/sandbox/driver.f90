program driver

   use amr_parameters, only : dp
   use ISR_linalg
   use initial_state_refinement, only : ISR_setup, ISR_cleanup
   use initial_state_refinement, only : ISR_N, ISR_mype, ISR_npes
#  if WITHMPI==1
   use initial_state_refinement, only : ISR_comm
#  endif

   implicit none

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   ! ==========================================================================
   ! Declare variables

   ! MPI
   integer :: mype, npes
#  if WITHMPI==1
   integer :: stat(MPI_STATUS_SIZE)
   integer :: ierr
#  endif

   ! Matrices and vectors
   real(dp), allocatable :: A(:,:)
   real(dp), allocatable :: solution(:), b(:), x(:)
   real(dp) :: max_error

   ! Indices
   integer :: proc, idx
   integer :: row, col

   ! ==========================================================================
   ! Set up

#  if WITHMPI==1
   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
   if (mype == 0) then
      write(*,*) "Testing BiCGStab with MPI."
   end if
#  else
   write(*,*) "Testing BiCGStab without MPI."
   mype = 0
   npes = 1
#  endif

   call ISR_setup(mype, npes)
   allocate(A(ISR_N,ISR_N*ISR_npes))
   allocate(b(ISR_N))
   allocate(x(ISR_N))
   allocate(solution(ISR_N))

   ! Create a matrix
   call fill_matrix(A)

   do proc = 0, ISR_npes-1
      do idx = 1, ISR_N

         ! Construct b
         b(:) = 0.0d0
         if (proc == ISR_mype) then
            b(idx) = 1.0d0
         end if

         ! Construct initial guess
         call copy(x, b)
         !x(:) = 0.0d0

         ! Solve Ax=b
         call BiCGStab(A, x, b)

         ! Check the solution
         call matvec(solution, A, x)
         max_error = maxval(abs(solution(:)-b(:)))
         if (ISR_mype == 0) then
            write(*,"(3x,i5,3x,es13.6)") proc*ISR_N+idx, max_error
         end if

         if (max_error >= 1.0d-9) then
            x(:) = 0.0d0
            call BiCGStab(A, x, b)
            call matvec(solution, A, x)
            max_error = maxval(abs(solution(:)-b(:)))
            if (ISR_mype == 0) then
               write(*,"(3x,i5,3x,es13.6)") proc*ISR_N+idx, max_error
            end if
         end if

      end do
   end do

   ! ==========================================================================
   ! Clean up

   deallocate(A,b,x,solution)

   call ISR_cleanup()

   call MPI_Finalize(ierr)

   stop
end program driver
