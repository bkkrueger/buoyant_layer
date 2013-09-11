module ISR_linalg

   contains

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine computes the inverse of matrix A

   subroutine invert_matrix (Ainv, A)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_mype, ISR_npes
#     if WITHMPI==1
      use initial_state_refinement, only : ISR_comm, glb_mype
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(inout) :: Ainv(ISR_N, ISR_N*ISR_npes)
      real(dp), intent(in   ) :: A   (ISR_N, ISR_N*ISR_npes)

      ! Locals ----------------------------------------------------------------
      integer :: proc, idx, col
      real(dp), dimension(ISR_N) :: ident_col, tmp1, tmp2
      real(dp) :: error, max_error, mean_error
      real(dp) :: tmp_scalar
#     if WITHMPI==1
      integer :: ierr
#     endif

      ! =======================================================================
      ! Invert matrix

      ! Loop over all columns
      do proc = 0, ISR_npes-1
         do idx = 1, ISR_N

            ! Compute column number
            col = proc*ISR_N + idx

            ! Print a note
#           if WITHMPI==1
            if (glb_mype == 0) then
#           endif
               write(*,"(3x,'invert matrix (column ',i5,')')") col
#           if WITHMPI==1
            endif
#           endif

            ! Construct the corresponding column of the identity matrix
            ident_col(:) = 0.0d0
            if (proc == ISR_mype) then
               ident_col(idx) = 1.0d0
            end if

            ! Construct a guess for the inverse column
            tmp1(:) = ident_col(:)

            ! Solve A x = I(:,col)
            call BiCGStab(tmp1, A, ident_col)

            ! Copy into inverse
            Ainv(:,col) = tmp1(:)

         end do
      end do

      ! =======================================================================
      ! Check inversion

      max_err = 0.0d0
      mean_err = 0.0d0

      ! Loop over all the columns
      do proc = 0, ISR_npes-1
         do idx = 1, ISR_N
            col = proc*ISR_N + col

            ! Compute A * Ainv(:,col)
            tmp1(:) = Ainv(:,col)
            call matvec(tmp2, A, tmp1)

            ! Compare it with Icol
            ident_col(:) = 0.0d0
            if (proc == ISR_mype) then
               ident_col(idx) = 1.0d0
            end if
            call axpy(tmp1, -1.0d0, tmp2, ident_col)

            ! Accumulate max and mean errors
            do row = 1, ISR_N
               error = abs(tmp1(row))
               mean_error = mean_error + error
               max_error  = max(max_error, error)
            end do
         end do
      end do
#     if WITHMPI==1
      tmp_scalar = mean_error
      call MPI_Allreduce(tmp_scalar, mean_error, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, ISR_comm, ierr)
      tmp_scalar = max_error
      call MPI_Allreduce(tmp_scalar, max_error , 1, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, ISR_comm, ierr)
      max_error = max_error / (ISR_N*ISR_npes)**2
#     endif

      ! Print error
#     if WITHMPI==1
      if (glb_mype == 0) then
#     endif
         write(*,"(3x,'inversion error: ',2(x,es13.6))") max_error, mean_error
#     if WITHMPI==1
      endif
#     endif

      return
   end subroutine invert_matrix

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine solves the equation Ax=b where A is a known matrix, b is a
   ! known vector, and x is an unknown vector

   subroutine solve_Axb (x, A, b)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(inout) :: x(ISR_N)
      real(dp), intent(in   ) :: A(ISR_N, ISR_N*ISR_npes)
      real(dp), intent(in   ) :: b(ISR_N)

      ! =======================================================================
      ! Solve Ax = b

      ! Pick a guess for x
      x(:) = b(:)

      ! Call BiCGStab
      call BiCGStab(x, A, b)

      return
   end subroutine solve_Axb

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine uses the Bi-Conjugate Gradient (Stabilized) algorithm to
   ! solve an equation of the form Ax=b

   subroutine BiCGStab (x, A, b)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes, glb_mype

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(inout) :: x(ISR_N)
      real(dp), intent(in   ) :: A(ISR_N, ISR_N*ISR_npes)
      real(dp), intent(in   ) :: b(ISR_N)

      ! Locals ----------------------------------------------------------------

      ! =======================================================================
      ! Bi-Conjugate Gradient (Stabilized)

      return
   end subroutine BiCGStab

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine computes the inner product of two vectors and returns the
   ! scalar result

   subroutine inner_product (scalar, vector1, vector2)

      use amr_parameters, only : dp

      implicit none

      return
   end subroutine inner_product

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine computes the outer product of two vectors and returns the
   ! matrix result

   subroutine outer_product (matrix, vector1, vector2)

      use amr_parameters, only : dp

      implicit none

      return
   end subroutine outer_product

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine copies vector v_in to vector v_out

   subroutine copy (v_out, v_in)

      use amr_parameters, only : dp

      implicit none

      return
   end subroutine copy

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine computes a*x+y, where x and y are vectors and a is a scalar

   subroutine axpy (v_out, a, x, y)

      use amr_parameters, only : dp

      implicit none

      return
   end subroutine axpy

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine computes the transpose of a matrix

   subroutine matrix_transpose (A_dst, A_src)

      use amr_parameters, only : dp

      implicit none

      return
   end subroutine matrix_transpose

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine computes the product of a matrix with a vector: M*v

   subroutine matvec (v_out, A, v_in)

      use amr_parameters, only : dp

      implicit none

      return
   end subroutine matvec

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine print_vector (vec)
   end subroutine print_vector

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine print_matrix (mat)
   end subroutine print_matrix

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------
   ! This routine computes an update to the inverse Jacobian for the Broyden
   ! method

   subroutine Broyden_update_inverse (dss, dff, Jinv)

      use amr_parameters, only : dp

      implicit none

      return
   end subroutine Broyden_update_inverse

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

end module ISR_linalg
