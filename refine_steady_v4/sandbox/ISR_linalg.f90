! TODO : Fix the prototypes; I've got a terrible mix of func(dst, src1..N) and
! func(src1..N, dst) that makes no sense and is likely to cause problems.

module ISR_linalg

   contains

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine invert_matrix (Ainv, A)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_stop, ISR_mype, ISR_npes
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
      real(dp), intent(  out) :: Ainv(ISR_N, ISR_N*ISR_npes)
      real(dp), intent(in   ) :: A   (ISR_N, ISR_N*ISR_npes)

      ! Locals ----------------------------------------------------------------
      integer :: proc, row, col, idx
      real(dp) :: tmp, max_err, mean_err
      real(dp) :: Ic(ISR_N)
#     if WITHMPI==1
      integer :: ierr
#     endif
      real(dp) :: error
      real(dp), dimension(ISR_N) :: temp_vec1, temp_vec2

      ! =======================================================================
      ! Invert matrix

      ! NOTE: Gaussian elimination is simpler in concept, and therefore good
      ! for initial testing to make sure things work as expected.  However, it
      ! is likely to be inefficient.  Additionally, it is difficult to
      ! parallelize when the vectors and matrices are distributed across
      ! multiple processors.  Thus is is better to use BiCGStab to solve an
      ! Ax=b type of problem for each column of the inverse.  I will have to
      ! look into whether or not there is some generalization of BiCGStab to
      ! solve AX=B (X and B matrices instead of vectors).

!      ! Gaussian elimination
!      call get_full_matrix(A, A_full)
!      call row_reduce_invert(Ainv_full, A_full)
!      call slice_matrix(Ainv_full, Ainv)

      ! BiCGStab
      do proc = 0, ISR_npes-1
         do col = 1, ISR_N

            if (glb_mype == 0) then
               write(*,"(3x,'invert matrix: column ',i5)") proc*ISR_N+col
            end if

            ! Construct b (column of identity matrix)
            Ic(:) = 0.0d0
            if (proc == ISR_mype) then
               Ic(col) = 1.0d0
            end if

            ! Construct guess for x
            Ainv(:,proc*ISR_N+col) = Ic(:)

            ! Solve A x = I_c for column of A^{-1}
            call BiCGStab(A, Ainv(:,proc*ISR_N+col), Ic)

            ! TODO
            ! Try a second version if necessary
            call matvec(temp_vec1, A, Ainv(:,proc*ISR_N+col))
            ! We know that b = Ic, so magnitude of b = 1
            call axpy(temp_vec2, -1.0d0, temp_vec1, Ic)
            error = maxval(abs(temp_vec2))
            if (error >= 1.0d-9) then
               Ainv(:,proc*ISR_N+col) = 0.0d0
               call BiCGStab(A, Ainv(:,proc*ISR_N+col), Ic)
               call matvec(temp_vec1, A, Ainv(:,proc*ISR_N+col))
               call axpy(temp_vec2, -1.0d0, temp_vec1, Ic)
               error = maxval(abs(temp_vec2))
               if (error >= 1.0d-9) then
                  call ISR_stop("Failure in BiCGStab.")
               end if
            end if

         end do
      end do

      ! =======================================================================
      ! Check inverse

      max_err  = 0.0d0
      mean_err = 0.0d0
      do col = 1, ISR_N*ISR_npes

         ! Compute column [:,col] of (A*Ainv)
         call matvec(Ic(:), A(:,:), Ainv(:,col))

         ! Compute maximum error in column
         do row = 1, ISR_N

            ! Compute error (absolute difference between (A*Ainv) and I)
            if (ISR_mype*ISR_N+row == col) then
               tmp = abs(1.0d0 - Ic(row))
            else
               tmp = abs(Ic(row))
            end if

            ! Update errors within (A*Ainv)
            max_err  = max(max_err, tmp)
            mean_err = mean_err + tmp
         end do
#        if WITHMPI==1
         call MPI_Allreduce(max_err, tmp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                            ISR_comm, ierr)
         max_err = tmp
#        endif

      end do
#     if WITHMPI==1
      call MPI_Allreduce(mean_err, tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         ISR_comm, ierr)
      mean_err = tmp / (ISR_N*ISR_npes)**2
#     endif

      if (glb_mype == 0) then
         write(*,"(3x,'inverse errors:',2(x,es13.6))") max_err, mean_err
      end if

      return
   end subroutine invert_matrix

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine solve_Axb (A, x, b)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes
#     if WITHMPI==1
      use initial_state_refinement, only : glb_mype
#     endif

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(in   ) :: A(ISR_N, ISR_N*ISR_npes)
      real(dp), intent(  out) :: x(ISR_N)
      real(dp), intent(in   ) :: b(ISR_N)

      ! Locals ----------------------------------------------------------------

      ! =======================================================================
      ! Solve vector equation

      ! NOTE: Gaussian elimination is simpler in concept, and therefore good
      ! for initial testing to make sure things work as expected.  However, it
      ! is likely to be inefficient.  Additionally, it is difficult to
      ! parallelize when the vectors and matrices are distributed across
      ! multiple processors.  Thus is is better to use BiCGStab.

!      call row_reduce_Axb(A, x, b)
      x(:) = 0.0d0
      call BiCGStab(A, x, b)

      return
   end subroutine solve_Axb

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine row_reduce_invert (Ainv, A)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N
#     if WITHMPI==1
      use initial_state_refinement, only : glb_mype
#     endif

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(  out) :: Ainv(ISR_N, ISR_N)
      real(dp), intent(in   ) :: A   (ISR_N, ISR_N)

      ! Locals ----------------------------------------------------------------
      integer :: row, col
      real(dp) :: tmp, M(ISR_N,2*ISR_N)

      ! =======================================================================
      ! Ugly method: Gaussian elimination

      ! Construct matrix [ A | I ] --------------------------------------------
      do row = 1, ISR_N
         do col = 1, ISR_N
            M(row,col) = A(row,col)
            if (row == col) then
               M(row,col+ISR_N) = 1.0d0
            else
               M(row,col+ISR_N) = 0.0d0
            end if
         end do
      end do

      ! Row-reduce A to I to end with [ I | Ainv ] ----------------------------

      ! Row-reduce to eliminate all sub-diagonal entries
      do col = 1, ISR_N

         ! Scale row such that diagonal entry is unity
         tmp = M(col,col)
         if (tmp == 0.0d0) then
#           if WITHMPI==1
            if (glb_mype == 0) then
#           endif
               write(*,*) "ZERO ON DIAGONAL AT", col
#           if WITHMPI==1
            end if
#           endif
         end if
         M(col,:) = M(col,:) / tmp

         ! Use row [col] to eliminate M[row,col]
         do row = col+1, ISR_N
            tmp = M(row,col)
            if (tmp /= 0.0d0) then
               M(row,:) = M(row,:) - tmp*M(col,:)
            end if
         end do

      end do

      ! Row-reduce to eliminate all super-diagonal entries
      do col = ISR_N, 1, -1

         ! Scale row such that diagonal entry is unity
         tmp = M(col,col)
         if (tmp == 0.0d0) then
#           if WITHMPI==1
            if (glb_mype == 0) then
#           endif
               write(*,*) "ZERO ON DIAGONAL AT", col
#           if WITHMPI==1
            end if
#           endif
         end if
         M(col,:) = M(col,:) / tmp

         ! Use row [col] to eliminate M[row,col]
         do row = col-1, 1, -1
            tmp = M(row,col)
            if (tmp /= 0.0d0) then
               M(row,:) = M(row,:) - tmp*M(col,:)
            end if
         end do

      end do

      ! Extract Ainv ----------------------------------------------------------
      do row = 1, ISR_N
         do col = 1, ISR_N
            Ainv(row,col) = M(row,ISR_N+col)
         end do
      end do

      return
   end subroutine row_reduce_invert

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

!   subroutine row_reduce_Axb (A, x, b)
!
!      use amr_parameters, only : dp
!      use initial_state_refinement, only : ISR_N
!#     if WITHMPI==1
!      use initial_state_refinement, only : glb_mype
!#     endif
!
!      implicit none
!
!      ! =======================================================================
!      ! Declare variables
!
!      ! Arguments -------------------------------------------------------------
!      real(dp), intent(in   ) :: A(ISR_N, ISR_N)
!      real(dp), intent(  out) :: x(ISR_N)
!      real(dp), intent(in   ) :: b(ISR_N)
!
!      ! Locals ----------------------------------------------------------------
!      integer :: row, col
!      real(dp) :: tmp, M(ISR_N,ISR_N+1)
!
!      ! =======================================================================
!      ! Ugly method: Gaussian elimination
!
!      ! Construct matrix [ A | b ] --------------------------------------------
!      do row = 1, ISR_N
!         do col = 1, ISR_N
!            M(row,col) = A(row,col)
!         end do
!         M(row,ISR_N+1) = b(row)
!      end do
!
!      ! Row-reduce A to I to end with [ I | x ] -------------------------------
!
!      ! Row-reduce to eliminate all sub-diagonal entries
!      do col = 1, ISR_N
!
!         ! Scale row such that diagonal entry is unity
!         tmp = M(col,col)
!         if (tmp == 0.0d0) then
!#           if WITHMPI==1
!            if (glb_mype == 0) then
!#           endif
!               write(*,*) "ZERO ON DIAGONAL AT", col
!#           if WITHMPI==1
!            end if
!#           endif
!         end if
!         M(col,:) = M(col,:) / tmp
!
!         ! Use row [col] to eliminate M[row,col]
!         do row = col+1, ISR_N
!            tmp = M(row,col)
!            if (tmp /= 0.0d0) then
!               M(row,:) = M(row,:) - tmp*M(col,:)
!            end if
!         end do
!
!      end do
!
!      ! Row-reduce to eliminate all super-diagonal entries
!      do col = ISR_N, 1, -1
!
!         ! Scale row such that diagonal entry is unity
!         tmp = M(col,col)
!         if (tmp == 0.0d0) then
!#           if WITHMPI==1
!            if (glb_mype == 0) then
!#           endif
!               write(*,*) "ZERO ON DIAGONAL AT", col
!#           if WITHMPI==1
!            end if
!#           endif
!         end if
!         M(col,:) = M(col,:) / tmp
!
!         ! Use row [col] to eliminate M[row,col]
!         do row = col-1, 1, -1
!            tmp = M(row,col)
!            if (tmp /= 0.0d0) then
!               M(row,:) = M(row,:) - tmp*M(col,:)
!            end if
!         end do
!
!      end do
!
!      ! Extract x -------------------------------------------------------------
!      do row = 1, ISR_N
!         x(row) = M(row,ISR_N+1)
!      end do
!
!      return
!   end subroutine row_reduce_Axb

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine BiCGStab(A, x, b)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes, glb_mype

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(in   ) :: A(ISR_N, ISR_N*ISR_npes)
      real(dp), intent(  out) :: x(ISR_N)
      real(dp), intent(in   ) :: b(ISR_N)

      ! Locals ----------------------------------------------------------------

      ! Temporaries for the algorithm
      real(dp), dimension(ISR_N) :: r, r_tilde, v, t, q, p, p_hat, s, s_hat
      real(dp) :: rho_curr, rho_prev
      real(dp) :: alpha, beta, omega, mag_sq_b, d1, d2
real(dp) :: tmp1, tmp2
integer :: idx

      ! Iteration/convergence parameters
      integer :: max_iter, iter
      real(dp), parameter :: TOLERANCE = 1.0d-9
      real(dp) :: error

      ! =======================================================================
      ! Bi-conjugate gradient (stabilized) method

      rho_curr = 0.0d0
      rho_prev = 0.0d0

      max_iter = (ISR_N*ISR_npes)**2

      ! Setup
      call inner_product(mag_sq_b, b, b)
      if (mag_sq_b == 0.0d0) then
         call copy(x,b)
      else
         ! Construct initial residual: r = b - Ax
         call matvec(q, A, x)
         call axpy(r, -1.0d0, q, b)

         ! Choose r_tilde
         call copy(r_tilde, r)

         ! Loop
         do iter = 1, max_iter

            ! Update and check rho
            rho_prev = rho_curr
            call inner_product(rho_curr, r_tilde, r)
            if (rho_curr == 0.0d0) then
               call inner_product(d1, r, r)
               error = sqrt(d1 / mag_sq_b)
               exit
            end if

            ! Compute p
            if (iter == 1) then
               call copy(p, r)
            else
               beta = (rho_curr / rho_prev) * (alpha / omega)
               call axpy(q, -1.0d0*omega, v, p)
               call axpy(p,        beta , q, r)
            end if

            ! Precondition p --- skip
            call copy(p_hat, p)

            ! Compute v
            call matvec(v, A, p_hat)

            ! Compute alpha
            call inner_product(d1, r_tilde, v)
            alpha = rho_curr / d1
!do idx = 1, ISR_N
!if ((r_tilde(idx) /= 0.0d0) .and. (v(idx) /= 0.0d0)) then
!write(*,*) idx, r_tilde(idx), v(idx)
!end if
!end do
!call inner_product(tmp1, r_tilde, r_tilde)
!call inner_product(tmp2, v, v)
!if (glb_mype == 0) then
!write(*,"('ALPHA = ',10(3x,es13.6))") alpha, d1, tmp1, tmp2
!end if

            ! Compute s
            call axpy(s, -1.0d0*alpha, v, r)

            ! Check norm of s
            call inner_product(d1, s, s)
            error = sqrt(d1 / mag_sq_b)
            if (error < TOLERANCE) then
               ! Borrow s_hat as an intermediate
               call axpy(s_hat, alpha, p_hat, x)
               call copy(x, s_hat)
               exit
            end if

            ! Precondition s --- skip
            call copy(s_hat, s)

            ! Compute t
            call matvec(t, A, s_hat)

            ! Compute omega
            call inner_product(d1, t, s)
            call inner_product(d2, t, t)
            omega = d1 / d2

            ! Update x
            call axpy(q, alpha, p_hat, x)
            call axpy(x, omega, s_hat, q)

            ! Update r
            call axpy(r, -1.0d0*omega, t, s)

            ! Check convergence
            call inner_product(d1, r, r)
            error = sqrt(d1 / mag_sq_b)
            if (error <= TOLERANCE) then
               exit
            end if

            ! Check omega
            if (omega == 0.0d0) then
               exit
            end if

         end do

         if (error > TOLERANCE) then
            if (glb_mype == 0) then
               write(*,"('NON-CONVERGENCE IN BICGSTAB: ',es13.6)") error
            end if
         end if

      end if

      return
   end subroutine BiCGStab

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine inner_product (scl, vec1, vec2)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N
#     if WITHMPI==1
      use initial_state_refinement, only : ISR_comm
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp),                   intent(  out) :: scl
      real(dp), dimension(ISR_N), intent(in   ) :: vec1, vec2

      ! Locals ----------------------------------------------------------------
      real(dp) :: tmp
#     if WITHMPI==1
      integer :: ierr
#     endif

      ! =======================================================================
      ! Compute the inner product of vectors vec1 and vec2

      tmp = sum(vec1(:)*vec2(:))

#     if WITHMPI==1
      call MPI_Allreduce(tmp, scl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         ISR_comm, ierr)
#     else
      scl = tmp
#     endif

      return
   end subroutine inner_product

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine outer_product (mat, vec1, vec2)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes
#     if WITHMPI==1
      use initial_state_refinement, only : ISR_comm
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), dimension(ISR_N,ISR_N*ISR_npes), intent(  out) :: mat
      real(dp), dimension(ISR_N)               , intent(in   ) :: vec1, vec2

      ! Locals ----------------------------------------------------------------
      integer :: i, j
      real(dp), dimension(ISR_N*ISR_npes) :: vec2_full
#     if WITHMPI==1
      integer :: ierr
#     endif

      ! =======================================================================
      ! Compute the outer product of vectors vec1 and vec2

#     if WITHMPI==1
      call MPI_Allgather(vec2     , ISR_N, MPI_DOUBLE_PRECISION, &
                         vec2_full, ISR_N, MPI_DOUBLE_PRECISION, &
                         ISR_comm , ierr)
#     else
      call copy(vec2_full, vec2)
#     endif

      do i = 1, ISR_N
         do j = 1, ISR_N*ISR_npes
            mat(i,j) = vec1(i) * vec2_full(j)
         end do
      end do

      return
   end subroutine outer_product

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine copy (dst, src)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(  out) :: dst(ISR_N)
      real(dp), intent(in   ) :: src(ISR_N)

      ! Locals ----------------------------------------------------------------

      ! =======================================================================
      ! Copy vector src to vector dst

      dst(:) = src(:)

      return
   end subroutine copy

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine axpy (dst, a, x, y)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), dimension(ISR_N), intent(  out) :: dst
      real(dp),                   intent(in   ) :: a
      real(dp), dimension(ISR_N), intent(in   ) :: x, y

      ! Locals ----------------------------------------------------------------

      ! =======================================================================
      ! Compute dst = a*x + y

      dst(:) = a*x(:) + y(:)

      return
   end subroutine axpy

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine mat_transpose (dst, src)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes
#     if WITHMPI==1
      use initial_state_refinement, only : ISR_comm, ISR_mype
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(  out) :: dst(ISR_N, ISR_N*ISR_npes)
      real(dp), intent(in   ) :: src(ISR_N, ISR_N*ISR_npes)

      ! Locals ----------------------------------------------------------------
      integer :: row, col
#     if WITHMPI==1
      integer :: proc
      integer :: sendtag, recvtag
      integer :: requests(2*ISR_npes), statuses(MPI_STATUS_SIZE, 2*ISR_npes)
      real(dp) :: B(ISR_npes, ISR_N, ISR_N)
      real(dp), dimension(ISR_N**2, ISR_npes) :: sendbufs, recvbufs
      integer :: ierr
#     endif

      ! =======================================================================
      ! Compute the transpose of the matrix

#     if WITHMPI==1
      ! parallel version ------------------------------------------------------

      ! Split into blocks
      do row = 1, ISR_N
         do col = 1, ISR_N
            do proc = 0, ISR_npes-1
               B(proc+1,row,col) = src(row,proc*ISR_N+col)
            end do
         end do
      end do

      ! Pack the send buffers
      do row = 1, ISR_N
         do col = 1, ISR_N
            do proc = 1, ISR_npes-1
               sendbufs((row-1)*ISR_N+col,proc+1) = B(proc+1,row,col)
            end do
         end do
      end do

      ! Send and receive
      do proc = 0, ISR_npes-1
         sendtag = ISR_mype * ISR_npes + proc
         call MPI_Isend(sendbufs(:,proc+1), ISR_N**2, MPI_INTEGER, proc, &
                        sendtag, ISR_comm, requests(proc+1         ), ierr)
         recvtag = proc     * ISR_npes + ISR_mype
         call MPI_Irecv(recvbufs(:,proc+1), ISR_N**2, MPI_INTEGER, proc, &
                        recvtag, ISR_comm, requests(proc+1+ISR_npes), ierr)
      end do
      call MPI_Waitall(2*ISR_npes, requests, statuses, ierr)

      ! Unpack the receive buffers
      do row = 1, ISR_N
         do col = 1, ISR_N
            do proc = 1, ISR_npes-1
               ! Don't overwrite the block that was not transferred
               if (proc /= ISR_mype) then
                  B(proc+1,row,col) = recvbufs((row-1)*ISR_N+col,proc+1)
               end if
            end do
         end do
      end do

      ! Transpose block and copy into matrix
      do row = 1, ISR_N
         do col = 1, ISR_N
            do proc = 0, ISR_npes-1
               dst(row,proc*ISR_N+col) = B(proc+1,col,row)
            end do
         end do
      end do

#     else
      ! serial version --------------------------------------------------------
      do row = 1, ISR_N
         do col = 1, ISR_N
            dst(row,col) = src(col,row)
         end do
      end do
#     endif

      return
   end subroutine mat_transpose

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine matvec (dst, mat, src)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes
#     if WITHMPI==1
      use initial_state_refinement, only : ISR_comm, ISR_mype
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(  out) :: dst(ISR_N)
      real(dp), intent(in   ) :: mat(ISR_N, ISR_N*ISR_npes)
      real(dp), intent(in   ) :: src(ISR_N)

      ! Locals ----------------------------------------------------------------
      integer :: i
#     if WITHMPI==1
      real(dp) :: src_full(ISR_N*ISR_npes)
      integer :: ierr
#     endif

      ! =======================================================================
      ! Compute the matrix-vector product: vec1 = mat vec2

#     if WITHMPI==1
      call MPI_Allgather(src     , ISR_N, MPI_DOUBLE_PRECISION, &
                         src_full, ISR_N, MPI_DOUBLE_PRECISION, &
                         ISR_comm , ierr)
#     else
      call copy(src_full, src)
#     endif

      do i = 1, ISR_N
         dst(i) = sum(mat(i,:)*src_full(:))
      end do

      return
   end subroutine matvec

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine Broyden_update_inverse (dss, dff, Jinv)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), dimension(ISR_N)               , intent(in   ) :: dss, dff
      real(dp), dimension(ISR_N,ISR_N*ISR_npes), intent(inout) :: Jinv

      ! Locals ----------------------------------------------------------------
      real(dp), dimension(ISR_N,ISR_N*ISR_npes) :: M
      real(dp), dimension(ISR_N) :: uu, vv, ww, numer
      real(dp) :: denom
      integer :: i

      ! =======================================================================
      ! Use Broyden's method to update the inverse Jacobian

      ! Choose vv
      ! --- Broyden's method one
      call mat_transpose(M, Jinv)
      call matvec(vv, M, dss)
      ! --- Broyden's method two
      !call copy(vv, dff)

      ! Compute uu
      call matvec(uu, Jinv, dff)
      call axpy(numer, -1.0d0, uu, dss)
      call inner_product(denom, vv, dff)
      denom = 1.0d0 / denom
      uu(:) = numer(:) * denom

      ! Compute rank-one update (outer product: uu x vv)
      call outer_product(M, uu, vv)

      ! Update inverse Jacobian
      Jinv(:,:) = Jinv(:,:) + M(:,:)

      return
   end subroutine Broyden_update_inverse

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine get_full_matrix (M, M_full)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes
#     if WITHMPI==1
      use initial_state_refinement, only : ISR_comm
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(in   ) :: M     (ISR_N         , ISR_N*ISR_npes)
      real(dp), intent(inout) :: M_full(ISR_N*ISR_npes, ISR_N*ISR_npes)

      ! Locals ----------------------------------------------------------------
#     if WITHMPI==1
      real(dp), dimension(ISR_N*ISR_N*ISR_npes) :: sendbuf
      real(dp), dimension((ISR_N*ISR_npes)**2)  :: recvbuf
      integer :: row, col
      integer :: ierr
#     endif

      ! =======================================================================
      ! Merge the distributed matrix so all processes have a full copy

#     if WITHMPI==1
      ! Pack
      do row = 1, ISR_N
         do col = 1, ISR_N*ISR_npes
            sendbuf((row-1)*ISR_N*ISR_npes+col) = M(row,col)
         end do
      end do

      ! Send
      call MPI_Allgather(sendbuf, ISR_N*ISR_N*ISR_npes, MPI_DOUBLE_PRECISION, &
                         recvbuf, ISR_N*ISR_N*ISR_npes, MPI_DOUBLE_PRECISION, &
                         ISR_comm, ierr)

      ! Unpack
      do row = 1, ISR_N*ISR_npes
         do col = 1, ISR_N*ISR_npes
            M_full(row,col) = recvbuf((row-1)*ISR_N*ISR_npes+col)
         end do
      end do

#     else
      M_full(:,:) = M(:,:)
#     endif

      return
   end subroutine get_full_matrix

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine get_full_vector (v, v_full)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes
#     if WITHMPI==1
      use initial_state_refinement, only : ISR_comm
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(in   ) :: v     (ISR_N         )
      real(dp), intent(inout) :: v_full(ISR_N*ISR_npes)

      ! Locals ----------------------------------------------------------------
#     if WITHMPI==1
      integer :: ierr
#     endif

      ! =======================================================================
      ! Merge the distributed matrix so all processes have a full copy

#     if WITHMPI==1
      ! Send
      call MPI_Allgather(v     , ISR_N         , MPI_DOUBLE_PRECISION, &
                         v_full, ISR_N*ISR_npes, MPI_DOUBLE_PRECISION, &
                         ISR_comm, ierr)
#     else
      v_full(:,:) = v(:,:)
#     endif

      return
   end subroutine get_full_vector

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine slice_matrix (M_full, M)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes
#     if WITHMPI==1
      use initial_state_refinement, only : ISR_comm, ISR_mype
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(in   ) :: M_full(ISR_N*ISR_npes, ISR_N*ISR_npes)
      real(dp), intent(inout) :: M     (ISR_N         , ISR_N*ISR_npes)

      ! Locals ----------------------------------------------------------------
#     if WITHMPI==1
      integer :: offset
#     endif

      ! =======================================================================
      ! Merge the distributed matrix so all processes have a full copy

#     if WITHMPI==1
      ! Calculate offset
      offset = ISR_mype * ISR_N

      ! Slice
      M(:,:) = M_full(offset+1:offset+ISR_N,:)
#     else
      M(:,:) = M_full(:,:)
#     endif

      return
   end subroutine slice_matrix

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine slice_vector (v_full, v)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_npes
#     if WITHMPI==1
      use initial_state_refinement, only : ISR_comm, ISR_mype
#     endif

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(in   ) :: v_full(ISR_N*ISR_npes)
      real(dp), intent(inout) :: v     (ISR_N         )

      ! Locals ----------------------------------------------------------------
#     if WITHMPI==1
      integer :: offset
#     endif

      ! =======================================================================
      ! Merge the distributed matrix so all processes have a full copy

#     if WITHMPI==1
      ! Calculate offset
      offset = ISR_mype * ISR_N

      ! Slice
      v(:) = v_full(offset+1:offset+ISR_N)
#     else
      v(:,:) = v_full(:,:)
#     endif

      return
   end subroutine slice_vector

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

end module ISR_linalg

