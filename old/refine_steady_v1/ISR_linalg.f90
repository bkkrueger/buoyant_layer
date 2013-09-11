module ISR_linalg

   contains

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine invert_matrix (Ainv, A)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, ISR_stop
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
      integer :: row, col, idx
      real(dp) :: tmp, max_err
      real(dp) :: Ic(ISR_N)

      ! =======================================================================
      ! Invert matrix

!      ! Gaussian elimination
!      call row_reduce_invert(Ainv, A)

      ! BiCGStab
      do col = 1, ISR_N
         Ic( : ) = 0.0d0
         Ic(col) = 1.0d0
         call BiCGStab(A, Ainv(:,col), Ic)
      end do

      ! =======================================================================
      ! Check inverse

      max_err = 0.0d0
      do row = 1, ISR_N
         do col = 1, ISR_N

            ! Compute (A*Ainv)[r][c]
            tmp = 0.0d0
            do idx = 1, ISR_N
               tmp = tmp + A(row,idx)*Ainv(idx,col)
            end do

            ! Compute error (absolute difference between (A*Ainv) and I)
            if (row == col) then
               tmp = abs(1.0d0 - tmp)
            else
               tmp = abs(tmp)
            end if

            ! Update maximum error within (A*Ainv)
            max_err = max(max_err, tmp)
         end do
      end do

      if (glb_mype == 0) then
         write(*,"('inverse error:',x,es13.6)") max_err
      end if

      return
   end subroutine invert_matrix

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine solve_Axb (A, x, b)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N
#     if WITHMPI==1
      use initial_state_refinement, only : glb_mype
#     endif

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(in   ) :: A(ISR_N, ISR_N)
      real(dp), intent(  out) :: x(ISR_N)
      real(dp), intent(in   ) :: b(ISR_N)

      ! Locals ----------------------------------------------------------------
      real(dp), dimension(ISR_N) :: x_reduce, b_reduce, e_reduce
      real(dp), dimension(ISR_N) :: x_krylov, b_krylov, e_krylov
      real(dp) :: error_reduce, error_krylov
      integer  :: i

      ! =======================================================================
      ! Solve vector equation

!      call row_reduce_Axb(A, x, b)
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

   subroutine row_reduce_Axb (A, x, b)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N
#     if WITHMPI==1
      use initial_state_refinement, only : glb_mype
#     endif

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(in   ) :: A(ISR_N, ISR_N)
      real(dp), intent(  out) :: x(ISR_N)
      real(dp), intent(in   ) :: b(ISR_N)

      ! Locals ----------------------------------------------------------------
      integer :: row, col
      real(dp) :: tmp, M(ISR_N,ISR_N+1)

      ! =======================================================================
      ! Ugly method: Gaussian elimination

      ! Construct matrix [ A | b ] --------------------------------------------
      do row = 1, ISR_N
         do col = 1, ISR_N
            M(row,col) = A(row,col)
         end do
         M(row,ISR_N+1) = b(row)
      end do

      ! Row-reduce A to I to end with [ I | x ] -------------------------------

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

      ! Extract x -------------------------------------------------------------
      do row = 1, ISR_N
         x(row) = M(row,ISR_N+1)
      end do

      return
   end subroutine row_reduce_Axb

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine BiCGStab(A, x, b)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N, glb_mype

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(in   ) :: A(ISR_N, ISR_N)
      real(dp), intent(  out) :: x(ISR_N)
      real(dp), intent(in   ) :: b(ISR_N)

      ! Locals ----------------------------------------------------------------

      ! Temporaries for the algorithm
      real(dp), dimension(ISR_N) :: r, r_tilde, v, t, q, p, p_hat, s, s_hat
      real(dp) :: rho_curr, rho_prev
      real(dp) :: alpha, beta, omega, mag_sq_b, d1, d2

      ! Iteration/convergence parameters
      integer :: max_iter, iter
      real(dp), parameter :: TOLERANCE = 1.0d-9
      real(dp) :: error

      ! =======================================================================
      ! Bi-conjugate gradient (stabilized) method

      rho_curr = 0.0d0
      rho_prev = 0.0d0

      max_iter = ISR_N

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
               error = sqrt(d1)
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

            ! Compute s
            call axpy(s, -1.0d0*alpha, v, r)

            ! Check norm of s
            call inner_product(d1, s, s)
            error = sqrt(d1 / mag_sq_b)
            if (error < TOLERANCE) then
               call axpy(x, alpha, p_hat, x)
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

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp),                   intent(  out) :: scl
      real(dp), dimension(ISR_N), intent(in   ) :: vec1, vec2

      ! Locals ----------------------------------------------------------------

      ! =======================================================================
      ! Compute the inner product of vectors vec1 and vec2

      scl = sum(vec1(:)*vec2(:))

      return
   end subroutine inner_product

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine outer_product (mat, vec1, vec2)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), dimension(ISR_N,ISR_N), intent(  out) :: mat
      real(dp), dimension(ISR_N)      , intent(in   ) :: vec1, vec2

      ! Locals ----------------------------------------------------------------
      integer :: i, j

      ! =======================================================================
      ! Compute the outer product of vectors vec1 and vec2

      do i = 1, ISR_N
         do j = 1, ISR_N
            mat(i,j) = vec1(i) * vec2(j)
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
      use initial_state_refinement, only : ISR_N

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(  out) :: dst(ISR_N, ISR_N)
      real(dp), intent(in   ) :: src(ISR_N, ISR_N)

      ! Locals ----------------------------------------------------------------
      integer :: i, j

      ! =======================================================================
      ! Compute the transpose of the matrix

      do i = 1, ISR_N
         do j = 1, ISR_N
            dst(i,j) = src(j,i)
         end do
      end do

      return
   end subroutine mat_transpose

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine matvec (dst, mat, src)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N
      use variables

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), intent(  out) :: dst(ISR_N)
      real(dp), intent(in   ) :: mat(ISR_N, ISR_N)
      real(dp), intent(in   ) :: src(ISR_N)

      ! Locals ----------------------------------------------------------------
      integer :: i, j
      real(dp) :: tmp

      ! =======================================================================
      ! Compute the matrix-vector product: vec1 = mat vec2

      do i = 1, ISR_N
         tmp = 0.0d0
         do j = 1, ISR_N
            tmp = tmp + mat(i,j) * src(j)
         end do
         dst(i) = tmp
      end do

      return
   end subroutine matvec

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine Broyden_update_inverse (dss, dff, Jinv)

      use amr_parameters, only : dp
      use initial_state_refinement, only : ISR_N

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      real(dp), dimension(ISR_N)      , intent(in   ) :: dss, dff
      real(dp), dimension(ISR_N,ISR_N), intent(inout) :: Jinv

      ! Locals ----------------------------------------------------------------
      real(dp), dimension(ISR_N,ISR_N) :: M
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

end module ISR_linalg

