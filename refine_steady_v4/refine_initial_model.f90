! Abandon hope, ye who value efficiency...

! This is a somewhat ugly hack aimed more at getting things up and running to
! test their value, and less at doing so in an elegant or efficient manner.
! Every processor does ALL of the linear algebra, instead of distributing it
! over the different processors.  If it works and is valuable, it may be worth
! the effort to rewrite this monstrosity to take advantage of the parallel
! nature of the code.

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine refine_initial_model (mype, npes)

   use heating_layer, only : do_refine
   use hydro_parameters, only : iu1, iu2
   use steady, only : ss0
   use variables, only : uin

   implicit none

   integer, intent(in) :: mype, npes

   integer :: idx

   if (do_refine) then
      call refine_actual(mype, npes)
!      call refine_vectest_norefine(mype, npes)
   end if

   ! Print steady state
   call write_ss_file(mype, npes)

   ! Set initial conditions
   do idx = iu1, iu2
      uin(1,idx,:,:,:) = 0.0d0
      uin(1,idx,:,:,1) = ss0(idx,1)
      uin(1,idx,:,:,2) = ss0(idx,2)
      uin(1,idx,:,:,5) = ss0(idx,3)
   end do

   return
end subroutine refine_initial_model

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

!subroutine refine_initial_model (mype, npes)
subroutine refine_actual (mype, npes)

   use amr_parameters, only : dp
   use hydro_parameters, only : iu1, iu2
   use initial_state_refinement, only : ISR_setup, ISR_cleanup, ISR_stop, &
                                      & ISR_N, ISR_npes
   use ISR_linalg, only : solve_Axb, invert_matrix, matvec, &
                        & Broyden_update_inverse, axpy, copy
   use steady, only : ss0
   use variables, only : uin

   implicit none

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   integer, intent(in) :: mype, npes

   ! Locals -------------------------------------------------------------------

   ! For calculations
   real(dp), allocatable, dimension(:,:) :: J, Jinv
   real(dp), allocatable, dimension(:) :: ss, dss, ss_old
   real(dp), allocatable, dimension(:) :: ff, dff, ff_old
   real(dp), allocatable, dimension(:) :: vv

   ! Convergence
   logical :: not_converged
   integer :: iter, iter_max
   integer :: loc
   real(dp) :: error, prev_error
   real(dp), parameter :: TOLERANCE = 1.0d-12
   real(dp) :: step
   real(dp), parameter :: min_step = 5.0d-2

   ! Misc
   integer :: idx
   integer :: row, col
   character(len=79) :: message, f_str

   ! ==========================================================================
   ! Set up

   call ISR_setup(mype, npes)

   allocate( ss    (ISR_N))
   allocate(dss    (ISR_N))
   allocate( ss_old(ISR_N))

   allocate( ff    (ISR_N))
   allocate(dff    (ISR_N))
   allocate( ff_old(ISR_N))

   allocate(J   (ISR_N, ISR_N*ISR_npes))
   allocate(Jinv(ISR_N, ISR_N*ISR_npes))

   allocate(vv(ISR_N))

   ! ==========================================================================
   ! Prepare for iteration

   if (mype == 0) then
      write(*,"('BEGINNING REFINEMENT ITERATION')")
   end if

   ! Construct a vector of the initial guess from Jerome's initial guess state
   ! NOTE: The way that ss0 is constructed does not guarantee that the guard
   !       cells are correct to machine precision.  Doing a pack here, then an
   !       unpack (inside refine_hydro_onestep) has the beneficial effect of
   !       correcting any discrepancies in the guard cells by using MPI
   !       communication to ensure that they are correct.
   call pack_vector(ss, ss0)

   ! Compute the change over a single time step
   call refine_hydro_onestep(ss, ff)

   call ISR_iteration_error(prev_error, 0, 0.0d0, ss, ff)

   ! Construct the initial Jacobian
   call construct_Jacobian(ss, ff, J)
   message = "Jacobian"
!call print_matrix(J, message)

   ! Invert the Jacobian
   call invert_matrix(Jinv, J)
   message = "Inverse Jacobian"
!call print_matrix(Jinv, message)

   ! ==========================================================================
   ! Refinement iteration
   ! -- A variation on Broyden's method (Broyden 1965) with some changes
   !    partially inspired by Kvaalen 1991.

   ! Loop
   not_converged = .true.
   iter_max = ISR_N*ISR_npes
   step = 1.0d0
   do iter = 1, iter_max

      ! Find corrections to steady state
      ! -- Taylor series: f(s*) = 0 ~ f(s) + J ds, so ds = - Jinv f
      call matvec(dss, Jinv, ff)
      dss(:) = -1.0d0 * step * dss(:)

      ! Cycle steady state
      call copy(ss_old, ss)
      call axpy(ss, 1.0d0, ss_old, dss)

      ! Cycle change in single time step
      call copy(ff_old, ff)
      call refine_hydro_onestep(ss, ff)
      call axpy(dff, -1.0d0, ff_old, ff)

      ! Check convergence
      call ISR_iteration_error(error, iter, step, ss, ff)
      if (error <= TOLERANCE) then
         not_converged = .false.
         exit
      end if

      ! Update the inverse of the Jacobian
      if (error < prev_error) then

         ! If we are converging, then all is well: keep using Broyden's method
         call Broyden_update_inverse(dss, dff, Jinv)

         ! Move back towards a full time step, because we may have escaped the
         ! "problem area" where we needed the small step
         step = 0.5d0 * (1.0d0 + step)

      else
         ! We are not converging, so we need to try something other than a
         ! simple Broyden's method; hopefully a smaller step will resolve the
         ! problem
         step = 0.5d0 * step

         ! If the problem is bad enough, we may have to recalculate the
         ! Jacobian from scratch (although we very strongly want to avoid that
         ! if possible)
         if (step >= min_step) then
            ! If the step is isn't too small, keep trying Broyden's method, but
            ! with a smaller step size
            call Broyden_update_inverse(dss, dff, Jinv)
         else
            ! If we are in a problem area (the step has dropped significantly),
            ! bite the bullet and recalculate the Jacobian from scratch; then
            ! return to Broyden's method with a full step
            call construct_Jacobian(ss, ff, J)
            call invert_matrix(Jinv, J)
            step = 1.0d0
         end if
      end if
      prev_error = error

   end do

   ! Verify convergence
   if (not_converged) then
      write(message,"('Unable to converge; error = ',es13.6,'.')") error
      call ISR_stop(message)
   end if

   ! Save final result to ss0
   call unpack_vector(ss0, ss, .true.)

   ! ==========================================================================
   ! Clean up

   deallocate(ss, dss, ss_old)
   deallocate(ff, dff, ff_old)
   deallocate(J, Jinv)
   deallocate(vv)

   call ISR_cleanup()

   return
end subroutine refine_actual
!end subroutine refine_initial_model

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine ISR_iteration_error (error, iter, step, ss, ff)

   use amr_parameters, only : dp
   use initial_state_refinement, only : ISR_N, glb_mype
   use variables, only : x
#  if WITHMPI==1
   use initial_state_refinement, only : ISR_comm, ISR_mype, ISR_npes
#  endif

   implicit none

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp),                   intent(inout) :: error
   integer ,                   intent(in   ) :: iter
   real(dp),                   intent(in   ) :: step
   real(dp), dimension(ISR_N), intent(in   ) :: ss, ff

   ! Locals -------------------------------------------------------------------
   real(dp), dimension(ISR_N) :: vv
   character(len=128) :: f_str
#  if WITHMPI==1
   integer, parameter :: ISR_ROOT = 0
   integer :: ierr
   real(dp) :: err_local, err_global
#  endif

   ! ==========================================================================
   ! Check convergence

   ! Determine the maximum error and its location
   vv = abs(ff(:) / ss(:))
#  if WITHMPI==1
   err_local = maxval(vv)
   call MPI_Allreduce(err_local, err_global, 1, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, ISR_comm, ierr)
   error = err_global
#  else
   error = maxval(vv)
#  endif

   ! Print the iteration error
   if (glb_mype == ISR_ROOT) then
      f_str = "('iteration ',i5"  // &
            & ",' error ',es13.6" // &
            & ",' step ',es13.6"  // &
            & ")"
      write(*,f_str) iter, error, step
   end if

   return
end subroutine ISR_iteration_error

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------
! NOTE: Because of the way ss0 is built, it is not guaranteed that the boundary
! conditions will be exact to machine precision; that is, ss(nx+1,var) on
! processor P is not guaranteed to be EXACTLY equal to ss(iu1,var) on processor
! P+1.  Therefore the guard cells may show a small amount of error in this
! routine.  That is a fault of the initial setup of ss0, but that is a small
! enough issue that we can safely ignore it.

subroutine refine_vectest_norefine (mype, npes)

   use amr_parameters, only : dp, nx
   use hydro_parameters, only : iu1, iu2
   use initial_state_refinement, only : ISR_setup, ISR_cleanup, &
                                      & ISR_N, ISR_mype, ISR_npes, ISR_comm
   use steady, only : ss0, SS_nv
   use variables, only : uin
#  if WITHMPI==1
   use mpi_var, only : nyslice, nzslice
#  endif

   implicit none

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   integer, intent(in) :: mype, npes

   ! Locals -------------------------------------------------------------------
   integer :: idx, var
   integer :: tmp1, tmp2, tmp3
   real(dp), allocatable :: vv(:)
   real(dp), dimension(iu1:iu2,SS_nv) :: ssf, ss_diff
   character(len=79) :: f_str, title
#  if WITHMPI==1
   integer :: ierr
#  endif

   ! ==========================================================================
   ! Set up

   call ISR_setup(mype, npes)

   allocate(vv(ISR_N))

   ! Write the pre-gather state
   write(title,"('Array on processor',x,i2,x,'/',x,i2)") ISR_mype+1, ISR_npes
   call print_state(ss0, .true., title)

   ! ==========================================================================
   ! Demonstrate a vector gather operation

   ! Gather ss0 to vv
   call pack_vector(vv, ss0)

   ! Write the vector
   write(title,"('Vector on processor',x,i2,x,'/',x,i2)") ISR_mype+1, ISR_npes
   call print_vector(vv, title)

   ! ==========================================================================
   ! Demonstrate a vector scatter operation

   ! Scatter vv to ssf
   ssf(:,:) = 0.0d0
   call unpack_vector(ssf, vv, .true.)

   ! Write the post-scatter state
   write(title,"('Array on processor',x,i2,x,'/',x,i2)") ISR_mype+1, ISR_npes
   call print_state(ssf, .true., title)

   ! ==========================================================================
   ! Verify that there is no change from the original to the final

   ! Calculate difference from ss0 to ssf
   do idx = iu1, iu2
      do var = 1, SS_nv
         ss_diff(idx,var) = abs(ssf(idx,var) - ss0(idx,var))
      end do
   end do

   ! Write the difference
   f_str = "('Difference on processor',x,i2,x,'/',x,i2)"
   write(title,f_str) ISR_mype+1, ISR_npes
   call print_state(ss_diff, .true., title)

   ! ==========================================================================
   ! Clean up

   deallocate(vv)

   call ISR_cleanup()

   return
end subroutine refine_vectest_norefine

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine print_state (vec, print_gc, title)

   use amr_parameters, only : dp, nx
   use hydro_parameters, only : iu1, iu2
   use steady, only : SS_nv
#  if WITHMPI==1
   use initial_state_refinement, only : ISR_mype, ISR_npes, &
                                      & ISR_left, ISR_right, ISR_comm
   use mpi_var, only : yposition, zposition
#  endif

   implicit none

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp)         , intent(in   ) :: vec(iu1:iu2,SS_nv)
   logical          , intent(in   ) :: print_gc
   character(len=79), intent(in   ) :: title

   ! Locals -------------------------------------------------------------------
   integer :: idx, ilo, ihi
   integer :: var
#  if WITHMPI==1
   integer :: junk, ierr
   integer :: stat(MPI_STATUS_SIZE)
#  endif

   ! ==========================================================================
   ! Write the vector state

   ! Only print one vertical slice
#  if WITHMPI==1
   call MPI_Barrier(ISR_comm, ierr)
   if ((yposition /= 0) .or. (zposition /= 0)) then
      return
   end if
#  endif

   ! Adjust loop limits depending on if guard cells will be printed
   if (print_gc) then
      ilo = iu1
      ihi = iu2
   else
      ilo = 1
      ihi = nx
   end if

   ! Wait for lower neighbor to finish
#  if WITHMPI==1
   call MPI_Barrier(ISR_comm, ierr)
   if (ISR_mype /= 0) then
      call MPI_Recv(junk, 1, MPI_INTEGER, ISR_left, 1, ISR_comm, stat, ierr)
   end if
#  endif

   ! Write
   write(*,"(a79)") title
   do idx = ilo, ihi
      do var = 1, SS_nv
         write(*,"(3x,es21.14)",advance='no') vec(idx,var)
      end do
      if (print_gc .and. (1 <= idx) .and. (idx <= nx)) then
         write(*,"(3x,a1)",advance='no') "]"
      end if
      write(*,*) ""
   end do

   ! Notify upper neighbor
#  if WITHMPI==1
   if (ISR_mype /= ISR_npes-1) then
      call MPI_Send(junk, 1, MPI_INTEGER, ISR_right, 1, ISR_comm, ierr)
   end if
#  endif

   return
end subroutine print_state

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine print_vector (vec, title)

   use amr_parameters, only : dp
   use initial_state_refinement, only : ISR_N
#  if WITHMPI==1
   use initial_state_refinement, only : ISR_mype, ISR_npes, &
                                      & ISR_left, ISR_right, ISR_comm
   use mpi_var, only : yposition, zposition
#  endif

   implicit none

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp)         , intent(in   ) :: vec(ISR_N)
   character(len=79), intent(in   ) :: title

   ! Locals -------------------------------------------------------------------
   integer :: idx, ilo, ihi
   integer :: var
#  if WITHMPI==1
   integer :: junk, ierr
   integer :: stat(MPI_STATUS_SIZE)
#  endif

   ! ==========================================================================
   ! Write the vector state

   ! Only print one vertical slice
#  if WITHMPI==1
   call MPI_Barrier(ISR_comm, ierr)
   if ((yposition /= 0) .or. (zposition /= 0)) then
      return
   end if
#  endif

   ! Wait for lower neighbor to finish
#  if WITHMPI==1
   call MPI_Barrier(ISR_comm, ierr)
   if (ISR_mype /= 0) then
      call MPI_Recv(junk, 1, MPI_INTEGER, ISR_left, 1, ISR_comm, stat, ierr)
   end if
#  endif

   ! Write
   write(*,"(a79)") title
   do idx = 1, ISR_N
      write(*,"(3x,es21.14)") vec(idx)
   end do

   ! Notify upper neighbor
#  if WITHMPI==1
   if (ISR_mype /= ISR_npes-1) then
      call MPI_Send(junk, 1, MPI_INTEGER, ISR_right, 1, ISR_comm, ierr)
   end if
#  endif

   return
end subroutine print_vector

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine print_matrix (M, title)

   use amr_parameters, only : dp
   use initial_state_refinement, only : ISR_N
#  if WITHMPI==1
   use initial_state_refinement, only : ISR_mype, ISR_npes, &
                                      & ISR_left, ISR_right, ISR_comm
   use mpi_var, only : yposition, zposition
#  endif

   implicit none

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp)         , intent(in   ) :: M(ISR_N, ISR_N*ISR_npes)
   character(len=79), intent(in   ) :: title

   ! Locals -------------------------------------------------------------------
   integer :: row, col
#  if WITHMPI==1
   integer :: junk, ierr
   integer :: stat(MPI_STATUS_SIZE)
#  endif

   ! ==========================================================================
   ! Write the vector state

   ! Only print one vertical slice
#  if WITHMPI==1
   call MPI_Barrier(ISR_comm, ierr)
   if ((yposition /= 0) .or. (zposition /= 0)) then
      return
   end if
#  endif

   ! Wait for lower neighbor to finish
#  if WITHMPI==1
   call MPI_Barrier(ISR_comm, ierr)
   if (ISR_mype /= 0) then
      call MPI_Recv(junk, 1, MPI_INTEGER, ISR_left, 1, ISR_comm, stat, ierr)
   end if
#  endif

   ! Write
   write(*,"(a79)") title
   do row = 1, ISR_N
      do col = 1, ISR_N*ISR_npes
         write(*,"(3x,es21.14)",advance="no") M(row,col)
      end do
      write(*,*) ""
   end do

   ! Notify upper neighbor
#  if WITHMPI==1
   if (ISR_mype /= ISR_npes-1) then
      call MPI_Send(junk, 1, MPI_INTEGER, ISR_right, 1, ISR_comm, ierr)
   end if
   call MPI_Barrier(ISR_comm, ierr)
#  endif

   return
end subroutine print_matrix

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------
! The vector gather operation collects all components from across the slice and
! packs them into a 1D array by index-major/variable-minor ordering.  It
! neglects boundary conditions, and only includes body cells.

subroutine pack_vector (dst, src)

   use amr_parameters, only : dp, nx
   use hydro_parameters, only : iu1, iu2
   use initial_state_refinement, only : ISR_N
   use steady, only : SS_nv
#  if WITHMPI==1
   use initial_state_refinement, only : ISR_comm, ISR_npes
#  endif

   implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp), intent(  out) :: dst(ISR_N)
   real(dp), intent(in   ) :: src(iu1:iu2, SS_nv)

   ! Locals -------------------------------------------------------------------
   integer :: idx, var

   ! ==========================================================================
   ! Pack the vector from the state array in index-major ordering

   do idx = 1, nx
      do var = 1, SS_nv
         dst((idx-1)*SS_nv+var) = src(idx,var)
      end do
   end do

   return
end subroutine pack_vector

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine unpack_vector (dst, src, is_ss0)

   use amr_parameters, only : dp, nx
   use hydro_parameters, only : iu1, iu2
   use initial_state_refinement, only : ISR_N, ss_bc, ISR_mype, ISR_npes
   use steady, only : SS_nv
#  if WITHMPI==1
   use initial_state_refinement, only : ISR_left, ISR_right, ISR_nguard
   use mpi_var, only : xposition
#  endif

   implicit none

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp), intent(inout) :: dst(iu1:iu2, SS_nv)
   real(dp), intent(in   ) :: src(ISR_N)
   logical , intent(in   ) :: is_ss0

   ! Locals -------------------------------------------------------------------
   integer :: idx, var
#  if WITHMPI==1
   real(dp), dimension(ISR_nguard*SS_nv) :: sendbuf_lo, recvbuf_lo
   real(dp), dimension(ISR_nguard*SS_nv) :: sendbuf_hi, recvbuf_hi
   integer, dimension(MPI_STATUS_SIZE,4) :: statuses
   integer, dimension(4) :: requests
   integer, parameter :: tag_pass_up   = 15
   integer, parameter :: tag_pass_down = 16
   integer :: ierr
#  endif

   ! ==========================================================================
   ! Unpack the vector into the state array

   do idx = 1, nx
      do var = 1, SS_nv
         dst(idx,var) = src((idx-1)*SS_nv + var)
      end do
   end do

   ! ==========================================================================
   ! Fill the guard cells

#  if WITHMPI==1
   ! MPI Communication - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Pack send buffers
   do idx = 0, ISR_nguard-1
      do var = 1, SS_nv
         if (ISR_mype /= 0) then
            sendbuf_lo(idx*SS_nv+var) = dst(iu1+  ISR_nguard+  idx,var)
         end if
         if (ISR_mype /= ISR_npes-1) then
            sendbuf_hi(idx*SS_nv+var) = dst(iu2-2*ISR_nguard+1+idx,var)
         end if
      end do
   end do

   ! Asynchronous receives
   if (ISR_mype /= 0) then
      call MPI_Irecv(recvbuf_lo, ISR_nguard*SS_nv, MPI_DOUBLE_PRECISION,    &
                     ISR_left , tag_pass_up  , MPI_COMM_WORLD, requests(1), &
                     ierr)
   end if
   if (ISR_mype /= ISR_npes-1) then
      call MPI_Irecv(recvbuf_hi, ISR_nguard*SS_nv, MPI_DOUBLE_PRECISION,    &
                     ISR_right, tag_pass_down, MPI_COMM_WORLD, requests(2), &
                     ierr)
   end if

   ! Asynchronous send
   if (ISR_mype /= 0) then
      call MPI_Isend(sendbuf_lo, ISR_nguard*SS_nv, MPI_DOUBLE_PRECISION,    &
                     ISR_left , tag_pass_down, MPI_COMM_WORLD, requests(3), &
                     ierr)
   end if
   if (ISR_mype /= ISR_npes-1) then
      call MPI_Isend(sendbuf_hi, ISR_nguard*SS_nv, MPI_DOUBLE_PRECISION,    &
                     ISR_right, tag_pass_up  , MPI_COMM_WORLD, requests(4), &
                     ierr)
   end if

   ! Wait for receives to finish (don't worry about sends)
   ! NOTE: Only wait for the receives that are started.  But because of the way
   !       Fortran MPI works, you can't just do e.g. requests(2) because that
   !       is a real instead of a pointer to a real array.  Thus you have to do
   !       requests(2:2) to get a one-element array.
   if (ISR_mype == 0) then
      call MPI_Waitall(1, requests(2:2), statuses(:,2:2), ierr)
   else if (ISR_mype == ISR_npes-1) then
      call MPI_Waitall(1, requests(1:1), statuses(:,1:1), ierr)
   else
      call MPI_Waitall(2, requests(1:2), statuses(:,1:2), ierr)
   end if
!   if (ISR_mype /= 0) then
!      call MPI_Wait(requests(1), statuses(:,1), ierr)
!   end if
!   if (ISR_mype /= ISR_npes-1) then
!      call MPI_Wait(requests(2), statuses(:,2), ierr)
!   end if
#  endif

   ! Filling the guard cells - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Lower guard cells
   if (ISR_mype == 0) then
      ! Copy values at idx == 1
      do idx = iu1, 0
         do var = 1, SS_nv
            dst(idx,var) = dst(1,var)
         end do
      end do
   else
      ! Unpack receive buffer
      do idx = 0, ISR_nguard-1
         do var = 1, SS_nv
            dst(iu1+idx,var) = recvbuf_lo(idx*SS_nv+var)
         end do
      end do
   end if

   ! Upper guard cells
   if (ISR_mype == ISR_npes-1) then
      if (is_ss0) then
         ! Get values from fixed inflow condition
         do idx = nx+1, iu2
            do var = 1, SS_nv
               dst(idx,var) = ss_bc(var)
            end do
         end do
      else
         ! Copy values at idx == nx
         do idx = nx+1, iu2
            do var = 1, SS_nv
               dst(idx,var) = dst(nx,var)
            end do
         end do
      end if
   else
      ! Unpack receive buffer
      do idx = 0, ISR_nguard-1
         do var = 1, SS_nv
            dst(iu2-ISR_nguard+1+idx,var) = recvbuf_hi(idx*SS_nv+var)
         end do
      end do
   end if

   return
end subroutine unpack_vector

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine refine_hydro_onestep (ss, ff)

   use amr_parameters, only : dp
   use hydro_parameters, only : iu1, iu2
   use initial_state_refinement, only : ISR_N, glb_mype, glb_npes, ISR_stop
   use steady, only : SS_nv
   use variables, only : uin, old_uin, gravin, &
                       & flux, flux_pre_tot, &
                       & emfx, emfy, emfz, &
                       & bval_type, fargo, nu, eta, &
                       & dx, dy, dz, ngrid, &
                       & time, dt, courant

   implicit none

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp), intent(in   ) :: ss(ISR_N)
   real(dp), intent(  out) :: ff(ISR_N)

   ! Locals -------------------------------------------------------------------
   real(dp), dimension(iu1:iu2,SS_nv) :: ss_local, ff_local
   integer :: idx

   ! ==========================================================================
   ! Distribute steady state and load the initial conditions to uin

   call unpack_vector(ss_local, ss, .true.)

   uin(:,:,:,:,:) = 0.0d0
   do idx = iu1, iu2
      uin(:,idx,:,:,1) = ss_local(idx,1)
      uin(:,idx,:,:,2) = ss_local(idx,2)
      uin(:,idx,:,:,5) = ss_local(idx,3)
   end do

   ! ==========================================================================
   ! Compute a single hydro step (modified from dumses.f90)

   ! Reset time
   time = 0.0d0

   ! Compute time step
   call compute_dt(uin, dx, dy, dz, dt, courant, ngrid)

   ! Run the MUSCL algorithm
   call umuscl(uin, gravin, flux, flux_pre_tot, emfx, emfy, emfz, &
              & dx, dy, dz, dt, ngrid)

   ! Shearing
   if (bval_type(1,1) == 3) then
      ! should be: call bval_shear_flux(flux, emfy, time+dt/2.,dy, ngrid, mype)
      ! not using a shearing box --- just skip this
      call ISR_stop("Bad BC in refine: shearing box not allowed.")
   end if

   ! Copy to old_uin
   old_uin = uin

   ! Update
   call update(uin, flux, flux_pre_tot, emfx, emfy, emfz, bval_type, ngrid, dt)

   ! Sources
   call source_term(glb_mype, glb_npes)

   ! Rotation
   if (fargo) then
      ! should be: call fargo_update(uin, emfx, emfz, dx, dy, dz, dt, ngrid, mype)
      ! not doing rotating disks --- just skip this
      call ISR_stop("Bad option in refine: fargo not allowed.")
   end if

   ! BC fill
   ! should be: call bval(uin, ngrid, time, dt, dx, dy, dz, bval_type, mype)
   ! don't need a BC fill here; I do that elsewhere --- just skip

   ! Dissipative processes
   if ((nu > 0.0d0) .or. (eta > 0.0d0)) then
      ! should be: call dissip(uin, emfx, emfy, emfz, ngrid, mype)
      ! not doing resistivity or viscosity --- just skip this
      call ISR_stop("Bad setting in refine: dissipation not allowed.")
   end if

   ! ==========================================================================
   ! Compute change and gather

   ! Output is the absolute change in each cell
   ff_local(:,1) = (uin(1,:,1,1,1) - ss_local(:,1))
   ff_local(:,2) = (uin(1,:,1,1,2) - ss_local(:,2))
   ff_local(:,3) = (uin(1,:,1,1,5) - ss_local(:,3))

   ! Gather
   call pack_vector(ff, ff_local)

   return
end subroutine refine_hydro_onestep

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------
! Assuming that the step will move somewhat in the direction of ff (as that is
! the correction over a single time step), we can try to build an asymmetric
! finite-difference, where we take ss and ss+max(ff,max_step*ss), which limits
! the step size to not take a huge step (e.g. in case dt is really big for some
! reason).

subroutine construct_Jacobian (ss, ff, J)

   use amr_parameters, only : dp
   use initial_state_refinement, only : ISR_N, ISR_stop, ISR_npes, ISR_mype, &
                                        glb_mype
   use ISR_linalg, only : copy
   use steady, only : SS_nv

   implicit none

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp), dimension(ISR_N)               , intent(in   ) :: ss, ff
   real(dp), dimension(ISR_N,ISR_N*ISR_npes), intent(inout) :: J

   ! Locals -------------------------------------------------------------------
   real(dp), parameter :: min_eps = 1.0d-6
   real(dp), parameter :: max_eps = 1.0d-3
   real(dp) :: delta
   real(dp), dimension(ISR_N) :: ss_p, ff_p
   integer :: proc, row, col
   real(dp) :: min_on, max_off
   real(dp), parameter :: SMALL_J = 1.0d-6
   character(len=75) :: f_str, message
#  if WITHMPI==1
   integer :: ierr
#  endif

   ! ==========================================================================
   ! Construct a finite-difference estimate for each element

   do proc = 0, ISR_npes-1
      do col = 1, ISR_N

         if (glb_mype == 0) then
            write(*,"(3x,'construct Jacobian: column ',i5)") proc*ISR_N+col
         end if

         ! Construct perturbed input state
         ! -- ss_p = ss, except for one element
         ! -- that element is perturbed by ff (which is the change in ss over
         !    one step of the hydro routine starting from ss)
         ! -- the perturbation is limited to be no more than some fraction of
         !    the value in ss (by max_eps)
         call copy(ss_p, ss)
         if (proc == ISR_mype) then
            delta = max(min_eps*ss(col), min(max_eps*ss(col), ff(col)))
            ss_p(col) = ss_p(col) + delta
         end if

#        if WITHMPI==1
         ! Broadcast delta
         call MPI_Bcast(delta, 1, MPI_DOUBLE_PRECISION, &
                        proc, MPI_COMM_WORLD, ierr)
#        endif

         ! Compute output state from perturbed input state
         call refine_hydro_onestep(ss_p, ff_p)

         ! Construct finite-difference approximations for current column
         ! -- J = (f(x+h) - f(x)) / h
         J(:,proc*ISR_N+col) = (ff_p(:) - ff(:)) / delta

      end do
   end do

   return
end subroutine construct_Jacobian

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine write_ss_file (mype, npes)

   use amr_parameters, only : nx
   use hydro_parameters, only : iu1, iu2
#  if WITHMPI==1
   use mpi_var,          only : xposition, yposition, zposition
   use initial_state_refinement, only : ISR_comm, ISR_npes, ISR_mype
#  endif
   use steady,           only : ss0, SS_nv
   use variables,        only : x

   implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   integer, intent(in) :: mype, npes

   ! Locals -------------------------------------------------------------------

   integer :: proc, idx, var
   integer, parameter :: LUN = 16
   character(len=*), parameter :: fname = "stable_state.dat"
   integer :: ierr

   ! ==========================================================================
   ! Write the stable state

   ! Write the header (only the master processor does this)
#  if WITHMPI==1
   if (mype == 0) then
#  endif
      open(unit=LUN, file=fname, status='unknown')
      write(LUN,"('#',x,4(3x,a30))") "radius", &
         "density", "x-momentum", "total_energy_dens"
      close(LUN)
#  if WITHMPI==1
   end if
   ! synchronise - only one processor writes at a time
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
#  endif

   ! Write the data in order (lowest value to highest)
#  if WITHMPI==1
   do proc = 0, ISR_npes-1
      if (proc == ISR_mype) then
         if ((yposition == 0) .and. (zposition == 0)) then
#  endif
            open(unit=LUN, file=fname, access='append', status='old')
            do idx = 1, nx
               write(LUN,"(5x,es30.23)",advance="no") x(idx)
               do var = 1, SS_nv
                  write(LUN,"(3x,es30.23)",advance="no") ss0(idx,var)
               end do
               write(LUN,*) ""
            end do
            close(LUN)
#  if WITHMPI==1
         end if
      end if
      ! synchronise - wait your turn to write so file is in order
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
   end do
#  endif

end subroutine write_ss_file

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

