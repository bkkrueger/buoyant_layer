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
      !call refine_vectest_norefine(mype, npes)
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
                                      & ISR_N
   use ISR_linalg, only : solve_Axb, invert_matrix, matvec, &
                        & Broyden_update_inverse
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

   allocate(J   (ISR_N, ISR_N))
   allocate(Jinv(ISR_N, ISR_N))

   allocate(vv(ISR_N))

   ! ==========================================================================
   ! Prepare for iteration

   ! Construct the initial guess
   ! -- This is done by merging all the components of ss0 from across the
   !    processors into a single local array
   call gather_vector(ss0, ss)

   ! Compute the change over a single time step
   call refine_hydro_onestep(ss, ff)

   vv         = abs(ff(:) / ss(:))
   prev_error = maxval(vv)
   loc        = maxloc(vv, 1)
   if (mype == 0) then
      f_str = "('iteration ',i5"  // &
            & ",' error ',es13.6" // &
            & ",' step ',es13.6"  // &
            & ",' location ',i5"  // &
            & ")"
      write(*,f_str) 0, prev_error, 0.0d0, loc
   end if

   ! Construct the initial Jacobian
   call construct_Jacobian(ss, ff, J)

   ! Invert the Jacobian
   call invert_matrix(Jinv, J)

   ! ==========================================================================
   ! Refinement iteration
   ! -- A variation on Broyden's method (Broyden 1965) with some changes
   !    partially inspired by Kvaalen 1991.

   ! Loop
   not_converged = .true.
   iter_max = ISR_N
   step = 1.0d0
   do iter = 1, iter_max

      ! Find corrections to steady state
      ! -- Taylor series: f(s*) = 0 ~ f(s) + J ds, so ds = - Jinv f
      call matvec(dss, Jinv, ff)
      dss(:) = -1.0d0 * step * dss(:)

      ! Cycle steady state
      ss_old(:) = ss(:)
      ss    (:) = ss(:) + dss(:)

      ! Cycle change in single time step
      ff_old(:) = ff(:)
      call refine_hydro_onestep(ss, ff)
      dff(:) = ff(:) - ff_old(:)

      ! Check convergence
      vv    = abs(ff(:) / ss(:))
      error = maxval(vv)
      loc   = maxloc(vv, 1)
      if (mype == 0) then
         f_str = "('iteration ',i5"  // &
               & ",' error ',es13.6" // &
               & ",' step ',es13.6"  // &
               & ",' location ',i5"  // &
               & ")"
         write(*,f_str) iter, error, step, loc
      end if
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
   call scatter_vector(ss, ss0, .true.)

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

subroutine refine_vectest_norefine (mype, npes)

   use amr_parameters, only : dp, nx
   use hydro_parameters, only : iu1, iu2
   use initial_state_refinement, only : ISR_setup, ISR_cleanup, &
                                      & ISR_N, ISR_mype, ISR_npes
   use steady, only : ss0, SS_nv
   use variables, only : uin

   implicit none

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   integer, intent(in) :: mype, npes

   ! Locals -------------------------------------------------------------------
   integer :: idx, var
   integer :: tmp1, tmp2, tmp3
   real(dp), allocatable :: ss_in(:)
   real(dp), dimension(iu1:iu2,SS_nv) :: ss_out, ss_diff
   character(len=79) :: f_str, title

   ! ==========================================================================
   ! Set up

   call ISR_setup(mype, npes)

   allocate(ss_in(ISR_N))

   ! Write the pre-gather state
   write(title,"('Array on processor',x,i2,x,'/',x,i2)") ISR_mype+1, ISR_npes
   call print_state(ss0, .true., title)

   ! ==========================================================================
   ! Demonstrate a vector gather operation

   ! Gather ss0 to ss_in
   call gather_vector(ss0, ss_in)

   ! Write the post-gather state (ss_in)
   if (mype == 0) then
      f_str = "('Gathered array on processor',x,i2,x,'/',x,i2)"
      write(*,f_str) ISR_mype+1, ISR_npes
      do idx = 1, ISR_N
         tmp3 = mod(idx, SS_nv)
         tmp1 = (idx - tmp3) / SS_nv
         tmp2 = mod(tmp1, nx)
         tmp1 = (idx - tmp2) / nx
         tmp2 = tmp2 - 1
         tmp1 = tmp1 / nx
         write(*,"(4(2x,i3),2x,es13.6)") idx, tmp1, tmp2, tmp3, ss_in(idx)
      end do
   end if

   ! ==========================================================================
   ! Demonstrate a vector scatter operation

  ! Scatter ss_in to ss_out
   call scatter_vector(ss_in, ss_out, .true.)

   ! Write the post-scatter state (ss_out)
   write(title,"('Array on processor',x,i2,x,'/',x,i2)") ISR_mype+1, ISR_npes
   call print_state(ss_out, .true., title)

   ! ==========================================================================
   ! Verify that there is no change from the original to the final

   ! Calculate difference between ss_out and ss0
   do idx = iu1, iu2
      do var = 1, SS_nv
         ss_diff(idx,var) = abs(ss_out(idx,var) - ss0(idx,var))
      end do
   end do

   ! Write the difference (ss_diff)
   f_str = "('Difference on processor',x,i2,x,'/',x,i2)"
   write(title,f_str) ISR_mype+1, ISR_npes
   call print_state(ss_diff, .true., title)

   ! ==========================================================================
   ! Clean up

   deallocate(ss_in)

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

#     if WITHMPI==1
#     include "mpif.h"
#     endif

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
         write(*,"(3x,es13.6)",advance='no') vec(idx,var)
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
! The vector gather operation collects all components from across the slice and
! packs them into a 1D array by index-major/variable-minor ordering.  It
! neglects boundary conditions, and only includes body cells.

subroutine gather_vector (src, dst)

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
   real(dp), intent(in   ) :: src(iu1:iu2, SS_nv)
   real(dp), intent(  out) :: dst(ISR_N)

   ! Locals -------------------------------------------------------------------
   integer :: idx, var
#  if WITHMPI==1
   integer :: ierr
   real(dp) :: sendbuf(         nx*SS_nv)
   real(dp) :: recvbuf(ISR_npes*nx*SS_nv)
#  endif

   ! ==========================================================================
   ! Gather the vector from across the processors

#  if WITHMPI==1
   ! Pack the send buffer
   do idx = 1, nx
      do var = 1, SS_nv
         sendbuf((idx-1)*SS_nv+var) = src(idx,var)
      end do
   end do

   ! Allgather (as noted above, the packing order is such that the Allgather
   ! correctly packs the full array)
   call MPI_Allgather(sendbuf, nx*SS_nv, MPI_DOUBLE_PRECISION, &
                      recvbuf, nx*SS_nv, MPI_DOUBLE_PRECISION, ISR_comm, ierr)

   ! Unpack the receive buffer
   dst(:) = recvbuf(:)
#  else
   ! Pack into index-major ordering
   do idx = 1, nx
      do var = 1, SS_nv
         dst((idx-1)*SS_nv+var) = src(idx,var)
      end do
   end do
#  endif

   return
end subroutine gather_vector

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine scatter_vector (src, dst, is_ss0)

   use amr_parameters, only : dp, nx
   use hydro_parameters, only : iu1, iu2
   use initial_state_refinement, only : ISR_N, ss_bc
   use steady, only : SS_nv
#  if WITHMPI==1
   use initial_state_refinement, only : ISR_mype, ISR_npes
   use mpi_var, only : xposition
#  endif

   implicit none

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp), intent(in   ) :: src(ISR_N)
   real(dp), intent(  out) :: dst(iu1:iu2, SS_nv)
   logical , intent(in   ) :: is_ss0

   ! Locals -------------------------------------------------------------------
   integer :: offset
   integer :: idx, var

   ! ==========================================================================
   ! Extract the appropriate section of the vector

   ! Calculate offset
#  if WITHMPI==1
   offset = xposition * nx * SS_nv
#  else
   offset = 0
#  endif

   do idx = 1, nx
      do var = 1, SS_nv
         dst(idx,var) = src(offset + (idx-1)*SS_nv + var)
      end do
   end do

   ! ==========================================================================
   ! Fill the guard cells

   ! Lower guard cells
#  if WITHMPI==1
   if (ISR_mype /= 0) then
      do idx = 0, iu1, -1
         do var = 1, SS_nv
            dst(idx,var) = src(offset + (idx-1)*SS_nv + var)
         end do
      end do
   else
#  endif
      do idx = 0, iu1, -1
         do var = 1, SS_nv
            dst(idx,var) = src(offset + var) ! copy values at idx == 1
         end do
      end do
#  if WITHMPI==1
   end if
#  endif

   ! Upper guard cells
#  if WITHMPI==1
   if (ISR_mype /= ISR_npes-1) then
      do idx = nx+1, iu2
         do var = 1, SS_nv
            dst(idx,var) = src(offset + (idx-1)*SS_nv + var)
         end do
      end do
   else
#  endif
      if (is_ss0) then
         do var = 1, SS_nv
            dst(nx+1:iu2,var) = ss_bc(var)
         end do
      else
         do idx = nx+1, iu2
            do var = 1, SS_nv
               ! copy values at idx == nx
               dst(idx,var) = src(offset + (nx-1)*SS_nv + var)
            end do
         end do
      end if
#  if WITHMPI==1
   end if
#  endif

   return
end subroutine scatter_vector

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

   call scatter_vector(ss, ss_local, .true.)

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
   call gather_vector(ff_local, ff)

   return
end subroutine refine_hydro_onestep

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine construct_Jacobian (ss, ff, J)

   use amr_parameters, only : dp
   use initial_state_refinement, only : ISR_N, ISR_stop
   use steady, only : SS_nv

   implicit none

   ! ==========================================================================
   ! Declare variables

   ! Arguments ----------------------------------------------------------------
   real(dp), dimension(ISR_N)      , intent(in   ) :: ss, ff
   real(dp), dimension(ISR_N,ISR_N), intent(inout) :: J

   ! Locals -------------------------------------------------------------------
   real(dp), parameter :: eps = 0.5d-3  ! +/-0.5% --> difference of 1%
   real(dp) :: delta
   real(dp), dimension(ISR_N) :: ss_p, ss_m, ff_p, ff_m
   integer :: row, col
   real(dp) :: min_on, max_off
   real(dp), parameter :: SMALL_J = 1.0d-6
   character(len=75) :: f_str, message

   ! ==========================================================================
   ! Construct a finite-difference estimate for each element

   do col = 1, ISR_N

      delta = eps * ss(col)

      ! Construct perturbed input states
      ss_p( : ) = ss( : )
      ss_p(col) = ss(col) + delta
      ss_m( : ) = ss( : )
      ss_m(col) = ss(col) - delta

      ! Compute output states from perturbed input states
      call refine_hydro_onestep(ss_p, ff_p)
      call refine_hydro_onestep(ss_m, ff_m)

      ! Construct finite-difference approximations for current column
      J(:,col) = (ff_p(:) - ff_m(:)) / (2.0d0 * delta)

   end do

   return
end subroutine construct_Jacobian

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine write_ss_file (mype, npes)

   use hydro_parameters, only : iu1, iu2
#  if WITHMPI==1
   use mpi_var,          only : xposition, yposition, zposition, nxslice
   use initial_state_refinement, only : ISR_comm
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

   integer :: xpos, idx, var
   integer, parameter :: LUN = 16
   character(len=*), parameter :: fname = "stable_state.dat"
   integer :: ierr

   ! ==========================================================================
   ! Write the stable state

! TODO --- I could change this to do a send-receive chain instead of using
!          MPI_Barrier, but which is better?

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
   do xpos = 0, nxslice-1
      if (xposition == xpos) then
         if ((yposition == 0) .and. (zposition == 0)) then
#  endif
            open(unit=LUN, file=fname, access='append', status='old')
            do idx = iu1+3, iu2-3
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

