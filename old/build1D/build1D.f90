!=============================================================================!
! This program mimics the DUMSES core in order to iteratively create an       !
! initial model very close to a steady state.                                 !
!                                                                             !
! I strip out a lot of functionality of DUMSES (especially the MPI as the     !
! iterations have to go one cell at a time as a projection method), but copy  !
! the hydro update core in order to match the same numerical errors.          !
!=============================================================================!

program build1D

   implicit none

   ! ==========================================================================
   ! Declare variables

   integer :: idx, iter
   integer, parameter :: MAX_ITER = 500

   ! ==========================================================================
   ! Initialize

   call init()

   ! ==========================================================================
   ! Create the initial model

   ! Loop over all cells
   do idx = 1, nx

      ! Use Jerome's method to create an initial guess for cell idx+1
      call compute_guess(idx+1)

      ! Iteratively refine the values of cell idx+1
      ! --> In order to achieve a steady state for cell idx
      do iter = 1, MAX_ITER

         ! Do a single-step hydro update on cell idx
         call hydro_onestep_cell(idx)

         ! Check if the change to cell idx is sufficiently small to declare
         ! cell idx+1 to be converged

         ! Compute the change vector for cell idx

         ! Invert the Jacobian to find the update to interface state idx+1/2

         ! Apply the update to cell idx+1
         ! --> Assuming that the reconstruction of the interface state is
         !     monotonic (including the creation of the Riemann input L & R
         !     states, as well as the Riemann solve itself), and assuming that
         !     the derivative of the interface state (idx+1/2) with respect to
         !     the idx+1 state is between 0 and 1 (both assumptions I think are
         !     true but I haven't proven either), then applying the update to
         !     cell idx+1 instead of to the interface means that we are
         !     effectively updating the interface state in the correct
         !     direction but with too small of a magnitude.  The wrong
         !     direction would not work, and too large of a magnitude may
         !     overshoot and oscillate, but too small of a magnitude simply
         !     implies extra steps are necessary to get to the correct
         !     solution.  Hence we cheat and apply the update to cell idx+1
         !     even though formally the update to be applied to the interface
         !     state.
         ! --> The other option is to attempt to modify the values in cell
         !     idx+1 in such a way as to update the interface state by the
         !     desired values, but that is challenging as the Riemann solve and
         !     the slope limiter can be changed by runtime parameters so that
         !     we cannot easily write down an expression for the interface
         !     state as a function of cells idx and idx+1.

      ! End iteration loop

   ! End the cell traversal loop
   end do

   ! ==========================================================================
   ! Print results

   ! Write model to a data file

   ! Write a parameter file for starting a DUMSES simulation

end program

