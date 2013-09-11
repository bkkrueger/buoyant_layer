!=============================================================================!
! This subroutine performs a single hydro step on a single specified cell.    !
!                                                                             !
! This subroutine packages up alternate versions of the routines that DUMSES  !
! uses to advance a single timestep, adapted to my specific case (no MPI, 1D, !
! only advance a single cell at a time, etc).                                 !
!=============================================================================!

subroutine hydro_onestep_cell(idx)

   use hydro_parameters
   use variables

   implicit none

   ! ==========================================================================
   ! Declare variables

   ! ==========================================================================
   ! DUMSES algorithm

   ! Get the time step
   call compute_dt(uin(idx,1), uin(idx,2), uin(idx,5), dx, dt, courant)

   ! MUSCL algorithm
   call umuscl(uin, gravin, flux, flux_pre_tot, dx, dt)

   if (bval_type(1) == 3) then
      ! TODO --- not doing shearing box, so throw an error here
   end if

   old_uin = uin  ! TODO : better way to save this result

