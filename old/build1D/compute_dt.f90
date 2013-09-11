!=============================================================================!
! This subroutine replaces the DUMSES compute_dt subroutine.                  !
!                                                                             !
! Because I am only advancing a single cell, in 1D, on a single-layer grid, I !
! can simplify and/or neglect a number of arguments.                          !
!=============================================================================!

subroutine compute_dt(dens_in, momx_in, Ener_in, dx, dt, courant)

   use hydro_parameters
   use variables, only : x, cartesian, cylindrical, spherical, fargo

   implicit none

   ! ==========================================================================
   ! Declare variables

   ! Arguments
   real(dp), intent(in   ) :: dens_in, momx_in, Ener_in
   real(dp), intent(in   ) :: dx, courant
   real(dp), intent(  out) :: dt

   ! Locals
   real(dp) :: dens, velx, Ekin, pres
   real(dp) :: vx_max

   ! ==========================================================================
   ! Compute the time step

   ! Get primitives
   dens = max(dens_in, smallr)
   velx = momx_in / dens
   Ekin = 0.5 * dens * velx**2
   pres = max(Ener_in - Ekin) * (gamma - 1.0d0), smallp)

   ! Compute fastest signal speed
   vx_max = sqrt(gamma*pres/dens) + abs(velx)

   ! Compute dt
   dt = courant * min(courant*dx/smallc, dx/vx_max)

   return
end subroutine compute_dt

