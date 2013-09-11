!=============================================================================!
! This subroutine computes an initial guess for a cell based on a Runge-Kutta !
! integration.                                                                !
!                                                                             !
! Jerome's method is so poorly-written (in terms of style and comments), that !
! I can't make heads or tails of what it's doing, and I'm hesitant to blindly !
! copy it.  Since I don't need the exact method in this particular case, I'm  !
! just going to code up a basic 4th-order RK method.                          !
!=============================================================================!

subroutine compute_guess(i)

   use hydro_parameters
   use variables

   implicit none

   ! ==========================================================================
   ! Declare variables

   real(dp) :: x0
   real(dp) :: dx0
   integer, parameter :: Nv = 2
   real(dp) :: y0(Nv)   ! 1: v; 2: c
   real(dp) :: dydx1(Nv), dydx2(Nv), dydx3(Nv), dydx4(Nv)

   ! ==========================================================================
   ! Just do RK4

   ! Integrate from cell i-1 to cell i
   x0 = x(i-1)
   y0(1) = uin(i-1,2) / uin(i-1,1)
   y0(2) = sqrt(gamma * (gamma-1.0d0) * (uin(i-1,5)/uin(i-1,1) - 0.5*u0**2))
   dx0 = x(i) - x0

   ! Construct k-values
   xx = x0
   yy = y0
   call derivatives(xx, yy, dydx1)
   xx = x0 + 0.5d0*dx0
   yy = y0 + 0.5d0*dx0*dydx1
   call derivatives(xx, yy, dydx2)
   xx = x0 + 0.5d0*dx0
   yy = y0 + 0.5d0*dx0*dydx2
   call derivatives(xx, yy, dydx3)
   xx = x0 +       dx0
   yy = y0 +       dx0*dydx3
   call derivatives(xx, yy, dydx4)

   ! Construct final values
   yy = y0 + (dydx1 + 2.0d0*dydx2 + 2.0d0*dydx3 + dydx4) * dx0 / 6.0d0

   ! Save the results
   uin(i,1) = rho_u / y0(1)
   uin(i,2) = rho_u
   uin(i,3) = 0.0d0
   uin(i,4) = 0.0d0
   uin(i,5) = (y0(2)**2 / (gamma*(gamma-1.0d0)*y0(1)) + 0.5 * y0(1)) * rho_u

end subroutine compute_guess
