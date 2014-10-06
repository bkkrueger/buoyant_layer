! =============================================================================

module buoyant_layer

   use amr_parameters, only : dp

   implicit none

   ! Heating layer
   character(len=20) :: layer_shape    ! The shape of the buoyant layer
   real(dp)          :: layer_limit    ! The width of the buoyant layer
   real(dp)          :: Kheat          ! The strength of heating
   real(dp)          :: heat_coef      ! Scaling constant
   real(dp)          :: Kgrav          ! The strength of gravity
   real(dp)          :: grav_coef      ! Scaling constant

contains

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   double precision function shape_function(x0)
      use amr_parameters, only : dp
      implicit none
      real(dp) :: x0

      ! Trapezoid:
      ! - f(x) = 1 if |x| < H/2
      ! - f(x) = 0 if |x| > H
      ! - linear connection between 0 and 1 if H/2 < |x| < H
      if (layer_shape .eq. 'trapezoid') then
         if (abs(x0) > layer_limit) then
            shape_function = 0.0d0
         else if (abs(x0) <= 0.5d0*layer_limit) then
            shape_function = 1.0d0
         else
            shape_function = 2.0d0 * (1.0d0-abs(x0)/layer_limit)
         end if

      ! Square:
      ! - f(x) = 1 if |x| < H
      ! - f(x) = 0 if |x| > H
      else if (layer_shape .eq. 'square') then
         if (abs(x0) > layer_limit) then
            shape_function = 0.0d0
         else
            shape_function = 1.0d0
         end if

      ! default:
      ! - f(x) = 0
      else
         shape_function = 0.0d0
      end if

      return
   end function shape_function

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine give_gravity(x0, gravity)
      use amr_parameters, only : dp
      implicit none
      real(dp), intent(in) :: x0
      real(dp), intent(out) :: gravity

       gravity = Kgrav * grav_coef * shape_function(x0)

       return
   end subroutine give_gravity

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine give_cool(x0, dens, pres, cool)
      use amr_parameters, only : dp
      use hydro_parameters, only : gamma
      implicit none
      real(dp), intent(in) :: x0, dens, pres
      real(dp), intent(out) :: cool

      cool = Kheat * heat_coef * dens * shape_function(x0)

      return
   end subroutine give_cool

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

end module buoyant_layer

