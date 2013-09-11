! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

module heating_layer
   use amr_parameters, only : dp
   implicit none
   real(dp)          :: Mup, Kheat, Kgrav
   character(len=20) :: layer_shape
   real(dp)          :: layer_limit
   real(dp)          :: amprand
   real(dp)          :: amp_min
   character(len=20) :: pert,spectrum
   real(dp)          :: xpmin, xpmax
   integer           :: nxmin, nxmax
   integer           :: nymin, nymax
   integer           :: pertseed
   logical           :: do_refine
end module heating_layer

