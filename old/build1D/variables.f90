module variables

   use amr_parameters

   implicit none

   real(dp), dimension(:,:), allocatable :: uin, gravin, old_uin
   real(dp), dimension(:,:), allocatable :: flux
   real(dp), dimension(:,:), allocatable :: flux_pre_tot
   !real(dp), dimension(:  ), allocatable :: emfx, emfy, emfz
   real(dp), dimension(:), allocatable :: x

   ! Physical model parameters
   real(dp) :: nu  = 0.0d0
   real(dp) :: eta = 0.0d0

   ! Numerical scheme parameters
   logical :: fargo

   real(dp) :: dx
   real(dp) :: dt, courant

   integer, dimension(2) :: bval_type

   ! Curvilinear coordinates parameters
   logical :: cartesian, cylindrical, spherical
   real(dp) :: dimension(:), allocatable :: dv
   real(dp) :: dimension(:), allocatable :: ds

end module variables
