module amr_parameters

   implicit none

   ! Define real types
   integer, parameter :: sp = kind(1.0e0)
   integer, parameter :: dp = kind(1.0d0)

   ! Mesh parameters
   integer :: geom = 1  ! 1: Cartesian; 2: cylindrical; 3: spherical
   integer :: nx   = 1  ! Number of cells in each dimension

end module amr_parameters
