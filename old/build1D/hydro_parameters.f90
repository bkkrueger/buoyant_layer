module hydro_parameters

   use amr_parameters

   implicit none

   ! Grid parameters
   real(dp) :: xmin, xmax

   ! Number of independent variables
   integer, parameter :: nvar = 5   ! removing B fields

   ! Size of hydro kernel
   integer, parameter :: nghost = 3
   integer :: iu1, iu2, if1, if2

   ! Hydro solver parameters
   integer :: slope_type = 2
   real(dp) :: gamma = 1.0001d0
   real(dp) :: courant_factor = 0.50d0
   real(dp) :: smallc = 1.0d-10
   real(dp) :: smallr = 1.0d-10
   character(len=10) :: riemann = "hlld"
   integer, parameter :: iroe = 1, illf = 2, ihll = 3, ihlld = 4, &
                         iupwind = 5, iacoustic = 6, ihllf = 7, ihlla = 8
   integer :: iriemann

end module hydro_parameters
