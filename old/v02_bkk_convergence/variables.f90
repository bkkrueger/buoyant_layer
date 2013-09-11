module variables
  use amr_parameters
  implicit none

  real(dp),dimension(:,:,:,:,:  ),allocatable ::uin,gravin,old_uin 
  real(dp),dimension(:,:,:,:,:,:),allocatable ::flux
  real(dp),dimension(:,:,:,:,:  ),allocatable ::flux_pre_tot
  real(dp),dimension(:,:,:,:    ),allocatable ::emfx,emfy,emfz
  real(dp),dimension(:),allocatable ::x,y,z

  real(dp),dimension(:),allocatable ::dens0,velx0,csnd0

  !Physical model parameters
  real(dp)::nu=0.d0,eta=0.d0

  !Numerical scheme parameters
  logical::fargo

  real(dp) ::dx,dy,dz
  real(dp) ::dt,tlim,time,courant

  integer,parameter::ngrid=1
  integer ::restart
  integer, dimension(ndim,2) :: bval_type

  ! Curvilinear coordinates parameters
  logical :: cartesian, cylindrical, spherical
  real(dp),dimension(:,:,:  ),allocatable :: dv
  real(dp),dimension(:,:,:,:),allocatable :: ds

  ! Output frequency parameters
  real(dp) :: dthist,dtdump,dtspec

  ! Collective output parameters
  logical :: use_collective_io
  logical :: inline_io
  logical :: store_ghost
  character(len=80) :: store_sides

end module variables
