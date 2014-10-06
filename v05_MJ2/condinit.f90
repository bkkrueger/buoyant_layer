! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine condinit(mype)
   use amr_parameters, only : dp
   use variables, only : restart
   use hydro_parameters, only : gamma, iu1, iu2
   use buoyant_layer
   use perturb
   use steady_state
   implicit none 

   ! ==========================================================================
   ! Declare variables --------------------------------------------------------

   ! Arguments ----------------------------------------------------------------
   integer :: mype

   ! Locals -------------------------------------------------------------------

   integer :: i, j, k, l                  ! loop indices
   real(dp) :: v, c, dens, ener
   integer, parameter :: LUN = 16

   !  RK looping
   integer :: irk
   integer :: Nskip
   integer, parameter :: NRK_2 = 5000
   integer, parameter :: NRK = NRK_2 * 2  ! because NRK must be even
   integer :: ncount                      ! RK looping
   double precision :: dx0, x0
   double precision :: y0(2)
   character(len=17) :: fname
   character(len=40) :: f_str, temp_str

   ! Name lists for problem-specific parameters
   namelist /buoyant_params/ Kheat, Kgrav, layer_shape, layer_limit
   namelist /physics_params/ mach_up, gamma
   namelist /perturb_params/ pert_type, rand_amplitude, pertseed,             &
                             pert, spectrum, pert_amplitude, amp_min,         &
                             xpmin, xpmax, nxmin, nxmax, nymin, nymax,        &
                             bubble_profile,                                  &
                             bubble_dx, bubble_dy, bubble_dz,                 &
                             bubble_x, bubble_y, bubble_z,                    &
                             pert_wavelength, n_horiz, L_damp, BC_drive_pert

   ! ==========================================================================
   ! Set the parameter values

   ! Default values -----------------------------------------------------------

   ! Background state
   mach_up        = 1.0d-1
   gamma          = 4.0d0 / 3.0d0

   ! Buoyant layer
   Kheat          = 1.0d-2
   Kgrav          = 1.0d0
   layer_shape    = 'trapezoid'
   layer_limit    = 1.0d0  ! layer extends from -layer_limit to +layer_limit

   ! Perturbations
   pert_type      = 'none'    ! Type of perturbation
   pert_amplitude = 0.0d0     ! Amplitude of perturbation
   rand_amplitude = 0.0d0     ! Amplitude of random noise
   pertseed       = 2         ! Random number seed
   xpmin          = 0.0d0     ! Lower limit of perturbed region
   xpmax          = 0.0d0     ! Upper limit of perturbed region

   ! Fourier
   pert           = 'rho'
   spectrum       = 'flat'
   amp_min        = 0.0d0
   nxmin          = 0
   nxmax          = 0
   nymin          = 0
   nymax          = 0

   ! Bubble
   bubble_profile = 'tophat'
   bubble_x       = 0.0d0
   bubble_y       = 0.0d0
   bubble_z       = 0.0d0
   bubble_dx      = layer_limit
   bubble_dy      = layer_limit
   bubble_dz      = layer_limit
  
   ! SASI-like (large-scale, smooth)
   pert_wavelength = 1.0d0    ! Wavelength of perturbation
   n_horiz         = 1        ! Number of oscillations to fit horizontally
   L_damp          = 1.0d100  ! Length over which to damp perturbation
   BC_drive_pert   = .false.  ! Continuously drive perturbation from upper BC?

   ! Read the parameters ------------------------------------------------------

   ! Read from the input file
   open(unit=LUN,file='input' ,status='old')
   read(LUN,physics_params)
   read(LUN,buoyant_params)
   read(LUN,perturb_params)
   close(LUN)

   ! Compute a couple support parameters
   grav_coef = -1.0d0 * csnd_up**2 / layer_limit
   ! Note that heat_coef also should have a dens_up factor to make everything
   ! properly dimensioned, but the dens_up in heat_coef cancels with the
   ! denominator of the unitless (dens/dens_up) term.  Therefore instead of
   ! having the units of the heat function, heat_coef has the units of the heat
   ! function divided by density.
   heat_coef = mach_up * csnd_up * abs(grav_coef) / gamma

   ! Print the resulting values
   if (mype == 0) then
      f_str  = "(3x,a15,x,'=',x,a20)"
      write(*,*) "physics parameters"
      write(temp_str,"(es13.6)") gamma
      write(*,f_str) "gamma", trim(temp_str)
      write(temp_str,"(es13.6)") mach_up
      write(*,f_str) "mach_up", trim(temp_str)
      write(temp_str,"(es13.6)") csnd_up
      write(*,f_str) "csnd_up", trim(temp_str)
      write(temp_str,"(es13.6)") dens_up
      write(*,f_str) "dens_up", trim(temp_str)
      write(*,*) ""
      write(*,*) "buoyant layer parameters"
      write(temp_str,"(es13.6)") Kheat
      write(*,f_str) "Kheat", trim(temp_str)
      write(temp_str,"(es13.6)") Kgrav
      write(*,f_str) "Kgrav", trim(temp_str)
      write(*,f_str) "layer_shape", trim(layer_shape)
      write(temp_str,"(es13.6)") layer_limit
      write(*,f_str) "layer_limit", trim(temp_str)
      write(*,*) ""
      if (trim(pert_type) == 'none') then
         write(*,*) "perturbation parameters: no perturbation"
      else if (trim(pert_type) == 'fourier') then
         write(*,*) "perturbation parameters: Fourier perturbation"
         write(*,f_str) "pert", trim(pert)
         write(*,f_str) "spectrum", trim(spectrum)
         write(temp_str,"(es13.6)") pert_amplitude
         write(*,f_str) "pert_amplitude",   trim(temp_str)
         write(temp_str,"(es13.6)") amp_min
         write(*,f_str) "amp_min",   trim(temp_str)
         write(temp_str,"(es13.6)") xpmin
         write(*,f_str) "xpmin",   trim(temp_str)
         write(temp_str,"(es13.6)") xpmax
         write(*,f_str) "xpmax",   trim(temp_str)
         write(temp_str,"(i13)") nxmin
         write(*,f_str) "nxmin",   trim(temp_str)
         write(temp_str,"(i13)") nxmax
         write(*,f_str) "nxmax",   trim(temp_str)
         write(temp_str,"(i13)") nymin
         write(*,f_str) "nymin",   trim(temp_str)
         write(temp_str,"(i13)") nymax
         write(*,f_str) "nymax",   trim(temp_str)
         write(temp_str,"(i13)") pertseed
         write(*,f_str) "pertseed",   trim(temp_str)
      else if (trim(pert_type) == 'bubble') then
         write(*,*) "perturbation parameters: bubble perturbation"
         write(*,f_str) "profile", bubble_profile
         write(temp_str,"(es13.6)") pert_amplitude
         write(*,f_str) "pert_amplitude", trim(temp_str)
         write(temp_str,"(es13.6)") bubble_x
         write(*,f_str) "center (x)", trim(temp_str)
         write(temp_str,"(es13.6)") bubble_y
         write(*,f_str) "center (y)", trim(temp_str)
         write(temp_str,"(es13.6)") bubble_z
         write(*,f_str) "center (z)", trim(temp_str)
         write(temp_str,"(es13.6)") bubble_dx
         write(*,f_str) "radius (x)", trim(temp_str)
         write(temp_str,"(es13.6)") bubble_dy
         write(*,f_str) "radius (y)", trim(temp_str)
         write(temp_str,"(es13.6)") bubble_dz
         write(*,f_str) "radius (z)", trim(temp_str)
      else if (trim(pert_type) == 'sasi') then
         write(*,*) "perturbation parameters: SASI-like perturbation"
         write(temp_str,"(es13.6)") pert_amplitude
         write(*,f_str) "pert_amplitude", trim(temp_str)
         write(temp_str,"(es13.6)") rand_amplitude
         write(*,f_str) "rand_amplitude", trim(temp_str)
         write(temp_str,"(i13)") pertseed
         write(*,f_str) "pertseed",   trim(temp_str)
         write(temp_str,"(es13.6)") xpmin
         write(*,f_str) "xpmin",   trim(temp_str)
         write(temp_str,"(es13.6)") xpmax
         write(*,f_str) "xpmax",   trim(temp_str)
         write(temp_str,"(es13.6)") pert_wavelength
         write(*,f_str) "pert_wavelength", trim(temp_str)
         write(temp_str,"(i13)") n_horiz
         write(*,f_str) "n_horiz", trim(temp_str)
         write(temp_str,"(es13.6)") L_damp
         write(*,f_str) "L_damp", trim(temp_str)
         write(temp_str,"(l13)") BC_drive_pert
         write(*,f_str) "BC_drive_pert", trim(temp_str)
      else
         write(*,*) "ERROR : bad perturbation type selected"
         stop
      end if
      write(*,*) ""
   end if

   ! Something gravity related?
   call update_gravin

   ! ==========================================================================
   ! Set the initial conditions -----------------------------------------------

   allocate(ss0(iu1:iu2,SS_NVAR))
   if (restart <= 0) then
      call ss_build()
      !call ss_correct() ! TODO
      call ss_write_to_file()
      call apply_perturbation(mype)
   end if

end subroutine condinit

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

