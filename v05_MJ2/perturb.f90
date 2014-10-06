module perturb

   use amr_parameters, only : dp

   implicit none

   ! ==========================================================================
   ! Parameters

   ! Runtime parameters -------------------------------------------------------

   character(len=20) :: pert_type   ! Type of perturbation
   real(dp) :: pert_amplitude       ! Amplitude of structured perturbation
   real(dp) :: rand_amplitude       ! Amplitude of random perturbation
   integer  :: pertseed             ! Random number generator seed
   real(dp) :: xpmin, xpmax         ! Upper and lower extents of perturbation

   ! Fourier perturbation
   character(len=20) :: pert        ! What style of Fourier perturbation?
   character(len=20) :: spectrum    ! What spectrum in Fourier space?
   real(dp) :: amp_min              ! Minimum amplitude of non-seeded modes
   integer  :: nxmin, nxmax         ! Limits of vertical modes to be seeded
   integer  :: nymin, nymax         ! Limits of horizontal modes to be seeded

   ! Bubble perturbation
   character(len=20) :: bubble_profile ! Cross-sectional profile of bubble
   ! The bubble is an ellipsoid specified by the coordinates of its center and
   ! by its semi-principle axes.
   real(dp) :: bubble_x,  bubble_y,  bubble_z   ! Coordinates of bubble center
   real(dp) :: bubble_dx, bubble_dy, bubble_dz  ! Size of bubble

   ! SASI-like perturbation
   real(dp) :: pert_wavelength      ! Wavelength of perturbation
   integer :: n_horiz               ! Number of horizontal wavelengths
   real(dp) :: L_damp               ! Damping length
   logical :: BC_drive_pert         ! Driving from the upper BC or no

   ! For tracking bubble information in history.f90
   real(dp) :: x_mass_old, x_advc_old, time_old
   real(dp) :: dv_grav_old, dv_buoy_old, dv_drag_old

   ! Other variables ----------------------------------------------------------

   integer :: iseed  ! initialized to pertseed, changes with each call to ran2
   integer :: nx_glob, ny_glob, nz_glob

   ! Fourier perturbation derived parameters
   real(dp), dimension(:,:), allocatable :: fourier_amp, &
                                            fourier_phasex, fourier_phasey

   ! SASI-like perturbation derived parameters
   real(dp) :: osc_omega, osc_slope, inv_L_damp
   real(dp) :: L_horiz

   real(dp), parameter :: TWOPI = 4.0d0*asin(1.0d0)

contains

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine apply_perturbation(mype)

      implicit none

      integer, intent(in) :: mype

      ! Seed the random number generator for the perturbations
      ! -- The user-defined pertseed allows the same simulation to be run again
      !    in a controlled manner.
      ! -- Adding mype eliminates repetition across processors (e.g., without
      !    this if I have the SASI-like perturbation and nyslice == 4, then my
      !    simulation has four side-by-side copies of the same thing because the
      !    random number generator sets the same random perturbations on each
      !    slice instead of being completely random across the domain).
      iseed = pertseed + mype

      if (trim(pert_type) .eq. 'bubble') then
         call perturb_bubble()
      else if (trim(pert_type) .eq. 'fourier') then
         call perturb_fourier(mype)
      else if (trim(pert_type) .eq. 'sasi') then
         call perturb_sasi(mype)
      end if

   end subroutine apply_perturbation

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine perturb_bubble()
      use amr_parameters, only : dp, ndim
      use hydro_parameters, only : gamma, iu1, iu2, ju1, ju2, ku1, ku2
      use variables, only : ngrid, uin, x, y, z
      implicit none
      integer :: i, j, k, l
      real(dp) :: xx, yy, zz, dd, bubble, rand
      real(dp) :: dens_old, dens_new
      real(dp) :: velx, vely, velz, v_sq, pres

      do l = 1, ngrid
         do i = iu1, iu2
            xx = abs((x(i) - bubble_x) / bubble_dx)
            do j = ju1, ju2
               if (ndim .ge. 2) then
                  yy = abs((y(j) - bubble_y) / bubble_dy)
               else
                  yy = 0.0d0
               endif
               do k = ku1, ku2
                  if (ndim .ge. 3) then
                     zz = abs((z(k) - bubble_z) / bubble_dz)
                  else
                     zz = 0.0d0
                  end if

                  dd = xx**2 + yy**2 + zz**2

                  ! Choose the shape of the bubble's profile
                  bubble = 0.0d0
                  ! no bubble
                  if (bubble_profile .eq. 'none') then
                     bubble = 0.0d0
                  ! parabolic cross-section
                  else if (bubble_profile .eq. 'parabola') then
                     if (sqrt(dd) < 1.0d0) then
                        bubble = 1.0d0 - dd
                     end if
                  ! gaussian cross-section
                  else if (bubble_profile .eq. 'gauss') then
                     bubble = exp(-1.0d0*dd)
                  ! top hat (square) cross-section
                  else if (bubble_profile .eq. 'tophat') then
                     if (sqrt(dd) <= 1.0d0) then
                        bubble = 1.0d0
                     end if
                  ! unknown bubble type: turn off the bubble
                  else
                     bubble = 0.0d0
                  end if

                  ! Scale the bubble by the user-defined amplitude
                  bubble = bubble * pert_amplitude
                  ! small randomization on top of wave perturbation
                  call ran2(iseed, rand)
                  rand = rand_amplitude * (2.0d0 * rand - 1.0d0)
                  bubble = bubble * (1.0d0 + rand)

                  ! Change density but maintain velocity and pressure
                  dens_old = uin(l,i,j,k,1)
                  velx = uin(l,i,j,k,2) / dens_old
                  vely = uin(l,i,j,k,3) / dens_old
                  velz = uin(l,i,j,k,4) / dens_old
                  v_sq = velx**2 + vely**2 + velz**2
                  pres = (uin(l,i,j,k,5) - 0.5d0 * dens_old * v_sq) * &
                         (gamma - 1.0d0)

                  dens_new = dens_old * (1.0d0 + bubble)

                  uin(l,i,j,k,1) = dens_new
                  uin(l,i,j,k,2) = dens_new * velx
                  uin(l,i,j,k,3) = dens_new * vely
                  uin(l,i,j,k,4) = dens_new * velz
                  uin(l,i,j,k,5) = pres/(gamma-1.0d0) + 0.5d0*dens_new*v_sq

               end do
            end do
         end do
      end do

      return
   end subroutine perturb_bubble

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine perturb_fourier(mype)

      implicit none

      integer, intent(in) :: mype

      call perturb_fourier_setup(mype)
      call perturb_fourier_apply()

      return
   end subroutine perturb_fourier

   ! ==========================================================================

   subroutine perturb_fourier_setup(mype)

      use amr_parameters, only : dp
      use hydro_parameters, only : ymin, ymax

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      integer, intent(in) :: mype
      integer :: irand, jrand
      integer :: offset, idx
      integer :: ny2
      real(dp) :: numer, denom, rvalue, amplitude
#     if WITHMPI==1
      real(dp), dimension(:), allocatable :: mpi_buffer
      integer :: ierr
#     endif

      call get_global_size()
      ny2 = ny_glob / 2

      if (pert .eq. 'rho') then
         allocate(fourier_amp   (nxmin:nxmax, 1:ny2))
         allocate(fourier_phasex(nxmin:nxmax, 1:ny2))
         allocate(fourier_phasey(nxmin:nxmax, 1:ny2))

         do irand = nxmin, nxmax
            do jrand = 1, ny2

               if ((nymin <= jrand) .and. (jrand <= nymax)) then
                  amplitude = pert_amplitude
               else
                  amplitude = amp_min
               end if

               if (spectrum .eq. 'flat') then
                  fourier_amp(irand,jrand) = amplitude
               else
                  numer = (nxmin/(xpmax-xpmin))**2 + (nymin/(ymax-ymin))**2
                  denom = (irand/(xpmax-xpmin))**2 + (jrand/(ymax-ymin))**2
                  fourier_amp(irand,jrand) = amplitude * sqrt(numer/denom)
               end if

               call ran2(iseed, rvalue)
               fourier_phasex(irand,jrand) = TWOPI * rvalue
               call ran2(iseed, rvalue)
               fourier_phasey(irand,jrand) = TWOPI * rvalue

            end do
         end do
      end if

#     if WITHMPI==1
      allocate(mpi_buffer(2*(nxmax-nxmin+1)*ny2))

      if (mype == 0) then
         do irand = nxmin, nxmax
            offset = (irand - nxmin) * ny2
            do jrand = 1, ny2
               idx = offset + jrand
               mpi_buffer(idx) = fourier_phasex(irand, jrand)
            end do
         end do
         do irand = nxmin, nxmax
            offset = (nxmax-nxmin+1)*ny2 + (irand-nxmin) * ny2
            do jrand = 1, ny2
               idx = offset + jrand
               mpi_buffer(idx) = fourier_phasey(irand, jrand)
            end do
         end do
      end if

      call MPI_Bcast(mpi_buffer, 2*(nxmax-nxmin+1)*ny2, &
                     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      if (mype /= 0) then
         do irand = nxmin, nxmax
            offset = (irand - nxmin) * ny2
            do jrand = 1, ny2
               idx = offset + jrand
               fourier_phasex(irand, jrand) = mpi_buffer(idx)
            end do
         end do
         do irand = nxmin, nxmax
            offset = (nxmax-nxmin+1)*ny2 + (irand-nxmin) * ny2
            do jrand = 1, ny2
               idx = offset + jrand
               fourier_phasey(irand,jrand) = mpi_buffer(idx)
            end do
         end do
      end if

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      deallocate(mpi_buffer)
#     endif

      return
   end subroutine perturb_fourier_setup

   ! ==========================================================================

   subroutine get_global_size()
      ! Get n[xyz]_glob; requires reading the full scheme_params block of the
      ! input file, so create dummy values for the parameters I don't need here
      use amr_parameters, only : dp, ndim
      implicit none
      character(len=10) :: riemann, riemann2d
      integer :: slope_type
      real(dp), dimension(ndim) :: bval_in, bval_out
      real(dp) :: courant
      logical :: fargo
      namelist /scheme_params/ nx_glob, ny_glob, nz_glob, riemann, riemann2d, &
                               slope_type, bval_in, bval_out, courant, fargo
      open(unit=53, file='input', status='old')
      read(53, scheme_params)
      close(53)
      return
   end subroutine get_global_size

   ! ==========================================================================

   subroutine perturb_fourier_apply()

      use amr_parameters, only : dp, ndim
      use hydro_parameters, only : ymin, ymax, iu1, iu2, ju1, ju2, ku1, ku2
      use variables, only : ngrid, uin, x, y, z

      implicit none

      integer :: l, i, j, k
      integer :: irand, jrand
      real(dp) :: rvalue
      real(dp) :: theta1, theta2
      real(dp) :: pert_mag

      do l = 1, ngrid
         do i = iu1, iu2
            if ((xpmin < x(i)) .and. (x(i) < xpmax)) then
               do j = ju1, ju2
                  do k = ku1, ku2

                     if (pert .eq. 'rhogrid') then ! --------------------------
                        ! random density perturbations
                        call ran2(iseed, rvalue)
                        pert_mag = pert_amplitude * (rvalue - 0.5d0)
                        uin(l,i,j,k,1) = uin(l,i,j,k,1) + pert_mag

                     else if (pert .eq. 'rho') then ! -------------------------
                        ! random density perturbations with Fourier modes
                        do irand = nxmin, nxmax
                           do jrand = 1, ny_glob/2
                              theta1 = irand*(TWOPI*x(i)/(xpmax-xpmin) + &
                                              fourier_phasex(irand,jrand))
                              theta2 = jrand*(TWOPI*y(j)/( ymax- ymin) + &
                                              fourier_phasey(irand,jrand))
                              pert_mag = fourier_amp(irand,jrand) * &
                                          cos(theta1) * cos(theta2)
                              uin(l,i,j,k,1) = uin(l,i,j,k,1) + pert_mag
                           end do
                        end do

                     else if (pert .eq. 'vgrid') then ! -----------------------
                        ! random velocity perturbations
                        call ran2(iseed, rvalue)
                        pert_mag = uin(l,i,j,k,1) * &
                                   pert_amplitude * (rvalue-0.5d0)
                        uin(l,i,j,k,2) = uin(l,i,j,k,2) + pert_mag
                        if (ndim .ge. 2) then
                           call ran2(iseed, rvalue)
                           uin(l,i,j,k,3) = uin(l,i,j,k,1) * &
                                            pert_amplitude * (rvalue-0.5d0)
                           if (ndim .ge. 3) then
                              uin(l,i,j,k,4) = uin(l,i,j,k,1) * &
                                               pert_amplitude * (rvalue-0.5d0)
                           end if
                        end if

                     end if ! -------------------------------------------------

                  end do
               end do
            end if
         end do
      end do

      return
   end subroutine perturb_fourier_apply

   ! ==========================================================================

   subroutine ran2(iseed,rvalue)

      integer iseed
      real rvalue
      integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      parameter (IM1=2147483563, IM2=2147483399, AM=1./IM1, IMM1=IM1-1, &
                 IA1=40014, IA2=40692, IQ1=53668, IQ2=52774, IR1=12211, &
                 IR2=3791, NTAB=32, NDIV=1+IMM1/NTAB, EPS=1.2e-7, RNMX=1.-EPS)
      integer idum2,jj,kk,iv(NTAB),iy
      data idum2/123456789/, iv/NTAB*0/, iy/0/

      idum=iseed
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 jj=NTAB+8,1,-1
            kk=idum/IQ1
            idum=IA1*(idum-kk*IQ1)-kk*IR1
            if (idum.lt.0) idum=idum+IM1
            if (jj.le.NTAB) iv(jj)=idum
11          continue
         iy=iv(1)
      endif
      kk=idum/IQ1
      idum=IA1*(idum-kk*IQ1)-kk*IR1
      if (idum.lt.0) idum=idum+IM1
      kk=idum2/IQ2
      idum2=IA2*(idum2-kk*IQ2)-kk*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      jj=1+iy/NDIV
      iy=iv(jj)-idum2
      iv(jj)=idum
      if(iy.lt.1)iy=iy+IMM1
      rvalue=min(AM*iy,RNMX)
      iseed=idum
      return
   end subroutine ran2

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine perturb_sasi(mype)
      use amr_parameters, only : dp, ndim, nx
      use buoyant_layer, only : layer_limit
      use hydro_parameters, only : gamma, iu1, iu2, ju1, ju2, ku1, ku2, &
                                   xmax, ymin, ymax
      use steady_state, only : mach_up, csnd_up, ss0, SS_DENS, SS_MOMX, SS_ENER
      use variables, only : ngrid, uin, x, y, z, dx

      implicit none
      integer, intent(in) :: mype

      integer :: i, j, k, l
      real(dp) :: L_horiz, lambda_vert
      real(dp) :: wave, amplitude, rand, perturb, damp
      real(dp) :: dens_old, dens_new
      real(dp) :: velx, vely, velz, v_sq, csnd, Eint

      ! Convert amplitudes out of logscale -------------------------------------
      pert_amplitude = 10.0d0**pert_amplitude
      rand_amplitude = 10.0d0**rand_amplitude

      ! Calculate perturbation geometry parameters -----------------------------

      L_horiz = ymax - ymin

      ! Verify that n <= n_max
      if (n_horiz > floor(L_horiz / pert_wavelength)) then
         write(*,*) "ERROR : Cannot fit ", n_horiz, " wavelengths."
         write(*,*) "        n_max = ", floor(L_horiz / pert_wavelength)
         stop
      end if

      ! Compute vertical wavelength
      lambda_vert = 1.0 / sqrt(pert_wavelength**(-2) - (n_horiz / L_horiz)**2)

      ! Compute vertical spatial frequency
      osc_omega = TWOPI / lambda_vert

      ! Compute slope
      osc_slope = lambda_vert * n_horiz / L_horiz

      ! Apply perturbation -----------------------------------------------------
      inv_L_damp = 1.0d0 / L_damp

      do l = 1, ngrid
         do i = iu1, iu2

            ! Compute the damping envelope
            ! -- between xpmin and xpmax, no damping
            ! -- if below xpmin, damp
            ! -- if we are not continuously driving the perturbation from
            !    the upper BC, damp above xpmax
            ! -- if we are continuously driving the perturbation from the
            !    upper BC, no damping above xpmax
            damp = 1.0d0
            if (x(i) < xpmin) then
               damp = damp * exp(-1.0d0 * ((x(i) - xpmin) * inv_L_damp)**2)
            end if
            if ((x(i) > xpmax) .and. (.not. BC_drive_pert)) then
               damp = damp * exp(-1.0d0 * ((x(i) - xpmax) * inv_L_damp)**2)
            end if

            do j = ju1, ju2
               do k = ku1, ku2

                  ! Construct the perturbation --------------------------------

                  ! Compute the wave shape:
                  ! -- wave = sin[w * (x - x0 - velx_up*t)]
                  ! -- x0 = osc_slope * y
                  ! -- t = 0 for initial conditions
                  wave = sin(osc_omega * (x(i) - osc_slope * y(j)))

                  ! small randomization on top of wave perturbation
                  call ran2(iseed, rand)
                  rand = rand_amplitude * (2.0d0 * rand - 1.0d0)
                  amplitude = pert_amplitude * (1.0d0 + rand)

                  ! Assemble the full perturbation
                  perturb = amplitude * damp * wave

                  ! Apply the perturbation ------------------------------------
                  ! -- Change density but maintain velocity and pressure

                  ! Get unperturbed primitive variables
                  dens_old = uin(l,i,j,k,1)
                  velx = uin(l,i,j,k,2) / dens_old
                  vely = uin(l,i,j,k,3) / dens_old
                  velz = uin(l,i,j,k,4) / dens_old
                  ! note: Constant velocity, so constant sum(v_i^2)
                  v_sq = velx**2 + vely**2 + velz**2
                  ! note: We use Eint, which is the density multiplied by
                  !       the specific internal energy, because it is more
                  !       convenient.  With an ideal-gas EoS,
                  !       pres = Eint * (gamma - 1)
                  !       so that a constant Eint is equivalent to a
                  !       constant pressure.  Using Eint saves us from
                  !       multiplying by (gamma - 1) to get pressure, only to
                  !       divide by (gamma - 1) again in a few lines.
                  Eint = uin(l,i,j,k,5) - 0.5d0 * dens_old * v_sq

                  ! Perturb density
                  dens_new = dens_old * (1.0d0 + perturb)

                  ! Construct perturbed values of conserved variables
                  uin(l,i,j,k,1) = dens_new
                  uin(l,i,j,k,2) = dens_new * velx
                  uin(l,i,j,k,3) = dens_new * vely
                  uin(l,i,j,k,4) = dens_new * velz
                  uin(l,i,j,k,5) = Eint + 0.5d0*dens_new*v_sq

               end do
            end do
         end do
      end do

      return
   end subroutine perturb_sasi

end module perturb

