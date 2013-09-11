! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

module perturb

   use amr_parameters, only : dp

   implicit none

   integer :: iseed
   integer :: nx_glob, ny_glob, nz_glob
   double precision, dimension(:,:), allocatable :: fourier_amp, &
      fourier_phasex, fourier_phasey
   real(dp), parameter :: PI = 2.0d0*ASIN(1.0d0)

   contains

   ! ==========================================================================
   ! ==========================================================================

   subroutine compute_fourier_perturbations()

      use amr_parameters, only : dp
      use heating_layer
      use hydro_parameters, only : ymin, ymax

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------

      ! Locals ----------------------------------------------------------------
      integer :: irand, jrand
      real(dp) :: numer, denom
      real(dp) :: rvalue
      real(dp) :: amplitude

      ! =======================================================================
      ! Compute the Fourier components (amplitudes and phases)

      call get_global_size()

      if (pert .eq. 'rho') then
         iseed = pertseed
         allocate(fourier_amp   (nxmin:nxmax, 1:ny_glob/2))
         allocate(fourier_phasex(nxmin:nxmax, 1:ny_glob/2))
         allocate(fourier_phasey(nxmin:nxmax, 1:ny_glob/2))

         do irand = nxmin, nxmax
            do jrand = 1, ny_glob/2

               if ((nymin <= jrand) .and. (jrand <= nymax)) then
                  amplitude = amprand
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
               fourier_phasex(irand,jrand) = 2.d0*PI*rvalue
               call ran2(iseed, rvalue)
               fourier_phasey(irand,jrand) = 2.d0*PI*rvalue
            end do
         end do
      end if

      return
   end subroutine compute_fourier_perturbations

   ! ==========================================================================
   ! ==========================================================================

   subroutine get_global_size()
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
   ! ==========================================================================

   subroutine perturb_base_state(l, i, j, k)

      use amr_parameters, only : dp
      use heating_layer
      use hydro_parameters, only : ymin, ymax
      use variables

      implicit none

      ! =======================================================================
      ! Declare variables

      ! Arguments -------------------------------------------------------------
      integer, intent(in) :: l, i, j, k

      ! Locals ----------------------------------------------------------------
      integer :: irand, jrand
      real(dp) :: rvalue
      real(dp) :: theta1, theta2
      real(dp) :: pert_mag

      ! =======================================================================
      ! Perturb

      if (pert .eq. 'rhogrid') then ! -----------------------------------------
         ! random density perturbations
         call ran2(iseed, rvalue)
         pert_mag = amprand * (rvalue - 0.5)
         uin(l,i,j,k,1) = uin(l,i,j,k,1) + pert_mag
      else if (pert .eq. 'rho') then ! ----------------------------------------
         ! random density perturbations with pre-computed Fourier modes
         do irand = nxmin, nxmax
            do jrand = 1, ny_glob/2
               theta1 = irand*(2.d0*PI*x(i)/(xpmax-xpmin) + &
                        fourier_phasex(irand,jrand))
               theta2 = jrand*(2.d0*PI*y(j)/(ymax-ymin) + &
                        fourier_phasey(irand,jrand))
               pert_mag = fourier_amp(irand,jrand) * cos(theta1) * cos(theta2)
               uin(l,i,j,k,1) = uin(l,i,j,k,1) + pert_mag
            end do
         end do
      else if (pert .eq. 'vgrid') then ! --------------------------------------
         ! random velocity perturbations
         call ran2(iseed, rvalue)
         pert_mag = uin(l,i,j,k,1) * amprand * (rvalue - 0.5)
         uin(l,i,j,k,2) = uin(l,i,j,k,2) + pert_mag         ! rho * velx
         if (ndim .ge. 2) then
            call ran2(iseed, rvalue)
            pert_mag = uin(l,i,j,k,1) * amprand * (rvalue - 0.5)
            uin(l,i,j,k,3) = uin(l,i,j,k,2) + pert_mag      ! rho * vely
            if (ndim .eq. 3) then
               call ran2(iseed, rvalue)
               pert_mag = uin(l,i,j,k,1) * amprand * (rvalue - 0.5)
               uin(l,i,j,k,4) = uin(l,i,j,k,2) + pert_mag   ! rho * velz
            end if
         end if
      end if

      return
   end subroutine perturb_base_state

end module perturb

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine perturb_initial_model(mype, npes)

   use heating_layer, only : xpmin, xpmax
   use hydro_parameters, only : iu1, iu2, ju1, ju2, ku1, ku2
   use perturb
   use variables, only : ngrid, x

   implicit none

   integer, intent(in) :: mype, npes

   integer :: i, j, k, l

   ! Compute the amplitude and phase for the Fourier components
   call compute_fourier_perturbations()

   ! Perturb cells inside the heating layer
   do l = 1, ngrid
      do i = iu1, iu2
         if ((xpmin < x(i)) .and. (x(i) < xpmax)) then
            do j = ju1, ju2
               do k = ku1, ku2
                  call perturb_base_state(l,i,j,k)
               end do
            end do
         end if
      end do
   end do

end subroutine perturb_initial_model

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------
!  numerical recipes random number generator ran2
!    requires input seed value=iseed
!    returns real random number=rvalue
!    Also updates iseed for next call 
!
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
11       continue
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
end

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------
