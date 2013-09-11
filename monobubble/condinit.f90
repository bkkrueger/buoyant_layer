! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

module heating_layer
   use hydro_parameters
   implicit none
   real(dp)     :: Mup, Kheat, Kgrav
   character*20 :: layer_shape
   real(dp)     :: layer_limit
   real(dp)     :: bubble_x,  bubble_y,  bubble_z
   real(dp)     :: bubble_dx, bubble_dy, bubble_dz
   real(dp)     :: bubble_amplitude
   character*20 :: bubble_profile
end module heating_layer

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine add_bubble()
   use amr_parameters, only : dp
   use heating_layer
   use variables
   implicit none
   integer :: i, j, k, l
   real(dp) :: xx, yy, zz, dd, bubble
   real(dp) :: rho_old, rho_new, mom_sq, delta_E

   do l = 1, ngrid
      do i = iu1, iu2
         xx = abs(2.0d0 * (x(i) - bubble_x) / bubble_dx)
         do j = ju1, ju2
            if (ndim .gt. 1) then
               yy = abs(2.0d0 * (y(j) - bubble_y) / bubble_dy)
            else
               yy = 0.0
            end if
            do k = ku1, ku2
               if (ndim .gt. 2) then
                  zz = abs(2.0d0 * (z(k) - bubble_z) / bubble_dz)
               else
                  zz = 0.0d0
               end if

               dd = xx**2 + yy**2 + zz**2

               bubble = 0.0d0
               if (bubble_profile .eq. 'parabola') then
                  bubble = exp(-1.0d0*dd)
               else if (bubble_profile .eq. 'gauss') then
                  if (sqrt(dd) <= 1.0d0) then
                     bubble = 1.0d0 - dd
                  end if
               else ! 'tophat' is the catch-all in case of typo
                  if (sqrt(dd) <= 1.0d0) then
                      bubble = 1.0d0
                  end if
               end if
               rho_old = uin(l,i,j,k,1)
               rho_new = uin(l,i,j,k,1) + bubble_amplitude * bubble
               uin(l,i,j,k,1) = rho_new
               mom_sq  = uin(l,i,j,k,2)**2 + &
                         uin(l,i,j,k,3)**2 + &
                         uin(l,i,j,k,4)**2
               delta_E = -0.5 * mom_sq * (1.0 / rho_old - 1.0 / rho_new)
               uin(l,i,j,k,5) = uin(l,i,j,k,5) + delta_E

            end do
         end do
      end do
   end do

end subroutine add_bubble

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine condinit(mype)
   use variables
   use hydro_parameters
   use heating_layer
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
   double precision :: stable(iu1:iu2,nvar)
   character(len=17) :: fname
   character(len=20) :: temp_str
   character(len=80) :: f_str

   ! Physics Parameters namelist
   namelist /physics_params/ Mup, Kheat, Kgrav, gamma, layer_shape, &
      layer_limit, bubble_x, bubble_y, bubble_z, &
      bubble_dx, bubble_dy, bubble_dz, bubble_profile, bubble_amplitude

   ! ==========================================================================
   ! Fill the initial conditions

   ! Default values -----------------------------------------------------------

   ! Background state
   ! c_up = 1
   ! rho_up = 1
   Mup   = 1.0d-1
   Kheat = 1.0d-2
   Kgrav = 1.0d0
   gamma = 4.0 / 3.0

   ! Buoyant layer
   layer_shape = 'trapezoid'
   layer_limit = 1.0d0  ! layer extends from -layer_limit to +layer_limit

   ! Bubble
   bubble_x         = 0.0d0
   bubble_y         = 0.0d0
   bubble_z         = 0.0d0
   bubble_dx        = 0.5 * layer_limit
   bubble_dy        = 0.5 * layer_limit
   bubble_dz        = 0.5 * layer_limit
   bubble_amplitude = 0.0d0
   bubble_profile   = 'tophat'
  
   ! Read the physics parameters ----------------------------------------------
   open(unit=LUN,file='input' ,status='old')
   read(LUN,physics_params)
   close(LUN)
   if (mype == 0) then
      write(*,*) "physics parameters"
      f_str  = "(3x,a11,x,'=',x,a20)"
      write(temp_str,"(es13.6)") Mup
      write(*,f_str) "Mup",   trim(temp_str)
      write(temp_str,"(es13.6)") Kheat
      write(*,f_str) "Kheat", trim(temp_str)
      write(temp_str,"(es13.6)") Kgrav
      write(*,f_str) "Kgrav", trim(temp_str)
      write(temp_str,"(es13.6)") gamma
      write(*,f_str) "gamma", trim(temp_str)
      write(*,f_str) "layer_shape", trim(layer_shape)
      write(temp_str,"(es13.6)") layer_limit
      write(*,f_str) "layer_limit", trim(temp_str)
      write(*,*) "bubble parameters"
      write(*,f_str) "profile", trim(bubble_profile)
      write(temp_str,"(es13.6)") bubble_amplitude
      write(*,f_str) "amplitude", trim(temp_str)
      if (ndim .gt. 2) then
         f_str = "(3x,a11,x,'=',x,'(',es13.6,',',es13.6,',',es13.6,')')"
         write(*,f_str) "position", bubble_x, bubble_y, bubble_z
         write(*,f_str) "size", bubble_dx, bubble_dy, bubble_dz
      else if (ndim .gt. 1) then
         f_str = "(3x,a11,x,'=',x,'(',es13.6,',',es13.6,')')"
         write(*,f_str) "position", bubble_x, bubble_y
         write(*,f_str) "size", bubble_dx, bubble_dy
      else
         f_str = "(3x,a11,x,'=',x,es13.6)"
         write(*,f_str) "position", bubble_x
         write(*,f_str) "size", bubble_dx
      end if
   end if

   ! Something gravity related?
   call update_gravin

   ! Fill the unperturbed state -----------------------------------------------
   ! uin   (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)
   ! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)

   ! Integrate from xmax to x(iu1) (the lowest x-value in this slice)
   x0    = xmax
   y0(1) = -Mup   ! v
   y0(2) = 1.0d0  ! c
   ! NRK = number of RK steps between x(i) and x(i-1)
   dx0   = -1 * dx / NRK

   ! Integrate down to x(iu2)
   if (xmax > x(iu2)) then
      ! If xmax > x(iu2), i.e. we are in a slice that is not along the
      ! upper x-boundary of the domain, then we need to integrate down from
      ! xmax to x(iu2).  In this case (xmax-x(iu2))/dx is a half-integer
      ! (because the x(i) values are cell centers and xmax is the upper edge of
      ! the uppermost non-guard cells in the x-direction).  Therefore we can
      ! integrate NRK/2 steps to reach a cell center, then integrate NRK*Nskip
      ! times, where Nskip is the number of cells within the domain above
      ! x(iu2).  We can compute Nskip as floor((xmax - x(iu2))/dx), because
      ! (xmax - x(iu2))/dx = Nskip+0.5.
      Nskip = floor((xmax-x(iu2))/dx)
      do i = 1, Nskip*NRK + NRK_2
         call rkutta0(x0, y0, dx0, ncount)
         x0 = x0 + dx0
      end do
   else
      ! If x(iu2) > xmax and we assume that the heating layer is fully
      ! contained within the domain then x(iu2) can safely be set to the
      ! upstream values given by the user.
      x0 = x(iu2)
   end if

   ! Translate to easier-to-use variables
   v    = y0(1)
   c    = y0(2)
   dens = -Mup / v
   ener = dens * c**2 / (gamma*(gamma-1.0d0)) + 0.5d0 * dens * v**2

   ! Fill all cells at x = x(iu2)
   uin(:,iu2,:,:,:) = 0.0d0
   uin(:,iu2,:,:,1) = dens
   uin(:,iu2,:,:,2) = -Mup   ! rho*v not just v
   uin(:,iu2,:,:,5) = ener
   stable(iu2,:) = 0.0d0
   stable(iu2,1) = dens
   stable(iu2,2) = -Mup
   stable(iu2,5) = ener

   ! Integrate across the slice
   do i = iu2 - 1, iu1, -1

      ! Integrate from x(i) down to x(i-1)
      do irk = 1, NRK
         call rkutta0(x0, y0, dx0, ncount)
         x0 = x0 + dx0
      end do

      ! Translate to easier-to-use variables
      v    = y0(1)
      c    = y0(2)
      dens = -Mup / v
      ener = dens * c**2 / (gamma*(gamma-1.0d0)) + 0.5d0 * dens * v**2

      ! Fill all cells at x = x(i)
      uin(:,i,:,:,:) = 0.0d0
      uin(:,i,:,:,1) = dens
      uin(:,i,:,:,2) = -Mup   ! rho*v not just v
      uin(:,i,:,:,5) = ener
      stable(i,:) = 0.0d0
      stable(i,1) = dens
      stable(i,2) = -Mup
      stable(i,5) = ener

   end do

   ! Write the stable state ---------------------------------------------------

   write(fname,"('stable_',I0.6,'.dat')") mype  ! fname has length 17
   open(unit=LUN, file=fname, status='unknown')
   write(LUN,"('#',x,50(3x,a30))") "radius", "density", "x-momentum", &
      "y-momentum", "z-momentum", "total_energy_dens", "magnetic_field_x", &
      "magnetic_field_y", "magnetic_field_z"
   do i = iu1, iu2
      write(LUN,"(2x,50(3x,es30.23))") x(i), stable(i,:)
   end do
   close(LUN)

   ! Perturb the base state ---------------------------------------------------

   call add_bubble()

end subroutine condinit

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

double precision function shape_function(x0)
   use heating_layer, only : layer_shape, layer_limit
   implicit none
   double precision :: x0

   if (layer_shape .eq. 'trapezoid') then
      if (abs(x0) > layer_limit) then
         shape_function = 0.0d0
      else if (abs(x0) <= 0.5d0*layer_limit) then
         shape_function = 1.0d0
      else
         shape_function = 2.0d0 * (1.0d0-abs(x0)/layer_limit)
      end if
   else if (layer_shape .eq. 'square') then
      if (abs(x0) > layer_limit) then
         shape_function = 0.0d0
      else
         shape_function = 1.0d0
      end if
   else
      shape_function = 0.0d0
   end if

   return
end function shape_function

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine give_gravity(x0, gravity)
    use variables
    use hydro_parameters
    use heating_layer
    implicit none
    double precision, intent(in) :: x0
    double precision, Intent(out) :: gravity
    double precision :: shape_function

    gravity = -1.0d0 * Kgrav * shape_function(x0)

    return
end subroutine give_gravity

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine give_cool(x0,rho,P,cool)
   use variables
   use hydro_parameters
   use heating_layer
   implicit none
   double precision, intent(in) :: x0, rho, P
   double precision, Intent(out) :: cool
   double precision :: shape_function

   cool = (Kheat * rho * Mup / gamma) * shape_function(x0)

   return
end subroutine give_cool

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine derive0(x0,y0,dxy0)
   use variables
   use hydro_parameters
   use heating_layer
   implicit none
   double precision, intent(in) :: x0
   double precision, dimension(1:2), intent(in) :: y0
   double precision, dimension(1:2), intent(inout) :: dxy0  
   ! y(1) : v, velocity of the flow (NB : v< 0)
   ! y(2) : c, sound speed
   ! x : coordinate along the x axis

   !-----------variables locales---------	
   double precision :: rhov, rho, P, grav, cool
  
   call give_gravity(x0, grav)
   rhov = -Mup
   rho = rhov/y0(1)
   P = rho*y0(2)**2/gamma
   call give_cool(x0,rho,P,cool)
  

   dxy0(1) = y0(1)/(y0(2)**2-y0(1)**2)*( -grav + (gamma-1.d0)*cool/rhov )
   dxy0(2) = (gamma-1)/2.d0 * y0(2)/(y0(2)**2-y0(1)**2) * &
      (grav + (1.d0 - gamma*y0(1)**2/y0(2)**2)*cool/rhov )

   return 
end subroutine derive0

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine rkutta0(x0,y0,dx0,ncount)
   use hydro_parameters
   real(dp) :: y0(2)
   real(dp) ::x0,dx0
   integer*4 ::ncount,nvar0
   real(dp) ::a(2,2),b(2),c(2)
   !--------------------variables locales--------------	
   integer*4 :: i,n,icount,j
   real(dp) :: yst(2),dxy(2)
   real(dp) :: k(2,2),knew(2,2)
   real(dp) :: dk(2,2),sup,xst(2)
   !-----------------------------------------------
  
   call initrk(a,b,c)
   nst = 2
   nvar0 = 2
  
   do i = 1,nst
      xst(i) = x0+c(i)*dx0
      do n = 1,nvar0
         k(i,n) = 0.d0
      enddo
   enddo
  
   icount = 0
3  sup = 0.d0
   do i = 1,nst
      do n = 1,nvar0
         yst(n) = y0(n)
      enddo
      do n = 1,nvar0
         do j = 1,nst
            yst(n) = yst(n)+a(i,j)*k(j,n)
         enddo
      enddo
      call derive0(xst(i),yst,dxy)
     
      do n = 1,nvar0
         knew(i,n) = dx0*dxy(n)
         dk(i,n)   = abs(knew(i,n)-k(i,n))
         k(i,n)    = knew(i,n)
         if (dk(i,n) .gt. sup) sup = dk(i,n)
      enddo
   enddo
  
   icount = icount + 1
   if ((sup .ge. 1.d-10) .and. (icount .lt. maxcount)) goto 3
  
   ncount = icount
   do n = 1,nvar0
      do i = 1,nst
         y0(n) = y0(n)+b(i)*k(i,n)
      enddo
   enddo

   return
end subroutine rkutta0

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine initRK(a,b,c)
   use hydro_parameters
   real(dp) :: a(2,2),b(2),c(2)
   !------variables locales------------------	
   real(dp) :: p,q,r,t
   !-------------------------------------	
  
   ! coeffs. de la methode implicite d'ordre 4,nst = 2 
   a(1,1) = 0.25d0
   a(1,2) = 0.25d0 - dsqrt(3.d0)/6.d0
   a(2,1) = 0.25d0 + dsqrt(3.d0)/6.d0
   a(2,2) = 0.25d0
  
   b(1) = 0.5d0
   b(2) = 0.5d0
  
   c(1) = 0.5d0 - dsqrt(3.d0)/6.d0
   c(2) = 0.5d0 + dsqrt(3.d0)/6.d0

   return
end subroutine initRK

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

