! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

subroutine condinit(mype)
   use variables
   use hydro_parameters
   use heating_layer
   use steady
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

   ! Physics Parameters namelist
   namelist /physics_params/ Mup, Kheat, Kgrav, gamma, layer_shape, &
      amprand, pert, spectrum, nxmin, nxmax, nymin, nymax, pertseed, &
      xpmin, xpmax, layer_limit, amp_min, do_refine

   ! ==========================================================================
   ! Fill the initial conditions

   allocate(ss0(iu1:iu2,SS_nv))

   ! Default values -----------------------------------------------------------

   ! Background state
   ! c_up = 1
   ! rho_up = 1
   Mup   = 1.0d-1
   Kheat = 1.0d-2
   Kgrav = 1.0d0
   gamma = 4.0d0 / 3.0d0

   ! Buoyant layer
   layer_shape = 'trapezoid'
   layer_limit = 1.0d0  ! layer extends from -layer_limit to +layer_limit

   ! Perturbations
   pertseed  = 2
   xpmin     = -1.0d0 * layer_limit
   xpmax     =          layer_limit
   amprand   = 0.0d0
   amp_min   = 0.0d0
   pert      = 'rho'
   spectrum  = 'flat'
   nxmin     = 0
   nxmax     = 0
   nymin     = 0
   nymax     = 0
   do_refine = .true.
  
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
      write(*,*) ""
      write(*,*) "perturbation parameters"
      write(temp_str,"(l13)") do_refine
      write(*,f_str) "do_refine", trim(temp_str)
      write(temp_str,"(es13.6)") amprand
      write(*,f_str) "amprand",   trim(temp_str)
      write(temp_str,"(es13.6)") amp_min
      write(*,f_str) "amp_min",   trim(temp_str)
      write(*,f_str) "pert", trim(pert)
      write(*,f_str) "spectrum", trim(spectrum)
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
      write(temp_str,"(es13.6)") xpmin
      write(*,f_str) "xpmin",   trim(temp_str)
      write(temp_str,"(es13.6)") xpmax
      write(*,f_str) "xpmax",   trim(temp_str)
      write(*,*) ""
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
   ss0(iu2,1) = dens
   ss0(iu2,2) = -Mup
   ss0(iu2,3) = ener

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
      ss0(i,1) = dens
      ss0(i,2) = -Mup
      ss0(i,3) = ener

   end do

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

