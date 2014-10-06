! =============================================================================

module steady_state

   use amr_parameters, only : dp

   implicit none

   ! Input parameters ---------------------------------------------------------
   ! The upstream density and sound speed are irrelevant to the problem; they
   ! change the values that the code works with, but the physical problem is
   ! unchanged.  Therefore, instead of making them values to be changed by the
   ! inputs file, I am hard-coding them as parameters.  But for generality and
   ! clarity of the expressions, I still include them.
   real(dp)            :: mach_up            ! The upstream Mach number
   real(dp), parameter :: dens_up = 1.0d0    ! The upstream density
   real(dp), parameter :: csnd_up = 1.0d0    ! The upstream sound speed

   ! The steady state (1D array) ----------------------------------------------
   real(dp), dimension(:,:), allocatable :: ss0    ! the storage array
   integer, parameter :: SS_DENS = 1   ! index of density
   integer, parameter :: SS_MOMX = 2   ! index of velocity
   integer, parameter :: SS_ENER = 3   ! index of energy
   integer, parameter :: SS_NVAR = 3   ! number of variables

contains

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine ss_build()
      use amr_parameters, only : dp
      use variables, only : uin, x, dx
      use hydro_parameters, only : gamma, xmax, iu1, iu2
      implicit none 

      ! =======================================================================
      ! Declare variables -----------------------------------------------------

      ! Arguments -------------------------------------------------------------

      ! Locals ----------------------------------------------------------------

      integer :: i
      real(dp) :: v, c, dens, ener, rhov
      integer, parameter :: LUN = 16

      !  RK looping
      integer :: irk
      integer :: Nskip
      integer, parameter :: NRK_2 = 5000
      integer, parameter :: NRK = NRK_2 * 2  ! because NRK must be even
      integer :: ncount                      ! RK looping
      real(dp) :: dx0, x0
      real(dp) :: y0(2)
      character(len=17) :: fname

      ! =======================================================================
      ! Fill the unperturbed state --------------------------------------------
      ! uin   (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)
      ! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)

      ! Integrate from xmax to x(iu1) (the lowest x-value in this slice)
      x0    = xmax
      y0(1) = -mach_up * csnd_up ! v
      y0(2) =            csnd_up ! c
      ! NRK = number of RK steps between x(i) and x(i-1)
      dx0   = -1.0d0 * dx / NRK

      ! Integrate down to x(iu2)
      if (xmax > x(iu2)) then
         ! If xmax > x(iu2), i.e. we are in a slice that is not along the
         ! upper x-boundary of the domain, then we need to integrate down from
         ! xmax to x(iu2).  In this case (xmax-x(iu2))/dx is a half-integer
         ! (because the x(i) values are cell centers and xmax is the upper edge
         ! of the uppermost non-guard cells in the x-direction).  Therefore we
         ! can integrate NRK/2 steps to reach a cell center, then integrate
         ! NRK*Nskip times, where Nskip is the number of cells within the
         ! domain above x(iu2).  We can compute Nskip as
         ! floor((xmax - x(iu2))/dx), because (xmax - x(iu2))/dx = Nskip+0.5.
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
      rhov = -dens_up * mach_up * csnd_up
      dens = rhov / v
      ener = dens * c**2 / (gamma*(gamma-1.0d0)) + 0.5d0 * dens * v**2

      ! Fill all cells at x = x(iu2)
      uin(:,iu2,:,:,:) = 0.0d0
      ss0(iu2,SS_DENS) = dens
      uin(:,iu2,:,:,1) = dens
      ss0(iu2,SS_MOMX) = rhov
      uin(:,iu2,:,:,2) = rhov
      ss0(iu2,SS_ENER) = ener
      uin(:,iu2,:,:,5) = ener

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
         rhov = -dens_up * mach_up * csnd_up
         dens = rhov / v
         ener = dens * c**2 / (gamma*(gamma-1.0d0)) + 0.5d0 * dens * v**2

         ! Fill all cells at x = x(i)
         uin(:,i,:,:,:) = 0.0d0
         ss0(i,SS_DENS) = dens
         uin(:,i,:,:,1) = dens
         ss0(i,SS_MOMX) = rhov
         uin(:,i,:,:,2) = rhov
         ss0(i,SS_ENER) = ener
         uin(:,i,:,:,5) = ener

      end do

      return
   end subroutine ss_build

   ! ==========================================================================

   subroutine derive0(x0,y0,dxy0)
      use amr_parameters, only : dp
      use variables
      use hydro_parameters
      use buoyant_layer, only : give_gravity, give_cool
      implicit none
      real(dp), intent(in) :: x0
      real(dp), dimension(1:2), intent(in) :: y0
      real(dp), dimension(1:2), intent(inout) :: dxy0  
      ! y(1) : v, velocity of the flow (NB : v< 0)
      ! y(2) : c, sound speed
      ! x : coordinate along the x axis

      !-----------variables locales---------	
      real(dp) :: rhov, rho, P, grav, cool

      call give_gravity(x0, grav)
      rhov = -dens_up * mach_up * csnd_up
      rho = rhov/y0(1)
      P = rho*y0(2)**2/gamma
      call give_cool(x0,rho,P,cool)


      dxy0(1) = y0(1)/(y0(2)**2-y0(1)**2)*( -grav + (gamma-1.d0)*cool/rhov )
      dxy0(2) = (gamma-1)/2.d0 * y0(2)/(y0(2)**2-y0(1)**2) * &
         (grav + (1.d0 - gamma*y0(1)**2/y0(2)**2)*cool/rhov )

      return 
   end subroutine derive0

   ! ==========================================================================

   subroutine rkutta0(x0,y0,dx0,ncount)
      use hydro_parameters
      real(dp) :: y0(2)
      real(dp) ::x0,dx0
      integer*4 ::ncount,nvar0
      real(dp) ::a(2,2),b(2),c(2)
      !--------------------variables locales--------------	
      integer*4 :: i,n,icount,j,nst
      integer*4, parameter :: maxcount = 1000000
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

   ! ==========================================================================

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

   ! --------------------------------------------------------------------------
   ! ==========================================================================
   ! --------------------------------------------------------------------------

   subroutine ss_write_to_file()

      use amr_parameters, only : nx
      use hydro_parameters, only : iu1, iu2
#     if WITHMPI==1
      use mpi_var, only : nxslice, xposition, yposition, zposition
#     endif
      use variables, only : x

      implicit none

#     if WITHMPI==1
#     include "mpif.h"
#     endif

      integer :: pos, idx, var
      integer, parameter :: LUN = 37
      character(len=*), parameter :: fname = "stable_state.dat"
      integer :: ierr

#     if WITHMPI==1
      if ((xposition == 0) .and. (yposition == 0) .and. (zposition == 0)) then
#     endif
         open(unit=LUN, file=fname, status='unknown')
         write(LUN,"('#',x,4(3x,a30))") "radius", "density", "x-momentum", &
                                        "total_energy_dens"
         close(LUN)
#     if WITHMPI==1
      end if
      ! synchronise - only one processor writes at a time
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#     endif

      ! Write the data in order (lowest value to highest)
#     if WITHMPI==1
      do pos = 0, nxslice
         if ((xposition == pos) .and. &
             (yposition == 0) .and. (zposition == 0)) then
#     endif
            open(unit=LUN, file=fname, access="append", status="old")
            do idx = 1, nx
               write(LUN,"(5x,es30.23)",advance="no") x(idx)
               do var = 1, SS_NVAR
                  write(LUN,"(3x,es30.23)",advance="no") ss0(idx,var)
               end do
               write (LUN,*) ""
            end do
            close(LUN)
#     if WITHMPI==1
         end if
         ! synchronize - wait your turn to write so file is in order
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
      end do
#     endif

      return
   end subroutine ss_write_to_file

end module steady_state
