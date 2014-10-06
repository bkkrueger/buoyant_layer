subroutine history(uin,x,time,dt,dx,dy,dz,mype,npes)

   use amr_parameters, only : dp, nvector, nx, ny, nz, verbose
   use hydro_parameters, only : iu1, iu2, ju1, ju2, ku1, ku2, nvar, gamma
   use variables, only : cartesian
   use buoyant_layer, only : give_gravity
   use perturb, only : pert_type, bubble_x, x_mass_old, x_advc_old,           &
                       time_old, dv_grav_old, dv_buoy_old, dv_drag_old
   use steady_state, only : ss0, SS_DENS, SS_MOMX, SS_ENER

   implicit none

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   ! ==========================================================================
   ! Declare variables --------------------------------------------------------

   ! Arguments ----------------------------------------------------------------
   integer :: mype, npes
   real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) :: uin
   real(dp), dimension(iu1:iu2) :: x
   real(dp) :: time, dt, dx, dy, dz

   ! Locals -------------------------------------------------------------------
   character(len=255) :: filename='history.txt'
   integer, parameter :: LUN = 17
#  if WITHMPI==1
   real(dp), dimension(5) :: sendbuf, recvbuf
   integer :: ierr
#  endif
   integer :: l, i, j, k
   integer :: ilo, ihi, jlo, jhi, klo, khi
   real(dp) :: cell_vol, xsection   ! cell volume and horizontal cross-section
   real(dp) :: bubble(iu1:iu2,ju1:ju2,ku1:ku2), bubble_max
   real(dp) :: temp, contrast
   real(dp), parameter :: BUBBLE_CUTOFF = 1.0d-1
   real(dp) :: sum_dV, sum_x_dV, sum_b_dV, sum_x_b_dV
   real(dp) :: x_geom, x_mass, x_advc
   real(dp) :: sum_dA_in, sum_dA_out, sum_v_dA
   real(dp) :: xsect_vol
   real(dp) :: v_bubble, v_bgnd, v_advc
   real(dp) :: sum_dens_dV, sum_dens_g_dV, sum_rho0_g_dV
   real(dp) :: dens, rho0, grav
   real(dp) :: a_grav, a_buoy, a_drag
   real(dp), parameter :: C_D = 0.5d0
   real(dp) :: dv_grav, dv_buoy, dv_drag

   ! ==========================================================================
   ! Compute the next line of the history file --------------------------------

   if (verbose) write (*,*) 'Entering history...'

   if (trim(pert_type) /= 'bubble') then
      return
   end if

   ! Some dimensional juggling ------------------------------------------------
   ! - We only want to sum over body cells, not guard cells (or we will
   !   double-count certain cells and also include the semi-artificial BCs)
   ! - The volume depends on the dimensionality

   if (.not. cartesian) then
      ! I have not yet implemented cylindrical or spherical for this file, so
      ! for now insist that the geometry be Cartesian
      write(*,*) "HISTORY : Only Cartesian geometry supported."
      stop
   end if

   ilo = 1 ;   jlo = 1 ;   klo = 1
   ihi = nx;   jhi = ny;   khi = nz

   cell_vol = dx
   xsection = 1.0d0
#if NDIM > 1
   cell_vol = cell_vol * dy
   xsection = xsection * dy
#if NDIM > 2
   cell_vol = cell_vol * dz
   xsection = xsection * dz
#endif
#endif

   ! Define the bubble --------------------------------------------------------
   do i = iu1, iu2
      !bubble(i,:,:) = abs((uin(1,i,:,:,1) - ss0(i,SS_DENS)) / &
      !                    ss0(i,SS_DENS))
      bubble(i,:,:) = uin(1,i,:,:,2)**2+uin(1,i,:,:,3)**2+uin(1,i,:,:,4)**2
      bubble(i,:,:) = bubble(i,:,:) / (2.0d0 * uin(1,i,:,:,1))
      bubble(i,:,:) = uin(1,i,:,:,5) - bubble(i,:,:)
      bubble(i,:,:) = bubble(i,:,:) / &
                      (ss0(i,SS_ENER) - 0.5d0*ss0(i,SS_MOMX)/ss0(i,SS_DENS))
      bubble(i,:,:) = bubble(i,:,:) * &
                      (uin(1,i,:,:,1)/ss0(i,SS_DENS))**(-1.0d0*gamma)
      bubble(i,:,:) = log(bubble(i,:,:)) / (gamma - 1.0d0)
   end do
   bubble_max = maxval(bubble)

   ! Track the peak density contrast
   contrast = 0.0d0
   do i = iu1, iu2
      temp = maxval(abs(uin(1,i,:,:,1)-ss0(i,SS_DENS))/ss0(i,SS_DENS))
      if (temp > contrast) then
         contrast = temp
      end if
   end do

#  if WITHMPI==1
   sendbuf(1) = bubble_max
   sendbuf(2) = contrast
   call MPI_Allreduce(sendbuf, recvbuf, 2, MPI_DOUBLE_PRECISION, MPI_MAX,     &
                      MPI_COMM_WORLD, ierr)
   bubble_max = recvbuf(1)
   contrast = recvbuf(1)
#  endif

   ! Positions ----------------------------------------------------------------

   ! Local summations
   sum_dV     = 0.0d0
   sum_x_dV   = 0.0d0
   sum_b_dV   = 0.0d0
   sum_x_b_dV = 0.0d0
   do k = klo, khi
      do j = jlo, jhi
         do i = ilo, ihi
            if (bubble(i,j,k) > BUBBLE_CUTOFF * bubble_max) then
               sum_dV     = sum_dV     + cell_vol
               sum_x_dV   = sum_x_dV   + cell_vol * x(i)
               sum_b_dV   = sum_b_dV   + cell_vol * bubble(i,j,k)
               sum_x_b_dV = sum_x_b_dV + cell_vol * bubble(i,j,k) * x(i)
            end if
         end do
      end do
   end do

   ! Communication
#  if WITHMPI==1
   sendbuf(1:4) = (/ sum_dV, sum_x_dV, sum_b_dV, sum_x_b_dV /)
   call MPI_Allreduce(sendbuf, recvbuf, 4, MPI_DOUBLE_PRECISION, MPI_SUM,     &
                      MPI_COMM_WORLD, ierr)
   sum_dV     = recvbuf(1)
   sum_x_dV   = recvbuf(2)
   sum_b_dV   = recvbuf(3)
   sum_x_b_dV = recvbuf(4)
#  endif

   ! Compute bubble positions
   x_geom = sum_x_dV   / sum_dV
   x_mass = sum_x_b_dV / sum_b_dV

   ! Compute bubble velocity
   if (time == 0.0d0) then
      v_bubble = 0.0d0
   else
      ! formally should probably be x_geom, but x_mass is smoother while x_geom
      ! is some sort of strange step function that gives terrible v_bubble
      v_bubble = (x_mass - x_mass_old) / (time - time_old)
   end if

   ! Compute advection position
   if (time == 0.0d0) then
      x_advc = bubble_x
   else
      x_advc = x_advc_old + v_advc(x_advc_old,x,dx) * (time - time_old)
   end if

   ! Cycle new->old
   x_mass_old = x_mass
   x_advc_old = x_advc

   ! Cross-section and mean background velocity -------------------------------

   ! Local summations
   sum_dA_in  = 0.0d0
   sum_v_dA   = 0.0d0
   sum_dA_out = 0.0d0
   do k = klo, khi
      do j = jlo, jhi
         do i = ilo, ihi
            if ((x(i)-0.5d0*dx <= x_geom) .and.  (x_geom < x(i)+0.5d0*dx)) then
               if (bubble(i,j,k) > BUBBLE_CUTOFF * bubble_max) then
                  sum_dA_in  = sum_dA_in  + xsection
               else
                  sum_v_dA   = sum_v_dA   + xsection * &
                                            uin(1,i,j,k,2) / uin(1,i,j,k,1)
                  sum_dA_out = sum_dA_out + xsection
               end if
            end if
         end do
      end do
   end do

   ! Communication
#  if WITHMPI==1
   sendbuf(1) = sum_dA_in
   sendbuf(2) = sum_v_dA
   sendbuf(3) = sum_dA_out
   call MPI_Allreduce(sendbuf, recvbuf, 3, MPI_DOUBLE_PRECISION, MPI_SUM,     &
                      MPI_COMM_WORLD, ierr)
   sum_dA_in  = recvbuf(1)
   sum_v_dA   = recvbuf(2)
   sum_dA_out = recvbuf(3)
#  endif

   ! Compute cross-section/volume ratio
   xsect_vol = sum_dA_in / sum_dV

   ! Compute background velocity at bubble elevation
   v_bgnd = sum_v_dA / sum_dA_out

   ! Accelerations ------------------------------------------------------------

   ! Local summations
   sum_dens_dV   = 0.0d0
   sum_dens_g_dV = 0.0d0
   sum_rho0_g_dV = 0.0d0
   do k = klo, khi
      do j = jlo, jhi
         do i = ilo, ihi
            dens = uin(1,i,j,k,1)
            rho0 = ss0(i,SS_DENS)
            call give_gravity(x(i), grav)
            if (bubble(i,j,k) > BUBBLE_CUTOFF * bubble_max) then
               sum_dens_dV   = sum_dens_dV   + cell_vol * dens
               sum_dens_g_dV = sum_dens_g_dV + cell_vol * dens * grav
               sum_rho0_g_dV = sum_rho0_g_dV + cell_vol * rho0 * grav
            end if
         end do
      end do
   end do

   ! Communication
#  if WITHMPI==1
   sendbuf(1:3) = (/ sum_dens_dV, sum_dens_g_dV, sum_rho0_g_dV /)
   call MPI_Allreduce(sendbuf, recvbuf, 3, MPI_DOUBLE_PRECISION, MPI_SUM,     &
                      MPI_COMM_WORLD, ierr)
   sum_dens_dV   = recvbuf(1)
   sum_dens_g_dV = recvbuf(2)
   sum_rho0_g_dV = recvbuf(3)
#  endif

   ! Compute accelerations from sums
   a_grav = sum_dens_g_dV / sum_dens_dV
   a_buoy = sum_rho0_g_dV / sum_dens_dV

   ! Compute drag acceleration
   a_drag = 0.5d0 * (v_bgnd - v_bubble)**2 * C_D * xsect_vol

   ! Accumulate velocity changes
   if (time == 0.0d0) then
      dv_grav = a_grav * time
      dv_buoy = a_buoy * time
      dv_drag = a_drag * time
   else
      dv_grav = dv_grav_old + a_grav * (time - time_old)
      dv_buoy = dv_buoy_old + a_buoy * (time - time_old)
      dv_drag = dv_drag_old + a_drag * (time - time_old)
   end if
   dv_grav_old = dv_grav
   dv_buoy_old = dv_buoy
   dv_drag_old = dv_drag

   ! Cycle the time -----------------------------------------------------------
   time_old = time

   ! Print the next line of the history file ----------------------------------
   if (mype == 0) then
      if (time == 0.d0) then
         open(unit=LUN,file=filename,status='unknown',position='append')
         write(LUN,"('#',x,a13,50(3x,a13))") "time", "bubble_max",            &
            "X-section/volume",                                               &
            "x_advc", "x_geom", "x_mass", "v_bubble", "v_bgnd",               &
            "a_grav", "a_buoy", "a_drag", "dv_grav", "dv_buoy", "dv_drag",    &
            "contrast"
      else
         open(unit=LUN,file=filename,status='old',position='append')
      endif
      write(LUN,"(x,x,es13.6,50(3x,es13.6))") time, bubble_max,               &
         xsect_vol,                                                           &
         x_advc, x_geom, x_mass, v_bubble, v_bgnd,                            &
         a_grav, a_buoy, a_drag, dv_grav, dv_buoy, dv_drag,                   &
         contrast
      close(LUN)
   endif

   return
end subroutine history

! -----------------------------------------------------------------------------
! =============================================================================
! -----------------------------------------------------------------------------

function v_advc(x_advc, x, dx) result(v)
   use amr_parameters, only : dp, nx
   use hydro_parameters, only : iu1, iu2
#  if WITHMPI==1
   use mpi_var, only : yposition, zposition
#  endif
   use steady_state, only : ss0, SS_DENS, SS_MOMX
   implicit none
   real(dp), intent(in) :: x_advc, x(iu1:iu2), dx
   real(dp)             :: v, tmp
   integer              :: i, ierr

#  if WITHMPI==1
#  include "mpif.h"
#  endif

   v = 0.0d0
#  if WITHMPI==1
   if ((yposition == 0) .and. (zposition == 0)) then
#  endif
      do i = 1, nx
         if ((x(i)-0.5d0*dx <= x_advc) .and. (x_advc < x(i)+0.5d0*dx)) then
            v = ss0(i,SS_MOMX) / ss0(i,SS_DENS)
         end if
      end do
#  if WITHMPI==1
   end if
#  endif

#  if WITHMPI==1
   call MPI_Allreduce(v, tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM,               &
                      MPI_COMM_WORLD, ierr)
   v = tmp
#  endif

   return
end function v_advc

