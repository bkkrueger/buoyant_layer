subroutine xinner_ana
   use amr_parameters, only : nx
   use hydro_parameters, only : iu1, iu2, ju1, ju2, ku1, ku2
   use steady_state, only : ss0, SS_DENS, SS_MOMX, SS_ENER
   use variables, only : uin

   implicit none

   ! Zero-gradient, fixed outflow set by steady-state
   uin(1,iu1:0,ju1:ju2,ku1:ku2,:) = 0.0d0
   uin(1,iu1:0,ju1:ju2,ku1:ku2,1) = ss0(nx,SS_DENS)
   uin(1,iu1:0,ju1:ju2,ku1:ku2,2) = ss0(nx,SS_MOMX)
   uin(1,iu1:0,ju1:ju2,ku1:ku2,5) = ss0(nx,SS_ENER)

   !integer :: j, k
   ! TODO : Do I want to force zero horizontal velocity?
   !do k=ku1,ku2
   !   do j=ju1,ju2      
   !      uin(1,iu1  ,j,k,:) = uin(1,iu1+3,j,k,:)
   !      uin(1,iu1+1,j,k,:) = uin(1,iu1+3,j,k,:)
   !      uin(1,iu1+2,j,k,:) = uin(1,iu1+3,j,k,:)
   !   end do
   !end do

   return
end subroutine xinner_ana
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE xouter_ana
   use amr_parameters, only : dp, nx
   use hydro_parameters, only : gamma, iu1, iu2, ju1, ju2, ku1, ku2
   use perturb, only : pert_type, pert_amplitude, BC_drive_pert, osc_omega, &
                       rand_amplitude, iseed, ran2, osc_slope
   use steady_state, only : mach_up, dens_up, csnd_up, &
                            ss0, SS_DENS, SS_MOMX, SS_ENER
   use variables, only : uin, time, x, y
   implicit none

   integer :: i,j,k
   real(dp) :: velx_up, pres_up
   real(dp) :: perturbation, rand, x_eff, dens

   ! Unperturbed state:
   ! dens_up is given
   velx_up = -1.0d0 * mach_up * csnd_up
   pres_up = dens_up * csnd_up**2 / gamma

   if ((BC_drive_pert) .and. (trim(pert_type) == 'sasi')) then
      ! Wave perturbation around fixed inflow
      do i = 1, iu2
      !do i = nx+1, iu2
         do j = ju1, ju2
            do k = ku1, ku2
               ! Perturbation:
               ! -- scaled by the perturbation amplitude
               perturbation = pert_amplitude
               ! -- with a sinusoidal component
               x_eff = x(i) - osc_slope*y(j) - velx_up*time
               perturbation = perturbation * sin(osc_omega * x_eff)
               ! -- with a random component scaled by the random amplitude
               call ran2(iseed, rand)
               rand = 2.0d0*rand - 1.0d0
               perturbation = perturbation * (1.0d0 + rand_amplitude*rand)
               ! -- use it to perturb the density
               dens = dens_up * (1.0d0 + perturbation)
               ! Fill the state
               uin(1,i,j,k,:) = 0.0d0
               uin(1,i,j,k,1) = dens
               uin(1,i,j,k,2) = dens * velx_up
               uin(1,i,j,k,5) = pres_up/(gamma-1.0d0) + 0.50d0*dens*velx_up**2
            end do
         end do
      end do
   else
      ! Fixed inflow
      uin(1,nx+1:iu2,ju1:ju2,ku1:ku2,:) = 0.0d0
      uin(1,nx+1:iu2,ju1:ju2,ku1:ku2,1) = ss0(1,SS_DENS)
      uin(1,nx+1:iu2,ju1:ju2,ku1:ku2,2) = ss0(1,SS_MOMX)
      uin(1,nx+1:iu2,ju1:ju2,ku1:ku2,5) = ss0(1,SS_ENER)
   end if

   RETURN
END SUBROUTINE xouter_ana
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE yinner_ana
#if NDIM>1
  USE hydro_parameters
  USE variables
  IMPLICIT NONE

  INTEGER ::k,i

  ! Periodic BC
  DO k=ku1,ku2
     DO i=iu1,iu2
        uin(1,i,ju1  ,k,:) = uin(1,i,ju1+3,k,:)
        uin(1,i,ju1+1,k,:) = uin(1,i,ju1+3,k,:)
        uin(1,i,ju1+2,k,:) = uin(1,i,ju1+3,k,:)
     END DO
  END DO
!
#endif
  RETURN
END SUBROUTINE yinner_ana
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE youter_ana
#if NDIM>1
  USE hydro_parameters
  USE variables
  IMPLICIT NONE

  INTEGER ::k,i

  DO k=ku1,ku2
     DO i=iu1,iu2
        uin(1,i,ju2-2,k,:) = uin(1,i,ju2-3,k,:)
        uin(1,i,ju2-1,k,:) = uin(1,i,ju2-3,k,:)
        uin(1,i,ju2  ,k,:) = uin(1,i,ju2-3,k,:)
     END DO
  END DO
!
#endif
  RETURN
END SUBROUTINE youter_ana
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE zinner_ana
#if NDIM > 2
  USE hydro_parameters
  USE variables
  IMPLICIT NONE

  INTEGER ::j,i

  ! Periodic BC
  DO j=ju1,ju2
     DO i=iu1,iu2
        uin(1,i,j,ku1  ,:) = uin(1,i,j,ku1+3,:)
        uin(1,i,j,ku1+1,:) = uin(1,i,j,ku1+3,:)
        uin(1,i,j,ku1+2,:) = uin(1,i,j,ku1+3,:)
     END DO
  END DO
!
#endif 
  RETURN
END SUBROUTINE zinner_ana
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE zouter_ana
#if NDIM > 2
  USE hydro_parameters
  USE variables
  IMPLICIT NONE

  INTEGER ::j,i

  ! Periodic BC
  DO j=ju1,ju2
     DO i=iu1,iu2
        uin(1,i,j,ku2-2,:) = uin(1,i,j,ku2-3,:)
        uin(1,i,j,ku2-1,:) = uin(1,i,j,ku2-3,:)
        uin(1,i,j,ku2  ,:) = uin(1,i,j,ku2-3,:)
     END DO
  END DO
!
#endif
  RETURN
END SUBROUTINE zouter_ana
