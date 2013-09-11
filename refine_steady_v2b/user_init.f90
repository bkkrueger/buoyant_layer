subroutine user_init(mype, npes)
   use variables, only : time, dt, restart
   implicit none
   integer, intent(in) :: mype, npes
   if (restart == 0) then
      call refine_initial_model(mype, npes)
      call perturb_initial_model(mype, npes)
      time = 0.0d0
      dt   = 0.0d0
   end if
   return
end subroutine user_init
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_eta(eta,r)
  use amr_parameters
  implicit none
  real(dp) :: eta,r
  eta=0.d0
  return
end subroutine get_eta
