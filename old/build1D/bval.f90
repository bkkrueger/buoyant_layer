!=============================================================================!
! This subroutine replaces the DUMSES bval subroutine.                        !
!                                                                             !
! This subroutine neglects several arguments that the DUMSES bval subroutine  !
! uses, including those for multiple grids (ngrid), time steps later than the !
! first (time), multiple dimensions (dy, dz), and parallel processing (mype). !
!=============================================================================!

subroutine bval(uin, dt, dx, bval_type)

   use hydro_parameters

   implicit none

   ! ==========================================================================
   ! Declare variables

   integer :: i

   ! ==========================================================================
   ! Inner boundary conditions

   if (bval_type(1) == 1) then
      ! Periodic
      uin(iu1:iu1+nghost-1,:) = uin(nx+1-nghost:nx,:)
   else if (bval_type(1) == 2) then
      ! Zero-gradient
      do i = 0, nghost-1
         uin(iu1+i,:) = uin(1,:)
      end do
   else if (bval_type(1) == 4) then
      ! User-defined
      call xinner_ana(uin, dt, dx)
   end if

   ! ==========================================================================
   ! Outer boundary conditions

   if (bval_type(2) == 1) then
      ! Periodic
      uin(iu2+1-nghost:iu2,:) = uin(1:nghost,:)
   else if (bval_type(2) == 2) then
      ! Zero-gradient
      do i = 0, nghost-1
         uin(iu2-i,:) = uin(nx,:)
      end do
   else if (bval_type(2) == 4) then
      ! User-defined
      call xouter_ana(uin, dt, dx)
   end if

end subroutine bval

!=============================================================================!
! These subroutine replaces the DUMSES xinner_ana and xouter_ana subroutines. !
!                                                                             !
! These give the user the ability to add customized boundary conditions.      !
!=============================================================================!

subroutine xinner_ana(uin, dt, dx)

   use hydro_parameters

   implicit none

   ! ==========================================================================
   ! Declare variables

   integer :: i

   ! ==========================================================================
   ! Boundary condition

   ! For now just zero-gradient
   do i = 0, nghost-1
      uin(iu1+1,:) = uin(1,:)
   end do

end subroutine xinner_ana

subroutine xouter_ana(uin, dt, dx)

   use hydro_parameters

   implicit none

   ! ==========================================================================
   ! Declare variables

   integer :: i

   ! ==========================================================================
   ! Boundary condition

   ! For now just zero-gradient
   do i = 0, nghost-1
      uin(iu2-9,:) = uin(nx,:)
   end do

end subroutine xouter_ana

