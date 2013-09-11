! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)
SUBROUTINE source_term
  USE hydro_parameters
  USE const
  USE variables
  USE heating_layer
  IMPLICIT NONE

  INTEGER :: i,j,k,l
  REAL(dp) :: rho, grav, rhov, P, Pold, cool,cool_old
#ifdef FORCED
  REAL(dp) :: alpha,c_over_u_sq,pi,dsx
#endif

  IF (verbose) write (*,*) 'Entering source_term...'
  
 DO l=1,ngrid
     DO k=1,nz
        DO j=1,ny
           DO i=1,nx
              ! gravity :
              CALL give_gravity(x(i),grav)
              rho=half*(old_uin(l,i,j,k,1)+uin(l,i,j,k,1))
              rhov=half*(old_uin(l,i,j,k,2)+uin(l,i,j,k,2)) + half*rho*grav*dt ! predictive step of gravity added here..           
              ! cooling :
              Pold = (gamma-1)*( old_uin(l,i,j,k,5) &
                   - 0.5d0*(old_uin(l,i,j,k,2)**2+old_uin(l,i,j,k,3)**2+old_uin(l,i,j,k,4)**2)/old_uin(l,i,j,k,1) &
                   - 0.5d0*(old_uin(l,i,j,k,6)**2+old_uin(l,i,j,k,7)**2+old_uin(l,i,j,k,8)**2) )
              CALL give_cool(x(i),old_uin(l,i,j,k,1),Pold,cool_old) ! compute the cooling function at t for the predictive step
              P = (gamma-1)*( uin(l,i,j,k,5) &
                   - 0.5d0*(uin(l,i,j,k,2)**2+uin(l,i,j,k,3)**2+uin(l,i,j,k,4)**2)/uin(l,i,j,k,1) &
                   - 0.5d0*(uin(l,i,j,k,6)**2+uin(l,i,j,k,7)**2+uin(l,i,j,k,8)**2) )
              P = half*(Pold+P) + (gamma-1.d0)*half*cool_old*dt    ! predictive step of cooling added here..
              CALL give_cool(x(i),rho,P,cool) ! returns the cooling function at t + dt/2 in the variable "cool"
              ! update uin :
              uin(l,i,j,k,2) = uin(l,i,j,k,2) + rho*grav*dt
              uin(l,i,j,k,5) = uin(l,i,j,k,5) + rhov*grav*dt + cool*dt
           END DO
        END DO
     END DO
  END DO
CALL update_gravin	

 

#ifdef STRATIFIED
  DO l=1,ngrid
     DO k=1,nz
        DO j=1,ny
           DO i=1,nx
!                 rho=half*(old_uin(l,i,j,k,1)+uin(l,i,j,k,1))
                 rho=uin(l,i,j,k,1)
                 uin(l,i,j,k,3) = uin(l,i,j,k,3) - rho*Omega0**2*y(j)*dt
           END DO
        END DO
     END DO
  END DO
  CALL update_gravin
#endif

#ifdef FORCED

  pi=2.d0*ASIN(1.d0)

  alpha=1.25d0
  c_over_u_sq=0.3

  DO l=1,ngrid
      DO k=1,nz
         DO j=1,ny
            DO i=1,nx

               dsx =half*dt*alpha*ciso**2/c_over_u_sq * ( &
     &                 old_uin(l,i,j,k,1)*COS(2.*pi*((x(i)+0.5)-ciso/sqrt(c_over_u_sq)* &
     &                 (time   ))+2.*pi*(-1.42830e-3-0.04545123)) + &
     &                     uin(l,i,j,k,1)*COS(2.*pi*((x(i)+0.5)-ciso/sqrt(c_over_u_sq)* &
     &                 (time+dt))+2.*pi*(-1.42830e-3-0.04545123)) )


               uin(l,i,j,k,2)=uin(l,i,j,k,2) + dsx

            END DO
         END DO
      END DO
  END DO
  CALL update_gravin
#endif

END SUBROUTINE source_term
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE update_gravin
  USE variables
  USE hydro_parameters
  USE heating_layer
  IMPLICIT NONE

  INTEGER :: l,i,j,k
  REAL(dp) :: alphabis,c_over_u_sq,pi
  REAL(dp) ::  grav
  
  
  ! Marche de potentiel :
  DO l=1,ngrid
      DO k=ku1,ku2
         DO j=ju1,ju2
            DO i=iu1,iu2
               CALL give_gravity(x(i), grav)
               gravin(l,i,j,k,1)= grav
            END DO
         END DO
      END DO
  END DO


#ifdef STRATIFIED
  DO l=1,ngrid
      DO k=ku1,ku2
         DO j=ju1,ju2
            DO i=iu1,iu2
               gravin(l,i,j,k,1)=-Omega0**2*y(j)
            END DO
         END DO
      END DO
  END DO
#endif

#ifdef FORCED
  pi=2.d0*ASIN(1.d0)

  alphabis=1.25d0
  c_over_u_sq=0.3

  DO l=1,ngrid
      DO k=ku1,ku2
         DO j=ju1,ju2
            DO i=iu1,iu2

               gravin(l,i,j,k,1)=alphabis*ciso**2/c_over_u_sq*COS(2.*pi* &
     &                 ((x(i)+0.5)-ciso/sqrt(c_over_u_sq)* &
     &                 (time+dt))+2.*pi*(-1.42830e-3-0.04545123))

            END DO
         END DO
      END DO
  END DO
#endif

END SUBROUTINE update_gravin
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine gravity_predictor(v,gravin,dt,igrav,jgrav,kgrav,ngrid)
  use hydro_parameters
  use const
  implicit none

  integer :: ngrid
  logical :: igrav,jgrav,kgrav
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3   ) :: v
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim) :: gravin
  real(dp) :: dt

  integer :: i,j,k,l
  integer :: ilo,ihi,jlo,jhi,klo,khi

  ilo = min(1,iu1+1) ; ihi = max(1,iu2-1)
  jlo = min(1,ju1+1) ; jhi = max(1,ju2-1)
  klo = min(1,ku1+1) ; khi = max(1,ku2-1)

  !v=v+gravin*dt*half
  do k=klo,khi
     do j=jlo,jhi
        do i=ilo,ihi
           do l=1,ngrid

              if (igrav)                v(l,i,j,k,1) = v(l,i,j,k,1) + half*dt*gravin(l,i,j,k,1)
              if ((jgrav).and.(ndim>1)) v(l,i,j,k,2) = v(l,i,j,k,2) + half*dt*gravin(l,i,j,k,2)
              if ((kgrav).and.(ndim>2)) v(l,i,j,k,3) = v(l,i,j,k,3) + half*dt*gravin(l,i,j,k,3)

           end do
        end do
     end do
  end do

  return
end subroutine gravity_predictor 
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cooling_predictor(density,press)
  use hydro_parameters
  use const
  use variables

  implicit none

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::density
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::press
  real(dp) :: lcool,rho,P

  integer :: i,j,k,l
  integer :: ilo,ihi,jlo,jhi,klo,khi

  if(verbose) write(*,*) 'entering cooling_predictor...'
  ilo = min(1,iu1+1) ; ihi = max(1,iu2-1)
  jlo = min(1,ju1+1) ; jhi = max(1,ju2-1)
  klo = min(1,ku1+1) ; khi = max(1,ku2-1)

  !!Pressure update in t+dt/2
  do k=klo,khi
     do j=jlo,jhi
        do i=ilo,ihi
           do l=1,ngrid
              
              
              !!Cooling function in t+dt/2 (use t instead ?)
              CALL give_cool(x(i),density(l,i,j,k),press(l,i,j,k),lcool)
              press(l,i,j,k) = press(l,i,j,k) + 0.5*(gamma-1.d0)*lcool*dt
              
           end do
        end do
     end do
  end do
  
  if(verbose) write(*,*) 'Leaving cooling_predictor...'
  
  return
end subroutine cooling_predictor
