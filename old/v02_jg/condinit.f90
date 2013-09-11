MODULE HEATING_LAYER
  USE hydro_parameters
  IMPLICIT NONE
  REAL(DP) :: Mup, Kheat, Kgrav
  REAL(DP) :: first
  REAL(dp) :: amprand
  CHARACTER*20   pert,spectrum
  INTEGER  :: nxmin,nxmax,nymin,nymax,pertseed

END MODULE HEATING_LAYER
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE condinit(mype)
  USE variables
  USE hydro_parameters
  USE heating_layer
  IMPLICIT NONE 

  INTEGER :: mype

  INTEGER ::i,j,k,l,irk,nrk,ncount,iseed,irand,jrand
  real :: rvalue
  REAL(dp) ::B0, v, c, rho
  REAL(dp) ::pi
  DOUBLE PRECISION :: dx0, x0
  DOUBLE PRECISION, DIMENSION(1:2) :: y0
  DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE :: fourrier_amp,fourrier_phasex,fourrier_phasey

  ! Normalisation :
  ! c_up = 1
  ! rho_up = 1
  ! H = 1 (hauteur de la marche)

  
  ! Heating layer ::
  ! Read the physical parametters of the problems :
  namelist /physics_params/Mup, Kheat, Kgrav, gamma,amprand,pert,spectrum,nxmin,nxmax,nymin,nymax,pertseed
  OPEN(unit=1,file='input' ,status='old')
  READ(1,physics_params)
  CLOSE(1)
  CALL update_gravin

  first = 1.d0
  pi = 3.14159265359
  IF (pert .EQ. 'rho') THEN
     iseed=pertseed
     ALLOCATE(fourrier_amp(nxmin:nxmax,nymin:nymax),fourrier_phasex(nxmin:nxmax,nymin:nymax),fourrier_phasey(nxmin:nxmax,nymin:nymax))
     DO irand = nxmin,nxmax
        DO jrand = nymin,nymax
           IF (spectrum .EQ. 'flat') THEN
              fourrier_amp(irand,jrand) = amprand
           ELSE
              fourrier_amp(irand,jrand) = amprand*sqrt((nxmin**2/4.d0 + nymin**2/(ymax-ymin)**2)/(irand**2/4.d0 + jrand**2/(ymax-ymin)**2))
           END IF
           call ran2(iseed,rvalue)
           fourrier_phasex(irand,jrand) = 2.d0*pi*rvalue
           call ran2(iseed,rvalue)
           fourrier_phasey(irand,jrand) = 2.d0*pi*rvalue
        END DO
     END DO
  ENDIF
  !WRITE(*,*) 'pi :',pi, 'pertseed :',pertseed, 'iseed :',iseed
  !WRITE(*,*) 'amp :',fourrier_amp
  !WRITE(*,*) 'phase x :',fourrier_phasex
  
  !WRITE(*,*) 'phase y :',fourrier_phasey
  
  
  ! Calculate the initial state for uin...
  ! uin   (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)
  ! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)
  nrk = 10000
  pi=2.D0*ASIN(1.d0)
  B0 = 0.d0
  WRITE(*,*) 'Mup:', Mup, '  Kheat:', Kheat, '  Kgrav:', Kgrav
  iseed=mype
  
  
  
  
  y0(1) = -Mup  ! v
  y0(2) = 1.d0  ! c
  DO i=iu2,iu1,-1
     IF (i .EQ. iu2) THEN
        x0 = xmax
        dx0 = (x(0)-x(1))/nrk
     ENDIF
     DO WHILE(x0>x(i))
	CALL rkutta0(x0, y0, dx0, ncount)
        x0 = x0 + dx0
     END DO
     
     v = y0(1)
     c = y0(2)
     rho = -Mup/v
     DO l=1,ngrid
        DO k=ku1,ku2
           DO j=ju1,ju2
              !density initialization
              uin(l,i,j,k,1) = rho
              !rho*vx
              uin(l,i,j,k,2) = -Mup
              !rho*vy
              uin(l,i,j,k,3) = 0.d0
              !rho*vz
              uin(l,i,j,k,4) = 0.d0
              IF (pert .EQ. 'rhogrid') THEN
                 ! add random density perturbations in the heating layer
                 IF ((x(i) .GT. -1.d0).AND.(x(i) .LT. 1.d0)) THEN
                    call ran2(iseed,rvalue)
                    uin(l,i,j,k,1) = uin(l,i,j,k,1) + amprand*(rvalue-0.5)
                 END IF
              ENDIF
              ! add the perturbations :
              IF (pert .EQ. 'rho') THEN
                 ! add random density perturbations in the heating layer
                 IF ((x(i) .GT. -1.d0).AND.(x(i) .LT. 1.d0)) THEN
                    DO irand = nxmin,nxmax
                       DO jrand = nymin,nymax
                          uin(l,i,j,k,1) = uin(l,i,j,k,1) + fourrier_amp(irand,jrand)*COS(irand*(2.d0*pi*x(i)/2.d0 + fourrier_phasex(irand,jrand)))*COS(jrand*(2.d0*pi*y(j)/(ymax-ymin) + fourrier_phasey(irand,jrand)))
                       END DO
                    END DO
                 END IF
              ENDIF
              IF (pert .EQ. 'vgrid') THEN
                 ! add random velocities in the heating layer
                 IF ((x(i) .GT. -1.d0).AND.(x(i) .LT. 1.d0)) THEN
                    call ran2(iseed,rvalue)
                    uin(l,i,j,k,2) = uin(l,i,j,k,2)+ uin(l,i,j,k,1)*amprand*(rvalue-0.5)
                    IF (ndim .GE. 2) THEN
                       call ran2(iseed,rvalue)
                       uin(l,i,j,k,3) = uin(l,i,j,k,1)*amprand*(rvalue-0.5)
                    ENDIF
                    IF (ndim .EQ. 3) THEN
                       call ran2(iseed,rvalue)
                       uin(l,i,j,k,4) = uin(l,i,j,k,1)*amprand*(rvalue-0.5)
                    END IF
                 END IF
              ENDIF
              !bx
              uin(l,i,j,k,6) = B0           
              !by
              uin(l,i,j,k,7) = 0.d0
              !bz
              uin(l,i,j,k,8) = 0.d0
              
              !Energy
              uin(l,i,j,k,5) = rho*c**2/(gamma-1.)/gamma + &
                   & B0**2/2. + rho*v**2/2.d0
              
           END DO
        END DO
     END DO
  END DO
  

END SUBROUTINE condinit
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE give_gravity(x0, gravity)
  USE variables
  USE hydro_parameters
  USE heating_layer
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: x0
  DOUBLE PRECISION, iNTENT(OUT) :: gravity
  DOUBLE PRECISION :: shape

  IF (ABS(x0).LE. 0.5d0) THEN
     shape = 1.d0
  ENDIF
  IF (x0.GE.0.5d0) THEN
     shape = -(x0-1.d0)/0.5d0
  ENDIF
  IF (x0.LE.-0.5d0) THEN
     shape = (x0+1.d0)/0.5d0
  ENDIF
  IF (ABS(x0).GT.1.d0) THEN
     shape = 0.d0
  ENDIF
  
  gravity = -shape*Kgrav

  RETURN
END SUBROUTINE  give_gravity
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE give_cool(x0,rho,P,cool)
  USE variables
  USE hydro_parameters
  USE heating_layer
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: x0, rho, P
  DOUBLE PRECISION, iNTENT(OUT) :: cool
  DOUBLE PRECISION :: shape

  IF (ABS(x0).LE. 0.5d0) THEN
     shape = 1.d0
  ENDIF
  IF (x0.GE.0.5d0) THEN
     shape = -(x0-1.d0)/0.5d0
  ENDIF
  IF (x0.LE.-0.5d0) THEN
     shape = (x0+1.d0)/0.5d0
  ENDIF
  IF (ABS(x0).GT.1.d0) THEN
     shape = 0.d0
  ENDIF

  cool = shape*Kheat*rho*Mup/gamma

  RETURN
END SUBROUTINE  give_cool
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE derive0(x0,y0,dxy0)
 USE variables
  USE hydro_parameters
  USE heating_layer
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)   :: x0
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN)   :: Y0
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(INOUT) :: dxy0  
  ! y(1) : v, velocity of the flow (NB : v< 0)
  ! y(2) : c, sound speed
  ! x : coordinate along the x axis

  !-----------variables locales---------	
  DOUBLE PRECISION :: rhov, rho, P, grav, cool
  
  CALL give_gravity(x0, grav)
  rhov = -Mup
  rho = rhov/y0(1)
  P = rho*y0(2)**2/gamma
  CALL give_cool(x0,rho,P,cool)
  

  dxy0(1) = y0(1)/(y0(2)**2-y0(1)**2)*( -grav + (gamma-1.d0)*cool/rhov )
  dxy0(2) = (gamma - 1)/2.d0*y0(2)/(y0(2)**2-y0(1)**2)*( grav + (1.d0-gamma*y0(1)**2/y0(2)**2)*cool/rhov )
  	
  RETURN 
END SUBROUTINE derive0
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE rkutta0(x0,y0,dx0,ncount)
  USE hydro_parameters
  REAL(DP)   :: y0(2)
  REAL(DP)   ::x0,dx0
  INTEGER*4	ncount,nvar0
  REAL(DP)   ::a(2,2),b(2),c(2)
  !--------------------variables locales--------------	
  INTEGER*4 	i,n,icount,j
  REAL(DP) :: yst(2),dxy(2)
  REAL(DP) :: k(2,2),knew(2,2)	
  REAL(DP) :: dk(2,2),sup,xst(2)
  !-----------------------------------------------
  
  CALL initrk(a,b,c)
  nst = 2
  nvar0 = 2
  
  DO i = 1,nst
     xst(i) = x0+c(i)*dx0
     DO n = 1,nvar0
        k(i,n) = 0.d0
     ENDDO
  ENDDO
  
  icount = 0
3 sup = 0.d0
  DO i = 1,nst
     DO n = 1,nvar0
        yst(n) = y0(n)
     ENDDO
     DO n = 1,nvar0
        DO j = 1,nst
           yst(n) = yst(n)+a(i,j)*k(j,n)
        ENDDO
     ENDDO
     CALL derive0(xst(i),yst,dxy)
     
     DO n = 1,nvar0
        knew(i,n) = dx0*dxy(n)
        dk(i,n)   = ABS(knew(i,n)-k(i,n))
        k(i,n)    = knew(i,n)
        IF (dk(i,n).GT.sup) sup = dk(i,n)
     ENDDO
  ENDDO
  
  icount = icount+1
  IF ((sup.GE.1.d-10).AND.(icount.LT.maxcount)) GOTO 3
  
  ncount = icount
  DO n = 1,nvar0
     DO i = 1,nst
        y0(n) = y0(n)+b(i)*k(i,n)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE rkutta0
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE initRK(a,b,c)
  USE hydro_parameters
  REAL(DP) :: a(2,2),b(2),c(2)
  !------variables locales------------------	
  REAL(DP) :: p,q,r,t			
  !-------------------------------------	
  
  
  ! coeffs. de la methode implicite d'ordre 4,nst = 2 
  a(1,1) = 0.25d0
  a(1,2) = 0.25d0-dsqrt(3.d0)/6.d0
  a(2,1) = 0.25d0+dsqrt(3.d0)/6.d0
  a(2,2) = 0.25d0
  
  b(1) = 0.5d0
  b(2) = 0.5d0
  
  c(1) = 0.5d0-dsqrt(3.d0)/6.d0
  c(2) = 0.5d0+dsqrt(3.d0)/6.d0
  
  
  RETURN
END SUBROUTINE initRK
!====================================================================
!====================================================================
!====================================================================
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
      parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
               & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
               & IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      integer idum2,jj,kk,iv(NTAB),iy
      data idum2/123456789/, iv/NTAB*0/, iy/0/
!
      idum=iseed
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 jj=NTAB+8,1,-1
          kk=idum/IQ1
          idum=IA1*(idum-kk*IQ1)-kk*IR1
          if (idum.lt.0) idum=idum+IM1
          if (jj.le.NTAB) iv(jj)=idum
11      continue
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
!====================================================================
!====================================================================
!====================================================================

