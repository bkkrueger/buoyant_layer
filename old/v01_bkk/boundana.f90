SUBROUTINE xinner_ana
  USE hydro_parameters
  USE variables
  IMPLICIT NONE
!
  INTEGER :: j,k
!
  ! Zero-gradient
  DO k=ku1,ku2
     DO j=ju1,ju2      
        uin(1,iu1  ,j,k,:) = uin(1,iu1+3,j,k,:)
        uin(1,iu1+1,j,k,:) = uin(1,iu1+3,j,k,:)
        uin(1,iu1+2,j,k,:) = uin(1,iu1+3,j,k,:)
     END DO
  END DO
!
  RETURN
END SUBROUTINE xinner_ana
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE xouter_ana
  USE hydro_parameters
  USE variables
  IMPLICIT NONE

  INTEGER ::k,j

  ! Zero-gradient
  DO k=ku1,ku2
     DO j=ju1,ju2
        uin(1,iu2-2,j,k,:) = uin(1,iu2-3,j,k,:)
        uin(1,iu2-1,j,k,:) = uin(1,iu2-3,j,k,:)
        uin(1,iu2  ,j,k,:) = uin(1,iu2-3,j,k,:)
     END DO
  END DO
!
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
