!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!   frc modfd2 duffingBN_forced :      A periodically forced system
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION WN,ZETA,F_M,K3_M,OMEG,U1,U2,X,Y
      
       WN   = PAR(1) 
       ZETA = PAR(2) 
       F_M  = PAR(3)
       K3_M = PAR(4) 
       OMEG = PAR(5)

       U1=U(1)
       U2=U(2)
       X =U(3)
       Y =U(4)

       F(1) =  U2
       F(2) = -WN*U1 - 2*ZETA*WN*U2 - K3_M*U1**3 + F_M*X
       F(3) =  X + OMEG*Y - X*(X**2 + Y**2)    
       F(4) = -OMEG*X + Y - Y*(X**2 + Y**2)

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      DOUBLE PRECISION WN,ZETA,F_M,K3_M,OMEG,TPI

       WN   = 10 
       ZETA = 0.01 
       F_M  = 0.0
       K3_M = 0.05  !! BN 0.01 tried worked.
       OMEG = 1
      
! Initialize the equation parameters
       PAR(1) = WN 
       PAR(2) = ZETA
       PAR(3) = F_M 
       PAR(4) = K3_M 
       PAR(5) = OMEG

       TPI=8*ATAN(1.0D0) !! = 2*pi 
       PAR(11)=TPI/OMEG   !! Period = 2*pi/omeg
            !! PAR(11) is always reserved for period in AUTO.

! Initialize the solution
       U(1)=0
       U(2)=0
       U(3)=SIN(TPI*T)
       U(4)=COS(TPI*T)

      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
