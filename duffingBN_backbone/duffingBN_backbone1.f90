! An attemp to simply eliminate the 3rd and 4th equations (the extra oscillators to track the response frequency).
SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,*), DFDP(NDIM,*)
  DOUBLE PRECISION WN,K3_M,U1,U2,X,Y,OMEG 

  WN   = PAR(1) 
  K3_M = PAR(2) 
  OMEG = PAR(3)

  U1 = U(1)
  U2 = U(2)
  X  = U(3)
  Y  = U(4)

  F(1) =  U2
  F(2) = -WN**2*U1 - K3_M*U1**3 !- 2*ZETA*WN*U2 + F_M*X
  ! F(3) =  X + OMEG*Y - X*(X**2 + Y**2)    
  ! F(4) = -OMEG*X + Y - Y*(X**2 + Y**2)
END SUBROUTINE FUNC


SUBROUTINE STPNT(NDIM,U,PAR,T) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T
  DOUBLE PRECISION WN,K3_M,U1,U2,TPI,X,Y,OMEG 
  
  WN   = 1.0
  K3_M = 1.0 
  OMEG = 3.0

  PAR(1) = WN 
  PAR(2) = K3_M 
  PAR(3) = OMEG

  TPI=8*ATAN(1.0D0) ! = 2*pi 
  PAR(11)= 0.0 !|IPS=1 is selected !|=TPI/OMEG  ; Period = 2*pi/omeg 
  ! X = SIN(TPI*T) 
  ! Y = COS(TPI*T)

  U(1) = 0
  U(2) = 0
  ! U(3) = X
  ! U(4) = Y
END SUBROUTINE STPNT


SUBROUTINE PVLS(NDIM,U,PAR)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

  ! PAR(9) = U(1)**2 + U(2)**2
END SUBROUTINE PVLS


SUBROUTINE BCND 
END SUBROUTINE BCND

SUBROUTINE ICND 
END SUBROUTINE ICND

SUBROUTINE FOPT 
END SUBROUTINE FOPT

