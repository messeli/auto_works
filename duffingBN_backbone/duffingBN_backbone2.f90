! An attemp to start from k3_m=0 and continue it with par(11), to eliminate F(3) and F(4) eqs, the extra oscillator.
SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,*), DFDP(NDIM,*)
  DOUBLE PRECISION WN,K3_M,X,Y

  WN   = PAR(1) 
  K3_M = PAR(2) 

  F(1) =  U(2)
  F(2) = -WN**2*U(1) - K3_M*U(1)**3 
END SUBROUTINE FUNC


SUBROUTINE STPNT(NDIM,U,PAR,T) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T
  DOUBLE PRECISION WN,K3_M,U1,U2,TPI,X,Y,WR 
  
  WN   = 1.0
  K3_M = 0.0 
  WR = 1.0

  PAR(1) = WN 
  PAR(2) = K3_M 

  TPI=8*ATAN(1.0D0) ! = 2*pi 
  PAR(11)= TPI/WR   ! Period = 2*pi/omeg

  U(1) = 2.0  * COS(TPI*T) ! 2*cos(2pi * t/period)
  U(2) =-2.0*WR*SIN(TPI*T) ! 2*cos(2pi * t/period)
END SUBROUTINE STPNT


SUBROUTINE PVLS(NDIM,U,PAR)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
  DOUBLE PRECISION WR,TPI

  TPI=8*ATAN(1.0D0) ! = 2*pi 
  WR = TPI/PAR(11)  ! Response freq 
  PAR(9) = WR
END SUBROUTINE PVLS


SUBROUTINE BCND 
END SUBROUTINE BCND

SUBROUTINE ICND 
END SUBROUTINE ICND

SUBROUTINE FOPT 
END SUBROUTINE FOPT
