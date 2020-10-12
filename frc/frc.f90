!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!   frc :      A periodically forced system
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION FHN,Z,A,B,R,EPS,BET,D,V,W,X,Y,SS

       FHN(Z) = Z*(Z-A)*(1-Z) 

       A  = PAR(1) 
       B  = PAR(2) 
       R  = PAR(3) 
       EPS= PAR(4) 
       BET= PAR(5) 
       D  = PAR(6)

       V=U(1)
       W=U(2)
       X=U(3)
       Y=U(4)
       SS = X**2 + Y**2    

       F(1) = ( FHN(V) - W )/EPS
       F(2) = V - D*W - ( B + R*X )
       F(3) =  X + BET*Y - X*SS    
       F(4) = -BET*X + Y - Y*SS    

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      DOUBLE PRECISION A,B,R,EPS,BET,D,TPI

       A  =0.5
       B  = A
       R  =0.
       EPS=0.005
       BET=100  !! BN Why is this given a high number?
       D  =1.0

       PAR(1)=A
       PAR(2)=B
       PAR(3)=R  !! BN Start forcing from zero amplitude. 
            !! 1st run "continues" with ICP=[3,11]: ie from zero forcing to some forcing amplitude of r=0.5
            !! During this run the BET is kept constant. 
            !! DESPITE, WHY IS c.frc FILE HAS SET USER OUTPUT VALUES FOR PAR11 ??? 
       PAR(4)=EPS
       PAR(5)=BET 
            !! BN 2nd run continues this beta parameter. So ICP=[5,11],
            !! BN Beta is now changing, BUT WHY SET CERTAIN USER OUTPUT VALUES FOR PAR11, WHILE IT IS CONTROLLED EXACTLY BY BETA ALREADY ??? 
       PAR(6)=D
       TPI=8*ATAN(1.0D0)  !! = 2*pi
       PAR(11)=TPI/BET    !! = 2pi/omeg
            !! PAR(11) is always reserved for period in AUTO.

       U(1)=B              !! BN When forcing is zero(r=0), stationary sol of U1 isDs.
       U(2)=B*(B-A)*(1-B)  !! BN When forcing is zero(r=0), stationary sol of U2 isDs.
       U(3)=SIN(TPI*T)     !! BN Starting point (STPNT) is sin(2*pi*T), why?: to start form zero angle BN
       U(4)=COS(TPI*T)     !! BN Starting point (STPNT) is cos(2*pi*T), why?: to start from zero angle BN

      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
