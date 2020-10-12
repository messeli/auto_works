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

      DOUBLE PRECISION ZETA,U1,U2,U3,U4,U5,S4,C4,R2,GAMMA,THETAP,THETAPP,MH,EPSH,JPH,OMEG,THETA !WN,F_M,K3_M,
      
       GAMMA = PAR(1) 
       THETAP = PAR(2)
       MH = PAR(3) 
       EPSH = PAR(4)
       ZETA = PAR(5) 
       JPH  = PAR(6)
       THETAPP = PAR(7)
       
       U1=U(1)
       U2=U(2)
       U3=U(3)
       U4=U(4)
       ! THETA=U(5) !ROTATION ANGLE OF THE ROTOR
       ! X=U(6)
       ! Y=U(7)

	   ! S4 = SIN(4*THETA) 
	   ! C4 = COS(4*THETA)

	   ! R2 = U1**2+U3**2 ;
       
       F(1)=U2
       F(2)=-THETAP*(JPH-2)*U4-2*ZETA*U2+(THETAP**2*(1-JPH)-1)*U1+2*ZETA*THETAP*U3+MH*EPSH*THETAP**2 !-GAMMA/4*(3*U1*R2+C4*U1*(R2-4*U3^2)+S4*U3*(R2-4*U1^2)) 
       F(3)=U4
       F(4)=+THETAP*(JPH-2)*U2-2*ZETA*U4+(THETAP**2*(1-JPH)-1)*U3-2*ZETA*THETAP*U1-MH*EPSH*THETAPP  !-GAMMA/4*(3*U3*R2+C4*U3*(R2-4*U1^2)-S4*U1*(R2-4*U3^2)) 
       ! F(5)=THETAP 
       ! F(6) =  X + OMEG*Y - X*(X**2 + Y**2)    
       ! F(7) = -OMEG*X + Y - Y*(X**2 + Y**2)

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      DOUBLE PRECISION ZETA,U1,U2,U3,U4 ,S4,C4,R2 ,GAMMA,THETAP,THETAPP,MH,EPSH,JPH !WN,F_M,K3_M,OMEG,

       GAMMA = 0 !0.25 Continue gamma from 0 to 0.25
       THETAP = 0 ! Continue THETAP from 0 to 7. 
       	! Find stat sols for thetaP=0; put into U() below. Then continue by incrs thetaP.
       MH = 0.9 !Continue MH from 0 to 0.9
       EPSH =0.353
       ZETA = 0.01 !0.01
       JPH  = 0.143
       THETAPP = 0

! Initialize the equation parameters
       PAR(1)=GAMMA 	
       PAR(2)=THETAP 
       PAR(3)=MH 
       PAR(4)=EPSH 
       PAR(5)=ZETA 
       PAR(6)=JPH  
       PAR(7)=THETAPP 
 

       ! TPI=8*ATAN(1.0D0) !! = 2*pi 
       ! PAR(11)=TPI/OMEG   !! Period = 2*pi/omeg
            !! PAR(11) is always reserved for period in AUTO.

! Initialize the solution
       U(1)=0
       U(2)=0
       U(3)=0
       U(4)=0
       ! U(5)=0
       ! U(6)=SIN(TPI*T)
       ! U(7)=COS(TPI*T)

      END SUBROUTINE STPNT


      SUBROUTINE PVLS(NDIM,U,PAR)

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: NDIM
       DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
       DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      
       DOUBLE PRECISION R2

        R2 = U(1)**2+U(3)**2 

        PAR(8) = R2**(0.5)

      END SUBROUTINE PVLS


      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

