!“I don't know what the language of the year 2000 will look like, but I know it will be called Fortran.” —Tony Hoare, winner of the 1980 Turing Award, in 1982.
! From an article by Lee Philips at arstechnica.com .
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!   frc modfd2 duffingBN_forced modfd2 ZilliCubic:      
!   ZILLI CUBIC ROT (NOT AUTONOMOUS - INDEPENDENT STIFFNESS)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION ZETA,Q1,Q2,Q3,Q4,U5,S4,C4,R2,GAMMA,THETA,THETAP,THETAPP,MH,EPSH,JPH,X,Y
      
        GAMMA = PAR(1) 
        THETAP = PAR(2)
        MH = PAR(3) 
 
        EPSH = PAR(4)
        ZETA = PAR(5) 
        JPH  = PAR(6)
        THETAPP = PAR(7)
        ! THETA = PAR(8)


        Q1=U(1)
        Q2=U(2)
        Q3=U(3)
        Q4=U(4)
        X =U(5)
        Y =U(6)
        ! THETA=U(5) !ROTATION ANGLE OF THE ROTOR
        
        
        ! THETA = ATAN2(X,Y)
        S4 = X !SIN(4*THETA) 
        C4 = Y !COS(4*THETA) 
 
        R2 = (Q1**2)+(Q3**2)
  
 
        F(1) = Q2
        F(2) = -THETAP*(JPH-2)*Q4          &
               -2*ZETA*Q2                  &
               +((THETAP**2)*(1-JPH)-1)*Q1 &
               +2*ZETA*THETAP*Q3           &
               +MH*EPSH*(THETAP**2)        &
               -GAMMA/4*( 3*Q1*R2 + C4*Q1*(R2-4*(Q3**2)) + S4*Q3*(R2-4*(Q1**2)) ) 
        F(3) = Q4
        F(4) = +THETAP*(JPH-2)*Q2          &
               -2*ZETA*Q4                  &
               +((THETAP**2)*(1-JPH)-1)*Q3 &
               -2*ZETA*THETAP*Q1           &
               -MH*EPSH*THETAPP            &
               -GAMMA/4*( 3*Q3*R2 + C4*Q3*(R2-4*(Q1**2)) - S4*Q1*(R2-4*(Q3**2)) ) 
        F(5) = +X+(4*THETAP)*Y-X*(X**2+Y**2)    
        F(6) = -(4*THETAP)*X+Y-Y*(X**2+Y**2)

      IF (IJAC.EQ.1) RETURN
        DFDU(1,1)=0
        DFDU(1,2)=1
        DFDU(1,3)=0
        DFDU(1,4)=0
        DFDU(1,5)=0
        DFDU(1,6)=0

        DFDU(2,1)=+((THETAP**2)*(1-JPH)-1) &
                  -GAMMA/4*( 3*R2 + 6*Q1**2 + C4*(R2-4*Q3**2) + 2*C4*Q1**2 - 6*S4*Q3*Q1 )
        DFDU(2,2)=-2*ZETA
        DFDU(2,3)=+2*ZETA*THETAP + &
                  -GAMMA/4*( +6*Q1*Q3 - 6*C4*Q1*Q3 + S4*(R2-4*Q1**2) + 2*S4*Q3**2  )
        DFDU(2,4)=-THETAP*(JPH-2)
        DFDU(2,5)=-GAMMA/4*Q3*(R2-4*Q1**2) 
        DFDU(2,6)=-GAMMA/4*Q1*(R2-4*Q3**2)  

        DFDU(3,1)=0
        DFDU(3,2)=0
        DFDU(3,3)=0
        DFDU(3,4)=1
        DFDU(3,5)=0
        DFDU(3,6)=0

        DFDU(4,1)=-2*ZETA*THETAP &
                  -GAMMA/4*( 6*Q3*Q1 - 6*C4*Q3*Q1 - S4*(R2-4*Q3**2) - 2*S4*Q1**2 )
        DFDU(4,2)=+THETAP*(JPH-2)
        DFDU(4,3)=+((THETAP**2)*(1-JPH)-1) &
                  -GAMMA/4*( 3*R2 + 6*Q3**2 + C4*(R2-4*Q1**2) + 2*C4*Q3**2 +6*S4*Q1*Q3 )
        DFDU(4,4)=-2*ZETA
        DFDU(4,5)=+GAMMA/4*Q1*(R2-4*Q3**2)  
        DFDU(4,6)=-GAMMA/4*Q3*(R2-4*Q1**2) 

        DFDU(5,1)=0
        DFDU(5,2)=0
        DFDU(5,3)=0
        DFDU(5,4)=0
        DFDU(5,5)=1- ( 1*(X**2 + Y**2) + X*(2*X) )    
        DFDU(5,6)=4*THETAP-X*(2*Y)    

        DFDU(6,1)=0
        DFDU(6,2)=0
        DFDU(6,3)=0
        DFDU(6,4)=0
        DFDU(6,5)=-4*THETAP-Y*(2*X)
        DFDU(6,6)=1- ( 1*(X**2 + Y**2) + Y*(2*Y))
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      DOUBLE PRECISION ZETA,Q1,Q2,Q3,Q4,S4,C4,R2,GAMMA,THETA,THETAP,THETAPP,MH,EPSH,JPH,X,Y,TPI!WN,F_M,K3_M,OMEG,

       GAMMA = 0.0 !2nd continue gamma from 0 to 0.25
       THETAP = 0.1 ! Continue THETAP from 0 to 7. 
       MH = 0.0 !1st continue MH from 0 to 0.9
       
       EPSH =0.353
       ZETA = 0.01 !0.07
       JPH  = 0.143
       THETAPP = 0.0

       TPI=8*ATAN(1.0D0) !! = 2*pi 
       PAR(11)=TPI/THETAP   !! Period = 2*pi/thetaP
            ! PAR(11) is always reserved for period in AUTO.
            ! Because we define the PAR(11) here explicitly, the T time variable becomes 
            ! scaled to time/period=T.
            ! THEREFORE: omeg*time =period*omeg*(time/period)=(2pi/omeg)*omeg*T = 2pi*T

       X = SIN(4*TPI*T) ! HERE THEN, sin(4theta)=sin(4thetaP*t)=sin(4*2pi*T) 
       Y = COS(4*TPI*T)
       
       ! THETA = THETAP*(TPI*T)
       THETA = ATAN2(X,Y)

       PAR(1)=GAMMA 
       PAR(2)=THETAP 
       PAR(3)=MH 
       PAR(4)=EPSH 
       PAR(5)=ZETA 
       PAR(6)=JPH  
       PAR(7)=THETAPP 
       ! PAR(8)=THETA
 
       U(1)=0
       U(2)=0
       U(3)=0
       U(4)=0
       U(5)= X
       U(6)= Y
      END SUBROUTINE STPNT


      SUBROUTINE PVLS(NDIM,U,PAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION R2,Q1,Q3
        Q1 = U(1)
        Q3 = U(3)

        R2 = Q1**2+Q3**2 
        ! Which one is this?
        ! ?? max( U(1) )^2 + max( U(3) )^2
        ! ?? max( U(1)^2 + U(3)^2 )

        PAR(9) = R2**(0.5)
      END SUBROUTINE PVLS


      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

