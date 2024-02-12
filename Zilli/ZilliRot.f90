!“I don't know what the language of the year 2000 will look like, but I know it will be called Fortran.” —Tony Hoare, winner of the 1980 Turing Award, in 1982.
! From an article by Lee Philips at arstechnica.com .
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!   frc modfd2 duffingBN_forced modfd2 ZilliCubic:      
!   ZILLI ROT
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION ZETA,Q1,Q2,Q3,Q4,BETA,THETAP,THETAPP,MH,EPSH,JPH,R,ISCONTACT
      
        BETA = PAR(1) 
        THETAP = PAR(2)
        MH = PAR(3) 

        EPSH = PAR(4)
        ZETA = PAR(5) 
        JPH  = PAR(6)
        THETAPP = PAR(7)

        Q1=U(1)
        Q2=U(2)
        Q3=U(3)
        Q4=U(4)
 
        R = SQRT(Q1**2+Q3**2)

        !| SIGMOID FUNCTION to MODEL: isContact = R >= 1  
        !| tanh is closer to heaviside step function than sigmoid > See #20 presentation
        ISCONTACT = 1/( 1+EXP(-10*(R-1)) )

        F(1) = Q2
        F(2) = -THETAP*(JPH-2)*Q4          &
               -2*ZETA*Q2                  &
               +((THETAP**2)*(1-JPH)-1)*Q1 &
               +2*ZETA*THETAP*Q3           &
               +MH*EPSH*(THETAP**2)        &
               +ISCONTACT*(1/R-1)*BETA*Q1 
        F(3) = Q4
        F(4) = +THETAP*(JPH-2)*Q2          &
               -2*ZETA*Q4                  &
               +((THETAP**2)*(1-JPH)-1)*Q3 &
               -2*ZETA*THETAP*Q1           &
               -MH*EPSH*THETAPP            &
               +ISCONTACT*(1/R-1)*BETA*Q3 
       END SUBROUTINE FUNC


      SUBROUTINE STPNT(NDIM,U,PAR,T)  
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION ZETA,Q1,Q2,Q3,Q4,BETA,THETAP,THETAPP,MH,EPSH,JPH

        BETA = 0.0 !2nd, continue gamma from 0 to 0.25
        THETAP = 0.1 !3rd, continue THETAP from 0 to 7. 
        MH = 0.0 !1st, continue MH from 0 to 0.9
       
        EPSH = 0.353
        ZETA = 0.01 !0.1
        JPH  = 0.143
        THETAPP = 0.0

        PAR(1)=BETA 
        PAR(2)=THETAP 
        PAR(3)=MH 
        PAR(4)=EPSH 
        PAR(5)=ZETA 
        PAR(6)=JPH  
        PAR(7)=THETAPP 

        U(1)=0
        U(2)=0
        U(3)=0
        U(4)=0
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

        PAR(9) = SQRT(R2)
      END SUBROUTINE PVLS


      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

