!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!                   dufa:      forced duffing oscillator, attempt a 
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION kcube, zeta, P, Omega, Q, QDOT, SS, X, Y, TPI

       kcube=PAR(1)
	   zeta=PAR(2)
	   P=PAR(3)
	   Omega=PAR(4)
	   
	   ! Duffing oscillator state terms
	   Q=U(1)
	   QDOT=U(2)
	   ! These state terms give periodic forcing
	   X = U(3)
	   Y=U(4)
	   SS=X**2 + Y**2
       
	   ! Duffing mechanism
	   F(1)= QDOT
       F(2)= -2*zeta*QDOT-Q-kcube*Q**3 + P*Y
	   
	   ! oscillator mechanism
	   F(3)= X+Omega*Y-X*SS
	   F(4)= -Omega*X+Y-Y*SS
    
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

	  DOUBLE PRECISION kcube, zeta, P, Omega,TPI

	  
       kcube=0.2
	   zeta=0.05
	   P=0.0
	   Omega=1.
	   
	   PAR(1)=kcube
       PAR(2)=zeta
       PAR(3)=P
       PAR(4)=Omega
       TPI=8*ATAN(1.0D0)
       PAR(11)=TPI/Omega
	   
       U(1)=0.
       U(2)=0.
	   U(3)=SIN(TPI*Omega*T)
       U(4)=COS(TPI*Omega*T)

       
      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
