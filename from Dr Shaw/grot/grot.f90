!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!                   grot:      SDOF gravity bearing rotor stator system
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION rr, rs, muc, Omega, THETA, THETADOT, gor ,v ,sgnv, mubrk, cv, ff, mu

	  mubrk=0.2
	  cv=0.3
	  ff=0.005
	  
       rr=PAR(1)
	   rs=PAR(2)
	   muc=PAR(3)
	   Omega=PAR(4)
	   
	   ! state terms
	   THETA=U(1)  ! BN this is position of the rotor wrt stator
	   THETADOT=U(2)
	   
	   !some preliminary calcs
	   gor=9.81/(rs-rr)
	   v=Omega*rr-THETADOT
	   sgnv=ATAN(50.0*v)/(2.0*ATAN(1.0))
	   
	   ! set time derivatives  
	   F(1)= THETADOT
       !F(2)= muc*sgnv*(THETADOT**2 + gor*COS(THETA)) -gor*SIN(THETA)
	   ! using stribeck params
	   mu=(muc+(mubrk-muc)*EXP(-1.0*CV*ABS(v)))*sgnv + ff*v
	   F(2)= mu*(THETADOT**2 + gor*COS(THETA)) -gor*SIN(THETA)

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

	  DOUBLE PRECISION rr, rs, muc, Omega, TPI, gor ,v ,sgnv, mubrk, cv, ff, mu

	  mubrk=0.2
	  cv=0.3
	  ff=0.005
	  
       rr=1.0
	   rs=1.05
	   muc=0.15
	   Omega=15.0
	   
	   PAR(1)=rr
       PAR(2)=rs
       PAR(3)=muc
       PAR(4)=Omega
       TPI=8*ATAN(1.0)
       PAR(11)=0.0
	   
	   v=Omega*rr
	   sgnv=ATAN(500.0*v)/(2.0*ATAN(1.0))
	   
       !U(1)=ATAN(muc*sgnv)
       ! using stribeck params
	   mu=(muc+(mubrk-muc)*EXP(-1.0*CV*ABS(v)))*sgnv + ff*v
	   U(1)=ATAN(mu)
	   U(2)=0.0
       
      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
