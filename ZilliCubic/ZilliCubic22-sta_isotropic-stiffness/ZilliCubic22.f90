!“I don't know what the language of the year 2000 will look like, but I know it will be called Fortran.” —Tony Hoare, winner of the 1980 Turing Award, in 1982.
! From an article by Lee Philips at arstechnica.com .
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!   frc modfd2 duffingBN_forced modfd2 ZilliCubic:      
!   ZILLI CUBIC STA (Independent stiffness)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION ZETA,R2,Q1,Q3,Q2,Q4,GAMMA,OMEG,OMEGP,MH,EPSH,JPH,X,Y
      
        GAMMA = PAR(1) 
        OMEG = PAR(2)
        MH = PAR(3) 
        JPH  = PAR(4)
        EPSH = PAR(5)
        ZETA = PAR(6) 
        OMEGP = PAR(7)
        ! KAPPA = PAR(8) !|Nonlinearity Homotopy. 0:cubic, 1:contact
        ! BETA = PAR(10) !|Snubber ring stiffness measure: ks/kr 
        ! K = PAR(13) 
        ! KSI = PAR(14)
        ! RHO = PAR(15)

        Q1=U(1)
        Q2=U(2)
        Q3=U(3)
        Q4=U(4)
        X =U(5)
        Y =U(6)

        R2 = Q1**2 + Q2**2

        ! 'ISOTROPIC' CUBIC STIFFNESS X-EOM
        F(1)=Q3
        F(2)=Q4
        F(3)=-JPH*OMEG*Q4 &
             -2*ZETA*Q3 &
             -Q1 & 
             +MH*EPSH*( OMEGP*Y-(OMEG**2)*X ) &
             -GAMMA*R2*Q1   ! Cubic stiffness term
        F(4)=+JPH*OMEG*Q3 &
             -2*ZETA*Q4 &
             -Q2 &
             +MH*EPSH*( OMEGP*X+(OMEG**2)*Y ) &
             -GAMMA*R2*Q2   ! Cubic stiffness term
        F(5)=+X+OMEG*Y-X*(X**2 + Y**2) !|To remember the nonlinear oscillator   
        F(6)=+Y-OMEG*X-Y*(X**2 + Y**2) !|...remember 1st one, 2nd swaps x and y, and makes Omeg negative. 

        !| Sta-EOM from func_ode45.m > f_sta function 
        ! dq1 = (q3 ) ;
        ! dq2 = (q4 ) ;
        ! dq3 = (- JpH * (Omeg) * q4 ...
        !     - 2 * (zeta) * (q3) ...
        !     - q1 ... 
        !     + mH * epsH * ( OmegP * cos(phi) - (Omeg)^2 * sin(phi) ) ...
        !     - gamma * r2 * q1 ) ;%|ISO Stiffness one.
        ! dq4 = (+ JpH * (Omeg) * (q3) ...
        !     - 2 * (zeta) * (q4) ...
        !     - q2 ...
        !     + mH * epsH * ( OmegP * sin(phi) + (Omeg)^2 * cos(phi) ) ...
        !     - gamma * r2 * q2 ) ;%|ISO Stiffness one.





      ! IF (IJAC.EQ.1) RETURN
      !   DFDU(1,1)=0
      !   DFDU(1,2)=1
      !   DFDU(1,3)=0
      !   DFDU(1,4)=0
      !   DFDU(1,5)=0
      !   DFDU(1,6)=0

      !   DFDU(2,1)=-1-3*GAMMA*Q1**2
      !   DFDU(2,2)=-2*ZETA
      !   DFDU(2,3)=0
      !   DFDU(2,4)=-JPH*OMEG
      !   DFDU(2,5)=-MH*EPSH*(OMEG**2) 
      !   DFDU(2,6)=+MH*EPSH*OMEGP

      !   DFDU(3,1)=0
      !   DFDU(3,2)=0
      !   DFDU(3,3)=0
      !   DFDU(3,4)=1
      !   DFDU(3,5)=0
      !   DFDU(3,6)=0

      !   DFDU(4,1)=0
      !   DFDU(4,2)=+JPH*OMEG
      !   DFDU(4,3)=-1-3*GAMMA*Q2**2
      !   DFDU(4,4)=-2*ZETA
      !   DFDU(4,5)=+MH*EPSH*OMEGP
      !   DFDU(4,6)=+MH*EPSH*(OMEG**2)

      !   DFDU(5,1)=0
      !   DFDU(5,2)=0
      !   DFDU(5,3)=0
      !   DFDU(5,4)=0
      !   DFDU(5,5)=1- ( 1*(X**2 + Y**2) + X*(2*X) )    
      !   DFDU(5,6)=OMEG-X*(2*Y)    

      !   DFDU(6,1)=0
      !   DFDU(6,2)=0
      !   DFDU(6,3)=0
      !   DFDU(6,4)=0
      !   DFDU(6,5)=-OMEG-Y*(2*X)
      !   DFDU(6,6)=1- ( 1*(X**2 + Y**2) + Y*(2*Y))
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION ZETA,Q1,Q3,Q2,Q4,GAMMA,THETA,OMEG,OMEGP,MH,EPSH,JPH,X,Y,PI

        GAMMA = 0.25d0 !2nd continue gamma from 0 to 0.25 Or start it with 0.25 already
        OMEG = 0.01d0 !|IMPORTANT changing 0.0d0 to 0.01d0 solved the "0 of 1 lapack etc" problem.
        MH = 0.9d0 !1st continue MH from 0 to 0.9
        
        EPSH =0.d0 !1d-10 !0.353d0 !first homotopy is for this.
        ZETA = 0.0d0 !0.1d0 !0.01
        JPH  = 0.0143d0 !|for a TRIAL chaged from 0.143d0 to 0.0143d0
        OMEGP = 0.d0
 
        PI = 4.D0*ATAN(1.0D0) !! = pi 
        PAR(11) = 2 * ( 2.D0*PI/OMEG )   !! PAR(11) is always reserved for period in AUTO. 
        !|:Period = 2*pi/OMEG
        !|::In the synchronous whirl the period of oscillation is the same as the rotor speed!!
        !|::In the rotating frame synch freq is 0; so here we have: w_sta=w_rot+Omeg = Omeg. 
        !|:In the sta-frame subharmonic response of order 1/2, the freq content has OMEG and OMEG/2. 
        !|::So the EBOB of the freqs is OMEG/2. Therefore, the period is 2*(2.D0*PI/OMEG) !BN
        !|:Because we define the PAR(11) here explicitly, the time variable becomes 
        !|...scaled to time/period=T in range 0-1.

        X = SIN(2.D0*PI*T) 
        Y = COS(2.D0*PI*T)
        
        ! THETA = OMEG*(TPI*T)
        THETA = ATAN2(X,Y) !|Position of eccentric mass (origY at x-axis) 
        !|:ATAN2(Y,X) angle of X+jY in correct quadrant
        !|:We have X=sin ,Y=cos; so  ATAN2(X,Y) should be correct !BN  

 
        PAR(1)=GAMMA 
        PAR(2)=OMEG 
        PAR(3)=MH 
        PAR(4)=JPH  
        PAR(5)=EPSH 
        PAR(6)=ZETA 
        PAR(7)=OMEGP 
        ! PAR(8)=THETA

        U(1)=0.d0
        U(2)=0.d0
        U(3)=0.d0
        U(4)=0.d0
        U(5)=X
        U(6)=Y
      END SUBROUTINE STPNT


      SUBROUTINE PVLS(NDIM,U,PAR)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM) !|Dsll take only 1st col of U (:1st datum f d "solution")
        DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
        DOUBLE PRECISION R2,Q1,Q2
        DOUBLE PRECISION, EXTERNAL :: GETP, GETU
        DOUBLE PRECISION :: DUMM
        INTEGER :: NDX,NCOL,NTST,i,j,k
        !|SEE the explanation of demo pvl.f90 : 
        !|:"For algebraic problems the argument U is, as usual, the state vector.
        !|...For differential equations the argument U represents the approximate 
        !|...solution on the entire interval [0,1]. In this case its values can
        !|...be accessed indirectly by calls to GETP, as illustrated below, or
        !|...by obtaining NDIM, NCOL, NTST via GETP and then dimensioning U as
        !|...U(NDIM,0:NCOL*NTST) in a seperate subroutine that is called by PVLS. "
        !|::See Keep #fortran "fortran assumed dimensions and how it is interpreted."

        NDX  = NINT(GETP('NDX',0,U))
        NTST = NINT(GETP('NTST',0,U))
        NCOL = NINT(GETP('NCOL',0,U))

        !|___METHOD-A/B> METHOD A/B  
        !|:CORRECT AMPLITUDE FOR IPS=2
        DUMM = GETU(i,j,U,NDX,NTST,NCOL,PAR(2)) 
        !|:Will get amp and print the states to dummy_write_file.txt
        PAR(9) = DUMM 
        !|___METHOD-A/B.
        
        !|___METHOD-C> METHOD C
        ! DO i=1,NDX
        !   DO j=0,NCOL*NTST
        !     UU(i,j) = GETU(i,j,U,NDX,NTST,NCOL)
        !   END DO
        ! END DO
        ! R2(:) = UU(1,:)**2 + UU(2,:)**2  !**El by El multip
        ! R(:) = R2**0.5 
        ! PAR(9) = MAXVAL(R)
        !|___METHOD-C.

        !|___METHOD-D> DEFAULT METHOD - PHASED AMP
        ! PAR(9) = (U(1)**2 + U(2)**2)**0.5 
        !|:The initial data of one period is sampled for amplitude.
        !|:See the qoute above.
        !|___METHOD-D. 
        
      END SUBROUTINE PVLS




      DOUBLE PRECISION FUNCTION GETU(i,j,U,NDX,NTST,NCOL,MYPAR)
        INTEGER, INTENT(IN) :: NDX,NCOL,NTST,i,j
        INTEGER :: p, I_MAX
        DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
        DOUBLE PRECISION, ALLOCATABLE :: R2(:),R(:)
        DOUBLE PRECISION :: M,N , MYPAR
        !| MYPAR - the parameter to print in the 1st col of dummy_write_file.txt

        !|___METHOD-A>
        M = 0.D0
        DO p=0,NCOL*NTST
          N = ( U(1,p)**2+U(2,p)**2 )**0.5D0
          IF (N.GT.M) THEN
            M = N
            I_MAX = p
          END IF 
        END DO 

        !| Print the amplitude point's all 4 states.
        !| ...E.g for time simulation initial data, 
        !| ...parse the following file.  
        ! OPEN (unit = 10, file = "dummy_write_file.txt",access = "append")
        ! WRITE(10,*) MYPAR,U(1,I_MAX),U(2,I_MAX),U(3,I_MAX),U(4,I_MAX),M
        ! CLOSE(10)

        GETU = M !|return the maximum amplitude.
        
        !|___METHOD-A.

        !|___METHOD-B> This is forced to handle large orbit data: Unnecessary, use A. 
        ! R2(:) = U(1,:)**2 + U(2,:)**2  !
        ! R(:) = R2**0.5 
        ! GETU = MAXVAL(R)
        !|___METHOD-B.

        !|___METHOD-C> Get the pure U, find the max in the PVLS subroutine.
        ! GETU = U(i,j) 
        !|___METHOD-C.
      END FUNCTION GETU












      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

