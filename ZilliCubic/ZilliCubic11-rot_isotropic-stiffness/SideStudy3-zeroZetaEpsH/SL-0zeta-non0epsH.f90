!“I don't know what the language of the year 2000 will look like, but I know it will be called Fortran.” —Tony Hoare, winner of the 1980 Turing Award, in 1982.
! From an article by Lee Philips at arstechnica.com .
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!   frc modfd2 duffingBN_forced modfd2 ZilliCubic:      
!   ZILLI CUBIC ROT (AUTONOMOUS - ISOTROPIC STIFFNESS)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
        DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
        DOUBLE PRECISION ZETA,C,Q1,Q3,Q2,Q4,R2,R,GAMMA,BETA,KAPPA,RHO,K , X, Y, TPI
        DOUBLE PRECISION OMEG,OMEGP,MH,EPSH,JPH,FSNUB_u,FSNUB_v,FCUBIC_u,FCUBIC_v,F5_u,F5_v,KSI
      
        GAMMA = PAR(1) 
        OMEG = PAR(2)
        MH = PAR(3) 
        JPH = PAR(4) 
        EPSH = PAR(5)
        ZETA = PAR(6) 
        OMEGP = PAR(7)
        KAPPA = PAR(8) !|Nonlinearity Homotopy. 0:cubic, 1:contact
        BETA = PAR(10) !|Snubber ring stiffness measure: ks/kr 
        K = PAR(13) 
        KSI = PAR(14)
        RHO = PAR(15)

        Q1=U(1)
        Q2=U(2)
        Q3=U(3)
        Q4=U(4)
        X =U(5)
        Y =U(6)

        R2 = (Q1**2)+(Q2**2)
        R = R2**0.5D0
        C = 1 !|NondimD wrt c-clearance.
        TPI=8*ATAN(1.0D0) ! = 2*pi 

        FSNUB_u  = -0.5D0*(tanh(K*(R-1))+1) *(1-1/R)*BETA*Q1 
        FSNUB_v  = -0.5D0*(tanh(K*(R-1))+1) *(1-1/R)*BETA*Q2 
        FCUBIC_u = -GAMMA*R2*Q1
        FCUBIC_v = -GAMMA*R2*Q2
        F5_u     = -KSI*R2**2*Q1
        F5_v     = -KSI*R2**2*Q2

        F(1) = Q3
        F(2) = Q4
        F(3) = -OMEG*(JPH-2)*Q4          &
               -2*ZETA*Q3                &
               +((OMEG**2)*(1-JPH)-1)*Q1 &
               +2*ZETA*OMEG*Q2           &
               +MH*EPSH*(OMEG**2)        &
               +KAPPA*FSNUB_u + RHO*F5_u + (1-KAPPA-RHO)*FCUBIC_u !|NONLINEARITY HOMOTOPY

        F(4) = +OMEG*(JPH-2)*Q3          &
               -2*ZETA*Q4                &
               +((OMEG**2)*(1-JPH)-1)*Q2 &
               -2*ZETA*OMEG*Q1           &
               -MH*EPSH*OMEGP            &
               +KAPPA*FSNUB_v + RHO*F5_v + (1-KAPPA-RHO)*FCUBIC_v !|NONLINEARITY HOMOTOPY

        F(5) = X + TPI/PAR(11)*Y - X*(X**2 + Y**2) !|Omeg is replaced by 2pi/per, as the per of osc is not Omeg in this problem.
        F(6) = Y - TPI/PAR(11)*X - Y*(X**2 + Y**2) !|Omeg is replaced by 2pi/per, as the per of osc is not Omeg in this problem.

      ! IF (IJAC.EQ.1) RETURN
      !   DFDU(1,1)=0
      ! ... 
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: T
        DOUBLE PRECISION ZETA,Q1,Q3,Q2,Q4,S4,C4,R2,GAMMA,KAPPA,K,OMEG
        DOUBLE PRECISION OMEGP,MH,EPSH,JPH,BETA,TPI,KSI,RHO, X,Y, TPI

        !| .DAT file data from MATLAB has gamma Omeg mh epsh zeta JpH OmegP as below:
        GAMMA = 0.25D0 !2nd continue gamma from 0 to 0.25 !|0-0.25: Cubic stiffness ratio, (k_3*r^2)/k_r
        OMEG  = 3.50D0!3.50D0 3.46D0 3.40D0 3.32D0 !The speed the .dat file is generated at, in Matlab ode45
        MH    = 0.9D0 !1st continue MH from 0 to 0.9
        EPSH  = 0.353D0 !was 0.353D0
        ZETA  = 0.0D0 !was 0.010D0
        JPH   = 0.143D0
        OMEGP = 0.0D0

        KAPPA = 0.0D0    ! 1.0D0 !|0-1: Cubic to ~contact homotopy parameter
        BETA = 10.D0  ! 1.5 3 5  10  |btw 1.25 to 125 see explanation in func_ode45_tanh.m or betlek.31.01.2021
        K = 150.D0    ! 30 100  |see explanation in func_ode45_tanh.m or betlek.31.01.2021

        RHO   = 0.0D0 ! 1.0D0 !|0-1: Cubic to Power5(quintic) homotopy parameter
        KSI   = 0.1D0 ! 0.1D0 !|0-0.2: Quintic stiffness ratio, (k_5*c^4)/k_r 

        PAR(1)=GAMMA 
        PAR(2)=OMEG 
        PAR(3)=MH 
        PAR(4)=JPH
        PAR(5)=EPSH 
        PAR(6)=ZETA 
        PAR(7)=OMEGP 
        PAR(8)=KAPPA
        PAR(10)=BETA
        PAR(13)=K  
        PAR(14)=KSI
        PAR(15)=RHO

        TPI=8*ATAN(1.0D0) ! = 2*pi 
        PAR(11)= 2.8172d0 !|Tkn frm autocorper.m function     !TPI/OMEG  ! Period = 2*pi/omeg
        !|:What s d period f d system at internal resonance <<< 2pi/GCF(FW,BW) Greatest common factor ; 
        !|::For near-zero amp orbits (like 1e-5 amp) FW and BW correspond to Campbell plot values. 
        !|::Bt higher amp moves d freqs away frm Campbell, So Per s tb calcD w autocorper.m func fr d non0epsH dat orbits (amps like 0.5)

        X = SIN(TPI*T) 
        Y = COS(TPI*T)

      END SUBROUTINE STPNT


      ! SUBROUTINE PVLS(NDIM,U,PAR)
      ! IMPLICIT NONE
      ! INTEGER, INTENT(IN) :: NDIM
      ! DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      ! DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      ! DOUBLE PRECISION R2,Q1,Q2
      !   Q1 = U(1)
      !   Q2 = U(2)

      !   ! R2 = max(U(1)**2 + U(2)**2)  #|!WORKj
      !   R2 = Q1**2 + Q2**2 
      !   ! Which one is this?
      !   ! ?? max( U(1) )^2 + max( U(2) )^2
      !   ! ?? max( U(1)^2 + U(2)^2 )

      !   PAR(9) = R2**(0.5)
      ! END SUBROUTINE PVLS

      SUBROUTINE PVLS(NDIM,U,PAR)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
        DOUBLE PRECISION, ALLOCATABLE :: UU(:,:),R2(:),R(:)
        DOUBLE PRECISION, EXTERNAL :: GETP, GETU
        DOUBLE PRECISION :: DUMM
        INTEGER :: NDX,NCOL,NTST,i,j,k
        !|SEE the explanation of demo pvl.f90 : 
        !|:"For algebraic problems the argument U is, as usual, the state vector.
        !|:For differential equations the argument U represents the approximate 
        !|:solution on the entire interval [0,1]. In this case its values can
        !|:be accessed indirectly by calls to GETP, as illustrated below, or
        !|:by obtaining NDIM, NCOL, NTST via GETP and then dimensioning U as
        !|:U(NDIM,0:NCOL*NTST) in a seperate subroutine that is called by PVLS. "

        NDX  = NINT(GETP('NDX',0,U))
        NTST = NINT(GETP('NTST',0,U))
        NCOL = NINT(GETP('NCOL',0,U))

        !|>>|METHOD-A/B|>> METHOD A/B  
        !|:CORRECT AMPLITUDE FOR IPS=2
        DUMM = GETU(i,j,U,NDX,NTST,NCOL,PAR(2)) 
        !|:Will get amp and print the states to dummy_write_file.txt
        PAR(9) = DUMM 
        !|__|METHOD-A/B|__
        
        !|..|METHOD-C|.. METHOD C
        ! DO i=1,NDX
        !   DO j=0,NCOL*NTST
        !     UU(i,j) = GETU(i,j,U,NDX,NTST,NCOL)
        !   END DO
        ! END DO
        ! R2(:) = UU(1,:)**2 + UU(2,:)**2  !El by El multip
        ! R(:) = R2**0.5 
        ! PAR(9) = MAXVAL(R)
        !|__|METHOD-C|__

        !|>>|METHOD-D|>> DEFAULT METHOD - PHASED AMP
        ! PAR(9) = (U(1)**2 + U(2)**2)**0.5 
        !|:The initial data of one period is sampled for amplitude.
        !|:See the qoute above.
        !|__|METHOD-D|__

      END SUBROUTINE PVLS


      DOUBLE PRECISION FUNCTION GETU(i,j,U,NDX,NTST,NCOL,MYPAR)
        INTEGER, INTENT(IN) :: NDX,NCOL,NTST,i,j
        INTEGER :: p, I_MAX
        DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
        DOUBLE PRECISION, ALLOCATABLE :: R2(:),R(:)
        DOUBLE PRECISION :: M,N,MYPAR
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

        GETU = M

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
