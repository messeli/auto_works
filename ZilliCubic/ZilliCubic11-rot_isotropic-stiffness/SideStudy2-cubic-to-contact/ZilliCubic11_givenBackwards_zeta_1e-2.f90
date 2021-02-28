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
        DOUBLE PRECISION ZETA,A,BETA,Q1,Q2,Q3,Q4,C,R2,R,GAMMA,KAPPA,K,OMEG,OMEGP,MH,EPSH,JPH,FSNUB_v,FSNUB_u,FCUBIC_u,FCUBIC_v
      
        GAMMA = PAR(1) 
        OMEG = PAR(2) !|Rotor speed
        MH = PAR(3) 
        JPH = PAR(4) 
        EPSH = PAR(5) !|Eccentricity
        ZETA = PAR(6) 
        OMEGP = PAR(7) !|Rotor acceletation
        KAPPA = PAR(8) !|Nonlinearity Homotopy. 0:cubic, 1:contact
        BETA = PAR(10) !|Snubber ring stiffness measure: ks/kr 
        K = PAR(13) 

        Q1=U(1)
        Q2=U(2)
        Q3=U(3)
        Q4=U(4)
 
        R2 = (Q1**2)+(Q2**2)
        R = R2**0.5D0
        C = 1 
        
        FSNUB_u = -0.5D0*(tanh(K*(R-1))+1)*(1-1/R)*BETA*Q1  
        FSNUB_v = -0.5D0*(tanh(K*(R-1))+1)*(1-1/R)*BETA*Q2  
        FCUBIC_u = -GAMMA*R2*Q1
        FCUBIC_v = -GAMMA*R2*Q2
  
        F(1) = Q3
        F(2) = Q4
        F(3) = -OMEG*(JPH-2)*Q4          &
               -2*ZETA*Q3                &
               +((OMEG**2)*(1-JPH)-1)*Q1 &
               +2*ZETA*OMEG*Q2           &
               +MH*EPSH*(OMEG**2)        &
               +KAPPA*FSNUB_u+(1-KAPPA)*FCUBIC_u !|NONLINEARITY HOMOTOPY
        F(4) = +OMEG*(JPH-2)*Q3          &
               -2*ZETA*Q4                &
               +((OMEG**2)*(1-JPH)-1)*Q2 &
               -2*ZETA*OMEG*Q1           &
               -MH*EPSH*OMEGP            &
               +KAPPA*FSNUB_v+(1-KAPPA)*FCUBIC_v !|NONLINEARITY HOMOTOPY        

      ! IF (IJAC.EQ.1) RETURN
      !   DFDU(1,1)=0
      ! ...
   
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: T
        DOUBLE PRECISION ZETA,BETA,Q1,Q2,Q3,Q4,S4,C4,R2,GAMMA,KAPPA,K,OMEG,OMEGP,MH,EPSH,JPH 

        ! GIVEN BACKWARDS (MATLAB LOWER AMPLITUDE END DIMENSIONS)
        GAMMA = 0.25D0 !0.25 !2nd continue gamma from 0 to 0.25
        OMEG = 5.0D0 ! Continue OMEG from 0 to 7. 
        MH = 0.9D0 !1st continue MH from 0 to 0.9
        EPSH =0.353D0
        ZETA = 0.01D0  !TRUE ONE IS 1e-2;BUT, TRY 0.05 for avoiding backward continuation from simulation datum
        JPH  = 0.143D0
        OMEGP = 0.D0
  
        KAPPA = 1.D0  ! 1.D0
        BETA = 10.D0  ! 1.5 3 5 10  |btw 1.25 to 125 see explanation in func_ode45_tanh.m or betlek.31.01.2021
        K = 150.D0    ! 30 100  |see explanation in func_ode45_tanh.m or betlek.31.01.2021

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

        ! !| Orig one used for orig givenBackwards, where it was zeta=0.01 
        ! ! U(1)=-0.380009000620448
        ! ! U(2)= 0.000022150720618
        ! ! U(3)=-0.001297860009550
        ! ! U(4)= 0.000206860716715

        ! !| New ones from "Google Drive\PhD\Zilli ISO cubic stiffness\Resutls ISO cubic\
        ! !| ...7p0 end lowAmp datum, differentZeta for AUTO givenBackwards.txt"
        ! U(1)=-0.380085547357884
        ! U(2)=0.000000023340008
        ! U(3)=-0.001299216246013
        ! U(4)=0.000000244301170
        
        !| Take the stpnt below from matlab datum at Omeg=7.01 (or whatever it is in MATLAB!)
        ! U(1) = -3.79737E-01
        ! u(2) = -2.59377E-03
        ! U(3) = 0.00000E+00
        ! U(4) = 0.00000E+00
        !| Take the stpnt below from matlab datum at Omeg=5.00 beta=1.5 K=30
        U(1) = -3.88852349524439D-01
        u(2) = -1.903803765019D-03
        U(3) = 1.2420423D-8
        U(4) = -1.06246196D-7

      END SUBROUTINE STPNT


      ! SUBROUTINE PVLS(NDIM,U,PAR)
      ! IMPLICIT NONE
      ! INTEGER, INTENT(IN) :: NDIM
      ! DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      ! DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      ! DOUBLE PRECISION R2,Q1,Q3
      !   Q1 = U(1)
      !   Q3 = U(3)

      !   R2 = Q1**2 + Q3**2 
      !   ! Which one is this?
      !   ! ?? max( U(1) )^2 + max( U(3) )^2
      !   ! ?? max( U(1)^2 + U(3)^2 )

      !   PAR(9) = R2**(0.5)
      ! END SUBROUTINE PVLS

      SUBROUTINE PVLS(NDIM,U,PAR)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
        DOUBLE PRECISION, ALLOCATABLE :: UU(:,:),R2(:),R(:)
        DOUBLE PRECISION, EXTERNAL :: GETP, GETU
        INTEGER :: NDX,NCOL,NTST,i,j,k
        !|SEE the explanation of demo pvl.f90 : 
        !|:"For algebraic problems the argument U is, as usual, the state vector.
        !|:For differential equations the argument U represents the approximate 
        !|:solution on the entire interval [0,1]. In this case its values can
        !|:be accessed indirectly by calls to GETP, as illustrated below, or
        !|:by obtaining NDIM, NCOL, NTST via GETP and then dimensioning U as
        !|:U(NDIM,0:NCOL*NTST) in a seperate subroutine that is called by PVLS. "

        ! NDX  = NINT(GETP('NDX',0,U))
        ! NTST = NINT(GETP('NTST',0,U))
        ! NCOL = NINT(GETP('NCOL',0,U))

        !|AAA/BBB METHOD A/B
        ! PAR(9) = GETU(i,j,U,NDX,NTST,NCOL)
        !|AAA/BBB END
        
        !|CCC METHOD C
        ! DO i=1,NDX
        !   DO j=0,NCOL*NTST
        !     UU(i,j) = GETU(i,j,U,NDX,NTST,NCOL)
        !   END DO
        ! END DO
        ! R2(:) = UU(1,:)**2 + UU(3,:)**2  !El by El multip
        ! R(:) = R2**0.5 
        ! PAR(9) = MAXVAL(R)
        !|CCC END

        !|DDD DEFAULT METHOD - PHASED AMP
        !|:CORRECT AMPLITUDE FOR IPS=1
        PAR(9) = (U(1)**2 + U(2)**2)**0.5 
        !|:The state vector is directly used to calc amplitude.
        !|:See the qoute above.
        !|DDD END 
      END SUBROUTINE PVLS


      DOUBLE PRECISION FUNCTION GETU(i,j,U,NDX,NTST,NCOL)
        INTEGER, INTENT(IN) :: NDX,NCOL,NTST,i,j
        INTEGER :: p
        DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
        DOUBLE PRECISION, ALLOCATABLE :: R2(:),R(:)
        DOUBLE PRECISION :: M,N

        !|AAA METHOD A
        M = 0.0
        DO p=0,NCOL*NTST
          N = ( U(1,p)**2+U(3,p)**2 )**0.5
          if (N.GT.M) THEN
            M = N
          END IF 
        END DO 

        GETU = M
        !|AAA END

        !|BBB METHOD B
        ! R2(:) = U(1,:)**2 + U(3,:)**2  !
        ! R(:) = R2**0.5 
        ! GETU = MAXVAL(R)
        !|BBB END

        !|CCC METHOD C
        ! GETU = U(i,j)
        !|CCC END
      END FUNCTION GETU



      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

