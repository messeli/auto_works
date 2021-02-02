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
        DOUBLE PRECISION ZETA,BETA,Q1,Q3,Q2,Q4,R2,R,C,GAMMA,KAPPA,K,OMEG,OMEGP,MH,EPSH,JPH,FSNUB_u,FSNUB_v,FCUBIC_v,FCUBIC_u
      
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

        Q1=U(1)
        Q2=U(2)
        Q3=U(3)
        Q4=U(4)

        R2 = (Q1**2)+(Q2**2)
        R = R2**0.5D0
        C = 1
        
        ! WRITE (*,*) R
        ! IF (R.LT.0.1D0) THEN
        !   FSNUB_u = 0
        !   FSNUB_v = 0
        ! ELSE 
        !   FSNUB_u = -0.5D0*(tanh(K*(R-C))+1) * ABS(1-1.D0/R)*BETA*Q1  
        !   FSNUB_v = -0.5D0*(tanh(K*(R-c))+1) * ABS(1-1.D0/R)*BETA*Q2  
        ! ENDIF



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
               +KAPPA*FSNUB_u+(1-KAPPA)*FCUBIC_u !|NONLINEARITY HOMOTOPYs
        F(4) = +OMEG*(JPH-2)*Q3          &
               -2*ZETA*Q4                &
               +((OMEG**2)*(1-JPH)-1)*Q2 &
               -2*ZETA*OMEG*Q1           &
               -MH*EPSH*OMEGP            &
               +KAPPA*FSNUB_v+(1-KAPPA)*FCUBIC_v !|NONLINEARITY HOMOTOPY

      ! IF (IJAC.EQ.1) RETURN
      !   DFDU(1,1)=0
      !   ...
      
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: T
        DOUBLE PRECISION ZETA,Q1,Q3,Q2,Q4,S4,C4,R2,R,GAMMA,KAPPA,K,OMEG,OMEGP,MH,EPSH,JPH,BETA,TPI!WN,F_M,K3_M,OMEG,

        ! FORWARDS (MATLAB SINE SWEEP GIVES THIS TOO)
        GAMMA = 0.0 !2nd continue gamma from 0 to 0.25
        OMEG = 0.1 ! Continue OMEG from 0 to 7. 
        MH = 0.0 !1st continue MH from 0 to 0.9

        EPSH =0.353
        ZETA = 0.02 !0.01
        JPH  = 0.143
        OMEGP = 0.0
        KAPPA = 0.0  ! 1.D0
        BETA = 0.D0  ! 10.D0
        K = 0.D0     ! 10.D0

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

        U(1)=0.01
        U(2)=0.01
        U(3)=0
        U(4)=0         
      END SUBROUTINE STPNT


      SUBROUTINE PVLS(NDIM,U,PAR)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
        DOUBLE PRECISION, ALLOCATABLE :: UU(:,:),R2(:),R(:)
        DOUBLE PRECISION, EXTERNAL :: GETP, GETU
        INTEGER :: NDX,NCOL,NTST,i,j,k
        !|:SEE the explanation of demo pvl.f90 : 
        !|::"For algebraic problems the argument U is, as usual, the state vector.
        !|::For differential equations the argument U represents the approximate 
        !|::solution on the entire interval [0,1]. In this case its values can
        !|::be accessed indirectly by calls to GETP, as illustrated below, or
        !|::by obtaining NDIM, NCOL, NTST via GETP and then dimensioning U as
        !|::U(NDIM,0:NCOL*NTST) in a seperate subroutine that is called by PVLS. "

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
        ! R2(:) = UU(1,:)**2 + UU(2,:)**2  !El by El multip
        ! R(:) = R2**0.5 
        ! PAR(9) = MAXVAL(R)
        !|CCC END

        !|DDD DEFAULT METHOD - PHASED AMP
        !|:CORRECT AMPLITUDE FOR IPS=1, which this .f90 is used for.
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
          N = ( U(1,p)**2+U(2,p)**2 )**0.5
          if (N.GT.M) THEN
            M = N
          END IF 
        END DO 

        GETU = M
        !|AAA END

        !|BBB METHOD B
        ! R2(:) = U(1,:)**2 + U(2,:)**2  !
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

