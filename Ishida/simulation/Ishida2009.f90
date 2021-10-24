! 2009Ishida et al works reproduce with tanh-contact 

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
        DOUBLE PRECISION OMEG,DELT,BETA,IP,E,F0,C1,C2,ALPH,REST,M_2,KB_N,CB_N,K
        DOUBLE PRECISION Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,X,Y
        DOUBLE PRECISION R_2, CA, F_n
      
        OMEG = PAR(1) !|rotor speed
        DELT = PAR(2) !|clearance
        BETA = PAR(3) !|bearing mass (m2) over disk mass (m1)
        IP = PAR(4)   !|polar moment of inertia
        E = PAR(5)    !|eccentricity of the disk mass centre.
        F0 = PAR(6)   !|bearing misalignment
        C1 = PAR(7)   !|disk damping coefficient
        C2 = PAR(8)   !|lower bearing damping coefficient
        ALPH = PAR(9) !|a/l position of the disk
        REST = PAR(10)!|coefficient of restitution
        !PAR(11) reserved for period
        !PAR(12) reserved for torus angle if it is encountered sw. 
        M_2 = PAR(13) !|bearing mass 
        !PAR(14) could be the independent time, as I gather from ivp demo. 
        KB_N = PAR(15)!|contact stiffness
        CB_N = PAR(16)!|contact damping (calcD from c.restitution)
        K = PAR(17)   !|tanh steppness


        Q1=U(1) !x1
        Q2=U(2) !y1
        Q3=U(3) !x2
        Q4=U(4) !y2
        Q5=U(5) !x1_dot
        Q6=U(6) !y1_dot
        Q7=U(7) !x2_dot
        Q8=U(8) !y2_dot
        X =U(9) !sin(OMEG*T)
        Y =U(10)!cos(OMEG*T)

        R_2 = sqrt(Q3**2 + Q4**2) !bearing position
        CA = atan2(Q4,Q3) !contact angle 
        
        F_n  = +0.5D0*(tanh(K*(R_2-DELT))+1) * (   KB_N*(R_2-DELT) + CB_N*(Q7*cos(CA)+Q8*sin(CA))   )

        F(1) = + Q5
        F(2) = + Q6
        F(3) = + Q7
        F(4) = + Q8
        F(5) = - IP*OMEG*Q6 - C1*Q5 - (Q1 - ALPH*Q3) + E*OMEG**2* Y
        F(6) = + IP*OMEG*Q5 - C1*Q6 - (Q2 - ALPH*Q4) + E*OMEG**2* X + F0 !X=sin(OMEG*T), Y=cos(OMEG*T)
        F(7) = - C2/BETA*Q7 - ALPH/BETA*(ALPH*Q3 - Q1) - Q3/R_2*F_n/BETA
        F(8) = - C2/BETA*Q8 - ALPH/BETA*(ALPH*Q4 - Q2) - Q4/R_2*F_n/BETA
        F(9) = + X + OMEG*Y - X*(X**2 + Y**2) !|To remember the nonlinear oscillator   | x = sin(OMEG*T) , y = cos(OMEG*T) 
        F(10)= + Y - OMEG*X - Y*(X**2 + Y**2) !|...remember 1st one, 2nd swaps x and y, and makes OMEG negative. 


      ! IF (IJAC.EQ.1) RETURN
      !   DFDU(1,1)=0
      ! ... 
      END SUBROUTINE FUNC
                
      SUBROUTINE STPNT(NDIM,U,PAR,T)  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*) !|U's assumed dimension is rank-1, so only one col of data it is: the starting point.
        DOUBLE PRECISION, INTENT(IN) :: T
        DOUBLE PRECISION OMEG,DELT,BETA,IP,E,F0,C1,C2,ALPH,REST,M_2,KB_N,CB_N,K
        DOUBLE PRECISION Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,X,Y
        DOUBLE PRECISION PI


        OMEG = 2.19D0
        DELT = 0.1D0 
        BETA = 0.019D0 
        IP   = 0.1D0   
        E    = 0.5D0    
        F0   = 0.5D0
        C1   = 0.005D0  
        C2   = 0.05D0   
        ALPH = 0.611D0 !|a/l position of the disk
        REST = 0.98D0 !|2009Ishida does not say what value they used. This is my pick (takes abit less time to simulate than 0.8)
        M_2  = 0.3D0 !My pick, not given in the paper. 
        KB_N = 1.0D5  !look up formulas from s90 of b 1985Johnson Contact Mechanics Herzian contact ??
        CB_N = sqrt(   (2 * log10(REST)**2 * M_2*KB_N) / (pi**2 + log10(REST)**2)   ) 
        !PAR(11) reserved for period
        !PAR(12) reserved for torus sth ? 
        K = 150.0D0   !|tanh steppness

        PAR(1) = OMEG
        PAR(2) = DELT
        PAR(3) = BETA
        PAR(4) = IP
        PAR(5) = E
        PAR(6) = F0
        PAR(7) = C1
        PAR(8) = C2
        PAR(9) = ALPH !|a/l position of the disk
        PAR(10) = REST 
        !PAR(11) reserved for period, see doc/sect11.1
        !PAR(12) reserved for torus angle, see doc/sect11.1 
        PAR(13) = M_2
        !PAR(14) stores the independent time, see ivp demo and doc/sect11.1
        PAR(15) = KB_N
        PAR(16) = CB_N
        PAR(17) = K   !|tanh steppness

        !| Initial conditions of the IMPLICIT EULER TIME INTEGRATION << The same as that in MATLAB. 
        U = (/0.01d0, -0.012d0, -0.020d0, 0.018d0  ,  0.013d0, -0.010d0, -0.008d0, 0.024d0  ,  0.d0, 1.d0/) 



        ! PI = 4.D0*ATAN(1.0D0) !! = pi 
        ! PAR(11) = 2.D0*PI / (OMEG/2.D0)   !! PAR(11) is always reserved for period in AUTO. 
        !PAR(11) = 5.73813D0 !|This is autom taken from the .dat file, but can spec here. 
        !|:Period = 2*pi/GCF(OMEG,OMEG/2.D0)
        !|::In the synchronous whirl the period of oscillation is the same as the rotor speed!!
        !|::In the rotating frame synch freq is 0; so here we have: w_sta=w_rot+Omeg = Omeg. 
        !|:In the sta-frame subharmonic response of order 1/2, the freq content has OMEG and OMEG/2. 
        !|::So the EBOB of the freqs is OMEG/2. Therefore, the period is 2.D0*PI/(OMEG/2.D0) !BN
        !|:Because we define the PAR(11) here explicitly, the time variable becomes 
        !|...scaled to time/period=T in range 0-1.

        ! X = SIN(2.D0*PI*T) 
        ! Y = COS(2.D0*PI*T)
        !|:No need to start any state var as all states' data will come from the .dat file already. 
        
        ! THETA = OMEG*(TPI*T)
        ! THETA = ATAN2(X,Y) !|Position of eccentric mass (origY at x-axis) 
        !|:ATAN2(Y,X) angle of X+jY in correct quadrant
        !|:We have X=sin ,Y=cos; so  ATAN2(X,Y) should be correct !BN  

        !The starting data will be the .dat file from MATLAB, so it should include the X,Y column.
        !...This is done by populating cols as [t,x1,y1,x2,y2,dx1,dy1,dx2,dy2,X,Y]from T=0 to 1 
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

        NDX  = NINT(GETP('NDX',0,U)) !|GETP func is built-in in AUTO. GETU is handwritten.
        NTST = NINT(GETP('NTST',0,U))
        NCOL = NINT(GETP('NCOL',0,U))

        !|___METHOD-A/B> METHOD A/B  
        !|:CORRECT AMPLITUDE FOR IPS=2
        DUMM = GETU(i,j,U,NDX,NTST,NCOL,PAR(2)) 
        !|:Will get amp and print the states to dummy_write_file.txt
        PAR(18) = DUMM 
        !|___METHOD-A/B.
        
        !|___METHOD-C> METHOD C
        ! DO i=1,NDX
        !   DO j=0,NCOL*NTST
        !     UU(i,j) = GETU(i,j,U,NDX,NTST,NCOL)
        !   END DO
        ! END DO
        ! R2(:) = UU(1,:)**2 + UU(2,:)**2  !El by El multip
        ! R(:) = R2**0.5 
        ! PAR(9) = MAXVAL(R)
        !|___METHOD-C.

        !|___METHOD-D> DEFAULT METHOD - PHASED AMP
        ! PAR(9) = (U(1)**2 + U(2)**2)**0.5 
        !|:The initial data of one period is sampled for amplitude.
        !|:See the qoute above.
        !|___METHOD-D. 

        PAR(19) = GETP('MAX',7,U) !xdot
        PAR(20) = GETP('MAX',8,U) !ydot
        PAR(21) = GETP('MAX',9,U) !X (ie sin(Omeg*t))
        PAR(22) = GETP('MAX',10,U)!Y (ie cos(Omeg*t))

      END SUBROUTINE PVLS


      DOUBLE PRECISION FUNCTION GETU(i,j,U,NDX,NTST,NCOL,MYPAR)
        INTEGER, INTENT(IN) :: NDX,NCOL,NTST,i,j
        INTEGER :: p, I_MAX
        DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST) 
        !|See at Keep "Fortran assumed dimensions and how it is interpreted; Some printing"
        DOUBLE PRECISION, ALLOCATABLE :: R2(:),R(:)
        DOUBLE PRECISION :: M,N , MYPAR
        !| MYPAR - the parameter to print in the 1st col of dummy_write_file.txt

        !|___METHOD-A>
        M = 0.D0
        DO p=0,NCOL*NTST
          N = ( U(1,p)**2+U(2,p)**2 )**0.5D0 !|Amplitude of the disk, not bearing!
          IF (N.GT.M) THEN
            M = N
            I_MAX = p
          END IF 
        END DO 

        !| Print the amplitude point's all 4 states.
        !| ...E.g for time simulation initial data, 
        !| ...parse the following file.  
        ! OPEN (unit = 10, file = "dummy_write_file.txt",access = "append")
        ! WRITE(10,*) MYPAR,U(1,I_MAX),U(2,I_MAX),U(3,I_MAX),U(4,I_MAX),M !|"10" means write into the file named fort.10 .
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

