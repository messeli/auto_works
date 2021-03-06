!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!                   duf1 :      free duffing oscillator
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
        DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
        DOUBLE PRECISION k

        k=PAR(1)
       
	      F(1)= U(2)
        F(2)=-k*U(1)-U(1)**3
      END SUBROUTINE FUNC


      SUBROUTINE STPNT(NDIM,U,PAR,T)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: T

        PAR(1)=-1.
       
        U(1)=0.
        U(2)=0.
      END SUBROUTINE STPNT


      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
