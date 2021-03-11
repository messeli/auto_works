#| Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
import numpy as n
from scipy.integrate import solve_ivp

#|___V1___|> NESTED FUNCS (more like in d MATLAB)
def func_ode45_v1(Omeg,qn,tend,tol):
   def func_rot(tn, qn):
      JpH = 0.143  ; epsH = 3.53e-1 ; zeta = 0.010 ;
      gamma = 0.25 ; mH    = 0.9    ;

      q1 = qn[0] ; q2 = qn[1] ; q3 = qn[2] ; q4 = qn[3] ;
      
      r2 = q1**2 + q2**2 ;

      dq1 = q3 ; 
      dq2 = q4 ;    
      dq3 = - Omeg * (JpH-2) * q4 \
               - 2 * zeta * q3 \
               + ( Omeg**2 * (1-JpH)  -  1 ) * q1 \
               + 2 * zeta * Omeg * q2 \
               + mH*epsH*Omeg**2 \
               - gamma * r2 * q1 ; 
      dq4 = + Omeg * (JpH-2) * q3 \
               - 2 * zeta * q4 \
               + ( Omeg**2 * (1-JpH)  -  1 ) * q2 \
               - 2 * zeta * Omeg * q1 \
               - gamma * r2 * q2 ;

      dq = n.array([dq1,dq2,dq3,dq4]).reshape(4) ;
      return dq

   def contactevent():
      pass
      pass

   print(Omeg)
   qn = n.array(qn).reshape(4)

   sol1 = solve_ivp(func_rot, [0,tend], qn, rtol=tol, atol=tol*1e-2) ;
   Q_new = sol1.y ;#|numpy.ndarray, Q_new.shape=(4, 23850), rows (23850,)
   T_new = sol1.t ;#|numpy.ndarray, T_new.shape=(23850,)
   return (T_new,Q_new)
#|___V1___|.



#|___V2___|> SEPARATE FUNCS
def func_ode45_v2(Omeg,qn,tend,tol):
   global Omeg_shared
   Omeg_shared = Omeg #|Here "Omeg_shared = float(Omeg)" not needD 
   print(Omeg_shared)
   qn = n.array(qn).reshape(4)


   sol1 = solve_ivp( func_rot, [0,tend], qn, rtol=tol, atol=tol*1e-2 , method="RK45", events=None) ;

   Q_new = sol1.y ;#|numpy.ndarray, Q.shape = (4, 23850)
   T_new = sol1.t ;#|numpy.ndarray, T.shape = (23850,)
   return (T_new,Q_new)

def func_rot(tn, qn):
   JpH = 0.143  ; epsH = 3.53e-1 ; zeta = 0.010 ;
   gamma = 0.25 ; mH    = 0.9    ;

   q1 = qn[0] ; q2 = qn[1] ; q3 = qn[2] ; q4 = qn[3] ;
   
   r2 = q1**2 + q2**2 ;

   dq1 = q3 ; 
   dq2 = q4 ;    
   dq3 = - Omeg_shared * (JpH-2) * q4 \
            - 2 * zeta * q3 \
            + ( Omeg_shared**2 * (1-JpH)  -  1 ) * q1 \
            + 2 * zeta * Omeg_shared * q2 \
            + mH*epsH*Omeg_shared**2 \
            - gamma * r2 * q1 ; 
   dq4 = + Omeg_shared * (JpH-2) * q3 \
            - 2 * zeta * q4 \
            + ( Omeg_shared**2 * (1-JpH)  -  1 ) * q2 \
            - 2 * zeta * Omeg_shared * q1 \
            - gamma * r2 * q2 ;

   dq = n.array([dq1,dq2,dq3,dq4]).reshape(4) ;
   return dq

def contactevent():
   pass
   pass
#|___V2___|.