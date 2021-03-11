#| Runnung this script tool 100 seconds. 
#| Plotting took 
import numpy as n
import matplotlib.pyplot as p
from mpl_toolkits.mplot3d import Axes3D 
import func_ode45
import time 
#|EEE> DEFINE SOLVER INPUTS
start_sim = time.time()

#| Inputs
qn_ = n.array([0,0.6,-1.2460,0]).reshape(4,1)
tend = 1000 # 3000
tol = 1e-6 # 1e-7
Omeg_range = n.arange(5.1,7.2,0.1)  # ".01:0.01:7.01" 
print("Omeg_range length is {}".format(len(Omeg_range)) )

iRange = range(len(Omeg_range))  

T=[]  
Q=[] 
for i in iRange: #|parfor
  #|___INIT___|> INITIAL CONDITION
  Omeg = Omeg_range[i] 
  qn = qn_ 
  
  T_X2rot = n.array([[1,0,0,0],[0,1,0,0],[0,-Omeg,1,0],[Omeg,0,0,1]]) #|q_X = T_X2rot . q_rot
  T_phi2X = n.array([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]]) #|q_phi = T_phi2X . q_X
  T_X2rot_INV = n.linalg.inv(T_X2rot)  
  T_phi2X_INV = n.linalg.inv(T_phi2X) 
  # Transform_INV = n.dot(T_X2rot_INV,T_phi2X_INV) #| or n.matmul(A,B)
  # qn = n.dot(Transform_INV,qn) #|q_u = T_X2rot^(-1) . T_phi2X^(-1) . q_phi
  qn = T_X2rot_INV @ T_phi2X_INV @ qn 
  #|___INIT___|.

  #|___SOL___|> Given solver is the ode45 solver
  (T_new,Q_new) = func_ode45.func_ode45_v1(Omeg,qn,tend,tol) 
  T.append(T_new) 
  Q.append(Q_new) 
  #|___SOL___|. 

#| Stats
end_sim = time.time() 
time_sim = end_sim - start_sim 
print("time_sim is {}".format(time_sim) )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
start_plot = time.time()
#|888> INDIVIDUAL PLOT  #| STILL MATLAB , NOT IN PYTHON SCRIPT
#   for Omeg_check = 2.91 # Omeg_range(1)
#     I = find( abs( Omeg_range-Omeg_check ) <= 0.001 );
#     Zilli_individualplot_ISOcubicStiffnessNoContact(T{I},Q{I},...
#                                             Omeg_check,[0.80 1.0],false); 
#     #|:TF - Animation
#   end
#|888.

#|NNN> 3D&2D PLOTS, POINCARE SECTIONS
#|OOO> PREPARE FIGURE
colour = (1.0,0.7,0.0) 
fig = p.figure(figsize=(10.0,3.0)) #figure
ax1 = fig.add_subplot(121,projection='3d')  #|D ",projection='3d'" part adds Axes3D type axes
ax2 = fig.add_subplot(1,2,2)  # ax2 = subplot(122) ; hold(ax2,'on') ; grid(ax2,'on') ; axis(ax2, 'square')

ax1.set_xlabel("u")  # xlabel(ax1,"$\it\hat{u}$","interpreter","latex","FontSize",14) 
ax1.set_ylabel("v")  # ylabel(ax1,"$\it\hat{v}$","interpreter","latex","FontSize",14)
ax1.set_zlabel("Omeg")  # zlabel(ax1,"Rotor speed","fontSize",10) #|Omega^H^a^t is Omeg
# title (ax1,[solver+" "+frame,"Omeg-range "+Omeg_range_txt," "],"FontSize",10)
ax1.set_xlim(-3,3)  # xlim(ax1,[-3 3]); 
ax1.set_ylim(-3,3)  # ylim(ax1,[-3 3])

ax2.set_xlabel("Omeg"     )  # xlabel(ax2,"$$\rm Rotor\,Speed,\,\it\hat{\Omega}$$",'Interpreter','latex',"FontSize",14) #|Omega^H^a^t=Omeg
ax2.set_ylabel("Amplitude")  # ylabel(ax2,"$\rm Amplitude,\,\it\hat{r}$","interpreter","latex","FontSize",14) 
# title(ax2,["phi-EOM initial conditions",num2str(qn_')," "],"fontSize",10)
ax1.set_ylim(-3,3)  # ylim(ax2,[0 5]) 
#|OOO. 

for i in range(0,len(Omeg_range)):
  #| FIND INDEXES TO PLOT FROM PERCENTAGES
  from_ = int(round( 0.99 * Q[i].shape[1] , 0 )) 
  to = Q[i].shape[1] 
  dataLength = len(range(from_,to))  
  Omeg = Omeg_range[i] 

  #|111> 3D PLOT - STACK
  ax1.plot( Q[i][0,from_:to], Q[i][1,from_:to], Omeg, 'k' )
  #|111.

  #|222> 2D PLOT - AMPLITUDE
  rH = ( Q[i][0,:]**2 + Q[i][1,:]**2 )**0.5 
  ax2.plot( Omeg*n.ones(dataLength) , rH[from_:to] , "k.", Omeg , max(rH[from_:to]), "bo")  #|Nt:(1, 104) and (104,) does not match but (104,1) and (104,) does match
  #|222.

  #|PPP> STROBOSCOPIC POINCARE SECTION 
  #|P_III> GRAB PERIOD-CLOSE INDEXES 
  TPeriod = 2*n.pi/Omeg 
  k = 1  
  r = 1  
  IPeriod = []  
  fromP = int(round( 0.60 * Q[i].shape[1] , 0 )) 
  toP = -1 # Q[i].shape[1] 
  while k*TPeriod <= T[i][-1] :
    if k*TPeriod > T[i][fromP] and k*TPeriod <= T[i][toP]:
      IPeriod.append( n.argmin( abs( T[i] - k*TPeriod ) ) )
      r = r+1 #|IPeriod index
    k = k+1 #|Period's number.
  #|P_III.
  
  ax1.plot( Q[i][0,IPeriod], Q[i][1,IPeriod], Omeg,
            ".",linewidth=1,color=colour)
  ax2.plot( Omeg*n.ones(len(IPeriod)) , rH[IPeriod], ".",linewidth=1,color=colour)

  #|PPP. POINCARE
  #|Continue to the next Omeg in Omeg_range
p.savefig('myplot.png')
p.show()
#|NNN. 3D&2D 
#| Stats
stop_plot = time.time()
time_plot = stop_time - start_time 
print("time_plot is {}".format(time_plot) )