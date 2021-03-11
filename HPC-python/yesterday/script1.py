#!/usr/bin/env python3
import numpy as n
import matplotlib.pyplot as p
from mpl_toolkits.mplot3d import Axes3D 
import func_ode45
import time 
#|EEE> DEFINE SOLVER INPUTS
start_sim = time.time()
qn_ = n.array([0,0.6,-1.2460,0]).reshape(4,1)

tend = 1000 
tol = 1e-6 
Omeg_range = n.arange(3.1,7.2,0.1)  # ".01:0.01:7.01" 
print("Omeg_range length is {}".format(len(Omeg_range)) )

iRange = range(len(Omeg_range))  

T=[]  
Q=[] 
for i in iRange:  
  Omeg = Omeg_range[i] 
  qn = qn_ 
  
  T_X2rot = n.array([[1,0,0,0],[0,1,0,0],[0,-Omeg,1,0],[Omeg,0,0,1]]) #|q_X = T_X2rot . q_rot
  T_phi2X = n.array([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]]) #|q_phi = T_phi2X . q_X
  T_X2rot_INV = n.linalg.inv(T_X2rot)  
  T_phi2X_INV = n.linalg.inv(T_phi2X) 
  qn = T_X2rot_INV @ T_phi2X_INV @ qn 

  
  (T_new,Q_new) = func_ode45.func_ode45_v1(Omeg,qn,tend,tol) 
  T.append(T_new) 
  Q.append(Q_new) 
  
end_sim = time.time() 
time_sim = start_sim - end_sim 
print("time_sim is {}".format(time_sim) )



#|FFF> PLOTTING CURRENT COORD RESULTS
start_plot = time.time()

colour = (1.0,0.7,0.0) 
fig = p.figure(figsize=(10.0,3.0)) 
ax1 = fig.add_subplot(121,projection='3d')  
ax2 = fig.add_subplot(1,2,2)  

ax1.set_xlabel("u") 
ax1.set_ylabel("v") 
ax1.set_zlabel("Omeg") 

ax1.set_xlim(-3,3)  
ax1.set_ylim(-3,3) 

ax2.set_xlabel("Omeg"     )  
ax2.set_ylabel("Amplitude")  

ax1.set_ylim(-3,3)

for i in range(0,len(Omeg_range)):
  from_ = int(round( 0.99 * Q[i].shape[1] , 0 )) 
  to = Q[i].shape[1] 
  dataLength = len(range(from_,to))  
  Omeg = Omeg_range[i] 

  ax1.plot( Q[i][0,from_:to], Q[i][1,from_:to], Omeg, 'k' )

  rH = ( Q[i][0,:]**2 + Q[i][1,:]**2 )**0.5 
  ax2.plot( Omeg*n.ones(dataLength) , rH[from_:to] , "k.", Omeg , max(rH[from_:to]), "bo")  #|Nt:(1, 104) and (104,) does not match but (104,1) and (104,) does match


  TPeriod = 2*n.pi/Omeg 
  k = 1  
  r = 1  
  IPeriod = []  
  fromP = int(round( 0.60 * Q[i].shape[1] , 0 )) 
  toP = -1 
  while k*TPeriod <= T[i][-1] :
    if k*TPeriod > T[i][fromP] and k*TPeriod <= T[i][toP]:
      IPeriod.append( n.argmin( abs( T[i] - k*TPeriod ) ) )
      r = r+1 
    k = k+1 
  
  ax1.plot( Q[i][0,IPeriod], Q[i][1,IPeriod], Omeg,
            ".",linewidth=1,color=colour)
  ax2.plot( Omeg*n.ones(len(IPeriod)) , rH[IPeriod], ".",linewidth=1,color=colour)

p.savefig('myplot.png')
p.show()
#|NNN. 3D&2D 
stop_plot = time.time()
time_plot = stop_plot - start_plot
print("time_plot is {}".format(time_plot) )
#|FFF. PLOTTING
