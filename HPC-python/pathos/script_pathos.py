from pathos.pools import ParallelPool
import numpy as n
import matplotlib.pyplot as p
from mpl_toolkits.mplot3d import Axes3D 
import func_ode45
import time 

start_sim = time.time()
tend = 1000 
tol = 1e-6
qn_ = n.array([0,0.6,-1.2460,0]).reshape(4,1)

Omeg_range = n.arange(5.1,7.2,0.1)  
N = len(Omeg_range) 
print("Omeg_range length is {}".format(N) )

pool = ParallelPool(nodes=4)

def one_Omeg_sim(Omeg,qn_,tend,tol):
  import numpy as n
  import func_ode45

  T_X2rot = n.array([[1,0,0,0],[0,1,0,0],[0,-Omeg,1,0],[Omeg,0,0,1]]) #|q_X = T_X2rot . q_rot
  T_phi2X = n.array([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]]) #|q_phi = T_phi2X . q_X
  T_X2rot_INV = n.linalg.inv(T_X2rot)  
  T_phi2X_INV = n.linalg.inv(T_phi2X) 
  qn = T_X2rot_INV @ T_phi2X_INV @ qn_

  (T_new,Q_new) = func_ode45.func_ode45_v1(Omeg,qn,tend,tol)
  
  return (T_new,Q_new)

# X = [ one_Omeg_sim(*args) for args  in zip(Omeg_range,[qn_]*N,[tend]*N,[tol]*N) ] #|Non parallel version, used for profiling with cProfile
X = pool.map(one_Omeg_sim,Omeg_range,[qn_]*N,[tend]*N,[tol]*N)

T = [ T_new for T_new,Q_new in X]
Q = [ Q_new for T_new,Q_new in X]

end_sim = time.time() 
time_sim = end_sim  - start_sim
print("time_sim is {}".format(time_sim) )



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
start_plot = time.time()

colour = (1.0,0.7,0.0) 
fig = p.figure(figsize=(10.0,3.0)) #figure
ax1 = fig.add_subplot(121,projection='3d')  #|D ",projection='3d'" part adds Axes3D type axes
ax2 = fig.add_subplot(1,2,2)  # ax2 = subplot(122) ; hold(ax2,'on') ; grid(ax2,'on') ; axis(ax2, 'square')

ax1.set_xlabel("u")  
ax1.set_ylabel("v")  
ax1.set_zlabel("Omeg")  
ax1.set_xlim(-3,3)  
ax1.set_ylim(-3,3)  

ax2.set_xlabel("Omeg"     )  
ax2.set_ylabel("Amplitude")  
ax1.set_ylim(-3,3) 


for i in range(0,len(Omeg_range)):
  # with open(f"res/T-Q-{Omeg}","r"):
  #   lines = f.readlines()
  #   T = []
  #   Q = []
  #   for line in lines: 
  #     T.append = line.split(",")[0]
  #     Q.append = line.split(",")[1:4].reshape(4,1)

  from_ = int(round( 0.99 * Q[i].shape[1] , 0 )) 
  to = Q[i].shape[1] 
  dataLength = len(range(from_,to))  
  Omeg = Omeg_range[i] 

  ax1.plot( Q[i][0,from_:to], Q[i][1,from_:to], Omeg, 'k' )

  rH = ( Q[i][0,:]**2 + Q[i][1,:]**2 )**0.5 
  ax2.plot( Omeg*n.ones(dataLength) , rH[from_:to] , "k.", Omeg , max(rH[from_:to]), "bo")  #|Nt:(1, 104) and (104,) does not match but (104,1) and (104,) does match

  #|PPP> STROBOSCOPIC POINCARE SECTION 
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
    k = k+1   #|Period's number.
  
  ax1.plot( Q[i][0,IPeriod], Q[i][1,IPeriod], Omeg,
            ".",linewidth=1,color=colour)
  ax2.plot( Omeg*n.ones(len(IPeriod)) , rH[IPeriod], ".",linewidth=1,color=colour)
  #|PPP. POINCARE

p.savefig('myplot.png')
p.show()

stop_plot = time.time()
time_plot = stop_plot - start_plot
print("time_plot is {}".format(time_plot) )