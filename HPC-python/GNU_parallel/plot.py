import numpy as n 
import matplotlib.pyplot as p
from mpl_toolkits.mplot3d import Axes3D 
import time
from os import listdir

#|___PREP>
#| READ result file names
resList = listdir("res/")
N = len(resList)
Omegs = [None] * N
for i in range(N):
  Omegs[i] = float(resList[i][4:])

#| READ T&Q>
T = [None] * N
Q = [None] * N

for i in range(N):

  with open(f"res/T-Q-{Omegs[i]}","r") as db:
    lines = db.readlines()  
    
    T_new = n.empty((len(lines),))
    Q_new = n.empty((4,len(lines)))

    for j,line in enumerate(lines):

      line_content = line.split(",")
      T_new[j]   = line_content[0]
      Q_new[:,j] = line_content[1:5]

  T[i] = T_new
  Q[i] = Q_new
#|___PREP.


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


for i in range(0,N):  # N simulations

  from_ = int(round( 0.99 * Q[i].shape[1] , 0 )) 
  to = Q[i].shape[1] 
  dataLength = len(range(from_,to))  
  Omeg = Omegs[i]

  ax1.plot( Q[i][0,from_:to], Q[i][1,from_:to], Omeg, 'k' )

  rH = ( Q[i][0,:]**2 + Q[i][1,:]**2 )**0.5 
  ax2.plot( Omeg*n.ones(dataLength) , rH[from_:to] , "k.", Omeg , max(rH[from_:to]), "bo")  #|Nt:(1, 104) and (104,) does not match but (104,1) and (104,) does match


  #|___POINCARE> STROBOSCOPIC POINCARE SECTION 
  TPeriod = 2*n.pi/Omeg 
  k = 1  
  r = 1  
  IPeriod = []  
  fromP = int(round( 0.60 * Q[i].shape[1] , 0 )) 
  toP = -1
  while k*TPeriod <= T[i][-1] :
    if k*TPeriod > T[i][fromP] and k*TPeriod <= T[i][toP]:
      IPeriod.append( n.argmin( abs( T[i] - k*TPeriod ) ) )
      r = r+1 #|IPeriod index
    k = k+1 #|Period's number.
  
  ax1.plot( Q[i][0,IPeriod], Q[i][1,IPeriod], Omeg,
            ".",linewidth=1,color=colour)
  ax2.plot( Omeg*n.ones(len(IPeriod)) , rH[IPeriod], ".",linewidth=1,color=colour)

  #|___POINCARE.
p.savefig('myplot.png')
# p.show()

stop_plot = time.time()
time_plot = stop_plot - start_plot
print("time_plot is {}".format(time_plot) )
