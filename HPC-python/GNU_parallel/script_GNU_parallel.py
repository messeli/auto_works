#| HOW TO: 
#| source: https://www.gnu.org/software/parallel/parallel_tutorial.html
#|-Method 1: Run in terminal the command below,
# parallel python3 script_GNU_parallel.py ::: 5.01 5.02
#|-Method 2: Run write-Omeg_range.py to gen Omeg_range.txt, dn run in terminal
# parallel python3 script_GNU_parallel.py :::: Omeg_range.txt
import numpy as n
from argparse import ArgumentParser
import func_ode45
import time 

parser = ArgumentParser()

parser.add_argument("Omeg")
args = parser.parse_args()
Omeg = float(args.Omeg)

start_sim = time.time()
qn_ = n.array([0,0.6,-1.2460,0]).reshape(4,1)

tend = 1000 # 3000
tol = 1e-6 # 1e-7

T_X2rot = n.array([[1,0,0,0],[0,1,0,0],[0,-Omeg,1,0],[Omeg,0,0,1]]) #|q_X = T_X2rot . q_rot
T_phi2X = n.array([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]]) #|q_phi = T_phi2X . q_X
T_X2rot_INV = n.linalg.inv(T_X2rot)  
T_phi2X_INV = n.linalg.inv(T_phi2X) 
qn = T_X2rot_INV @ T_phi2X_INV @ qn_

(T_new,Q_new) = func_ode45.func_ode45_v1(Omeg,qn,tend,tol) 

with open(f"res/T-Q-{Omeg}","w") as d:
  for i in range(len(T_new)): #OR Q.shape[1] #col count # time data count
    d.writelines(f"{T_new[i]},{Q_new[0,i]},{Q_new[1,i]},{Q_new[2,i]},{Q_new[3,1]}\n")  


end_sim = time.time() 
time_sim = start_sim - end_sim 
print("time_sim is {}".format(time_sim) )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
