Source: 
1) https://edbennett.github.io/high-performance-python/05-profiling/index.html
2) A past tutorial but needed: https://edbennett.github.io/SCW-tutorial/06-optimising-for-parallel-processing/ 

>My Notes after the course:

The 4th (last) day (11 Mar 2021) of the HPC-Python Supercomputing-Wales (SCW) tutorial was quite usefull for me, with Michele. 
He helped so much with the identification of the wrongs in my implementation of GNU-parallel and Pathos: they are now functioning.

Indeed, I could not get the same level of help from Viladimir, I think that he was not prepared to the course content: he tried to help me with methods that were not taught in the course material (as least I thought so). 

Ed was quite helpfull too. However, I think he teaches a bit fast, faster than my brain quite likes. 

Since this morning, I have been reviewing the course material. 

Today we have done 
- cProfile the non-parallel script (see the commented part in pathos sript, indeed we used this pathos script with its the non-parallel version to cProfile )
- GNU Parallel version was fixed. Now works well in both my local Ubuntu and remote Sunbird.
- Pathos script was fixed. Now works on Win10, Ubuntu Virtual machine, Sunbird. 

>Notes on parallelisation:
with the parameters: 
  qn_ = [0;0.6;-1.2460;0];
  tend = 1000 ;%|Switch tdir setting - FAIL:damping, See Mathworks question.
  tol = 1e-6; 
  Omeg_range_txt = "5.01:0.1:7.01" ;%"4.81" "0.01:0.1:7.01" ;%|Cld use "0.01:0.01:7.01" |"2.91"

to be the same in Win10-Matlab, UbuVM-Python, Sunbird-Python, Win10-Python;
with 4 processing threads in parallel;
it took the simul and plot to complete about:
-  6 sec and 4   sec for Win10-Matlab (22 sec and 4 sec without parfor)
- 28 sec and 3.6 sec for both GNU Parallel and Pathos on Ubuntu-Python
- 26 sec and 3.6 sec for GNU Parallel + Numba.jit on Ubuntu-Python
- 32 sec and 3   sec for Pathos on Win10-Python 
- 38 sec and 1   sec for Pathos on Sunbird-Python
- 21 for Numba+Pathos on Win10-Python (66 sec without pathos with numba)
- (102 sec for no no parallel=noPathos+noGNUParallel and no numba.jit)
(Take the plotting times with a grain of salt: I might have done p.show() and closed in 4 sec etc.)
Also;
- python's scipy.integrate.solve_ivp generated 20000 - 30000 data points 
- Matlab's ode45 generated 120000 - 170000 data points (could explain the higher plotting times BN)
Also;
- Michele said that the implementation of the code on Fortran/C would not make much difference since already the scipy routines use FOrtran codes and wrap python punctions (->BN).

BOTTOMLINE: MATLAB ode45 is good enough (better even w/o parfor than any inplementation w/ Python), only about 6 times more data points are generated compared to scipy.integrate.solve_ivp. 




