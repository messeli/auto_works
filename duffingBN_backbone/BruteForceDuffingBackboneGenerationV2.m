clc
clear
%close all 

N = 100*8 ;
options = odeset("RelTol",1e-7,"abstol",1e-9) ;
for k = 1:N
% u0 = -3+6*rand(2,1) ;% This is working too, except for the colour grad.
u0 = [k*0.01;0] ;%init cond ; amp 20th s 1.0.
[T,U{k}] = ode45(@func,[0 200],u0,options) ;
[~,peakTimes] = findpeaks(U{k}(:,1), T) ;
PER= mean(diff(peakTimes));
AMP(k) = max(U{k}(:,1)) ;
W(k) = 2*pi./PER ;
end

figure
ax1 = subplot(121),hold(ax1,"on")
ax2 = subplot(122),hold(ax2,"on")
set(ax2,"xdir", "reverse") ;
%axis(ax2,[0 3 0 3])
for k = 1:N
  colour = [k/N (N-k)/N 0] ;
  plot(ax1,U{k}(1,1),U{k}(1,2),"o","color",colour)
  plot(ax1,U{k}(:,1),U{k}(:,2),"color",colour)
  plot(ax2,2*pi/W(k),AMP(k),".","color",colour) 
  
end
plot(ax2,  (2*pi./W),smoothdata(AMP),"b")
plot(ax2,2*pi./W,AMP,"k")


%%% Niche!! Also get the single harmonic approximate analytical one 
%%% (As eq.1.22 of b2015WaggNeild)
wn = 1 ; k3_m = 0.5 ;
AMP = 0:0.01:8 ;
W = wn*sqrt(1 + 3*k3_m*AMP.^2/4/wn^2) ; 
plot(ax2,2*pi./W,AMP,"r")
legend(ax2,"numerical sim","num sim smoothed","One harmonic approximate")
%%%


function f = func(t,u)
  wn = 1;
  k3_m = 0.5 ; 
  f = [u(2);-wn^2*u(1) - k3_m*u(1)^3] ;
end

% Credit of period finding solution is due to 
% Kaushik Lakshminarasimhan's answer @
% https://uk.mathworks.com/matlabcentral/answers/365236-how-to-calculate-a-period-of-a-wave-when-you-don-t-have-the-equation-of-the-wave#answer_289555


