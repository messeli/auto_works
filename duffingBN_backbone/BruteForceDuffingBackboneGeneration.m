clc
clear
%close all 

N = 100 ;
options = odeset("RelTol",1e-10) ;
for k = 1:N
% u0 = -3+6*rand(2,1) ;% This is working too, except for the colour grad.
u0 = [k*0.05;0] ;%init cond ; amp 20th s 1.0.
[T{k},U{k}] = ode45(@func,[0 100],u0,options) ;
[~,peakLoc] = findpeaks(U{k}(:,1),T{k}) ; 
PER{k} = mean(diff(peakLoc));
AMP{k} = max(U{k}(:,1)) ;
W{k} = 2*pi./PER{k} ;
end

figure
ax1 = subplot(121),hold(ax1,"on")
ax2 = subplot(122),hold(ax2,"on")
axis(ax2,[0 3 0 3])
for k = 1:N
  colour = [k/N (N-k)/N 0] ;
  plot(ax1,U{k}(1,1),U{k}(1,2),"o","color",colour)
  plot(ax1,U{k}(:,1),U{k}(:,2),"color",colour)
  plot(ax2,W{k},AMP{k},".","color",colour) 
end

function f = func(t,u)
  wn = 1;
  k3_m = 0.5 ; 
  f = [u(2);-wn^2*u(1) - k3_m*u(1)^3] ;
end

% Credit of period finding solution is due to 
% Kaushik Lakshminarasimhan's answer @
% https://uk.mathworks.com/matlabcentral/answers/365236-how-to-calculate-a-period-of-a-wave-when-you-don-t-have-the-equation-of-the-wave#answer_289555


