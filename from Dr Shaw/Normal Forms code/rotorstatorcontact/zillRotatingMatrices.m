%Transform linear system matrices to the rotating coordinate system
function [G_,K_,Kc_] = zillRotatingMatrices( M,G,C,K,Omega )
  % M,G,C,K - mass, gyroscopic, damping and stiffness matrices respectively
  % Omega - shaft speed (rad/s) 
  
J=[ 0 -1; 
    1  0 ];  %skew symmetric matrix

% M stays as it is
G_ = 2*M*J + G;
K_ = K - Omega^2*M + Omega^2*G*J;
% C stays as it is, so does internal C
Kc_ = Omega*C*J; %skew symmetric stiffness-like component due to damping
