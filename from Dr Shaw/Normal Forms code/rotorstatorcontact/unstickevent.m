% detects the end of a period of 'stuck' motion following 
% A.B Nordmark and P.T Piiroinen.
% Simulation and stability analysis of impacting systems with complete chattering.
% Nonlinear dynamics, 58(1-2):85–106, 2009.
function [value,isterminal,direction] =unstickevent( y , AA_ , b_ , M , d, delta )

 Fx =  AA_*y +  [ 0; 0; M\b_ ];
 value =   [ -(-y(1)*y(3) -y(2)*y(4))*y(1)/delta^3-y(3)/delta;
                   -(-y(1)*y(3) -y(2)*y(4))*y(2)/delta^3-y(4)/delta;
                   -y(1)/delta;
                   -y(2)/delta; ]'*Fx;

 
isterminal=true;
direction =1;
end