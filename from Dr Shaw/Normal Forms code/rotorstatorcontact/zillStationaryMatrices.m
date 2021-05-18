%Define system matrices for a Zilli system in stationary coordinates
function [M,G,C,K,Kom,Cint] = zillStationaryMatrices( m, Jp, zeta, k , zeta_internal)

if nargin < 5 
    zeta_internal=0;
end

M = [ m 0;
    0 m ];  

G = [ 0 Jp;
    -Jp  0 ]; %Gyroscopic matrix

% C = [ 2*zeta 0;
%     0   2*zeta  ];
C= 2*zeta  *eye(2);
Cint= 2*   zeta_internal  *eye(2);

K = [ k 0 ;
    0 k ]; %stiffness=identity due to nondim
 
Kom= [ 0 , zeta_internal;
       -zeta_internal,  0 ];
