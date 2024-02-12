function zmdl = CreateZilliMdl( Meff, Jp, zeta, k , m, epsln, beta, rc , cs , zeta_internal )
% CreateZilliMdl( Meff, Jp, zeta, k , m, epsln, beta, rc )
%   Creates a structure comprising all of the model parameters and matrices
%   for an overhung contacting rotor (a 'Zilli system'). 
% Meff - effective mass of rotor
% Jp  -  polar moment of inertia
% zeta -  damping ratio
% k  - effective stiffness
% m -  out of balance mass
% epsln -  out of balance distance
% beta- contact stiffness (ratio to linear stiffness)
% rc -contact clearance
% cs - stator damping
% zeta_internal - internal damping coefficient

zmdl.Meff= Meff ;
zmdl.Jp= Jp ; 
zmdl.zeta= zeta ; 
zmdl.k = k ;
zmdl.m= m ; 
zmdl.epsln= epsln ; 
zmdl.f = epsln*m;
zmdl.beta=  beta ;
zmdl.rc= rc ;

if nargin > 8
    zmdl.cs = cs;
end
    
if nargin < 10
    zeta_internal=0;
end
zmdl.zeta_internal=zeta_internal;

[M,G,C,K,Kom,Cint] = zillStationaryMatrices( Meff, Jp, zeta, k , zeta_internal );%

zmdl.M= M ;
zmdl.G= G ;
zmdl.C= C ;
zmdl.K= K ;
zmdl.Kom= Kom ;
zmdl.Cint= Cint ;

end

