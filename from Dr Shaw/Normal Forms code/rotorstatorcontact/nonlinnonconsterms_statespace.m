function [ ny ] = nonlinnonconsterms_statespace( y , M , Ndof , nq0 , Kc_ , C , nqzeta , b_ )
%returns a derivative vector (in state space form) comprising all nonlinear
%and nonconserative terms for generic mdof rotor in rotating system: 
% y - current state space vector
% The folowing are given for the usual 2nd order representation in rotating
% system: 
%   M - mass matrix
%   ndof - size (of 2nd order system)
%   nq0 - (function handle)nonlinear function with  conservative terms
%  Kc_ - stiffness like matrix arising from the transformation of damping
%  terms 
% C - damping matrix
% nqzeta -(function handle)  nonconervative nonlinear function
% b_ - out of balance forcing vector (constant in rotating system)
 ny =   [   zeros(Ndof,size(y,2));  M\nq0(y(1:Ndof,:)) ]  ...
        + [ zeros(Ndof,2*Ndof) ; -M\Kc_,    -M\C    ] * y ...
        + [ zeros(Ndof,size(y,2));  M\nqzeta(y(1:Ndof,:),y((Ndof+1):end,:)) ] ...
        + [ zeros(Ndof,size(y,2));  M\b_*ones(1,size(y,2)) ];

end

