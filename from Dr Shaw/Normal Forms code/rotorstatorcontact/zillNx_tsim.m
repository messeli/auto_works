function [ Q  ] = zillNx_tsim( x ,   beta , cs)
% zillNx - nonlinear function for the Zilli system

if nargin<3
    cs=0;
end

r = sqrt(x(1:2)'*x(1:2));
Q  = beta* x(1:2) * ( 1/r  - 1 );

if length(x)>2
    %x includes velocity terms, we can do damping
    rdot = ( x(1)*x(3) + x(2)*x(4) ) /r ;
    Q = Q - cs*rdot*x(1:2)/r;
end


end

