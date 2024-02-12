%G_ and K_ are the conservative gyoscopic and stiffness matrices
%Kc_ is the stiffness contribution due to damping
function [G_ ,  K_  , Kc_] = rotating_matrices(M,G1,C,K,K1,Omega)
ndof_ = size(M,1);  %underscore to remind us/me that this number excludes zeroed dofs
jdiag=ones(1,ndof_-1);               %
jdiag(2:2:end)=0;                    %
J = diag(-jdiag,1) + diag(jdiag,-1) ; clear jdiag; % assemble skew symmetric matrix J



G=G1*Omega;
G_ = 2*Omega*M*J+G ;


K_ =  K-Omega^2*M+Omega*G*J; 

Kc_ = Omega *C*J +K1*Omega;

end

