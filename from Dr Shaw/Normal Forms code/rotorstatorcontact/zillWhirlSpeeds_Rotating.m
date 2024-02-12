% evaluate whirl speeds in a corotating coordinate system
function [ omfw_ , ombw_ ] = zillWhirlSpeeds_Rotating( M,G,C,K,Omvals,Cint )
  % M,G,C,K - mass, gyroscopic, damping and stiffness matrices respectively
  % Omvals - a 1d array of values of shaft speed (rad/s) to evalaute whirls
  % speeds at.
  % Cint - damping that is due to internal material effects in the shaft,
  % rather than external fluid interaction. 

if nargin< 6
    Cint=zeros(2);
end

ombw_=zeros(size(Omvals));
omfw_=ombw_;

for ii=1:length(Omvals)
    Omega=Omvals(ii);
    
    %get rotating system matrices for current Omega
    [G_,K_,Kc_] = zillRotatingMatrices( M,G,C,K,Omega );
    
    %state space matrix
    AA_ = [ zeros(2)    eye(2)     ;
        -M\(K_+Kc_)     -M\( C +  Cint + Omega*G_) ];
    
    %extract roots
    [ V, D]=eig(AA_);
    %pick the pair of modes where sign of om indicates whirling dir
    flt = imag( V(2,:)./V(1,:) )<0;%|pick [1,-j] mode & judge dir by d sign f corrspJ eigval: BW if - , FW if +, and greater one is FW_mode evenif it is BW looking in rot frame
    s=diag(D);
    oms= imag(s(flt));
    ombw_(ii)=min(oms);
    omfw_(ii)=max(oms);
end