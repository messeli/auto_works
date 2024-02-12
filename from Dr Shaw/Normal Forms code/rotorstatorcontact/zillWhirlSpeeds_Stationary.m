% calculate linear whirl speeds of an overhung rotor in stationary
% coordinates
function [ omfw , ombw ] = zillWhirlSpeeds_Stationary( M,G,C,K,Omvals,Kom,Cint )
  % M,G,C,K - mass, gyroscopic, damping and stiffness matrices respectively
  % Omvals - a 1d array of values of shaft speed (rad/s) to evalaute whirls
  % speeds at
  % Kom - this gives a shaft speed dependant 'stiffness-like' response,
  % that can arise when there is internal damping on the shaft ( usually a
  % skew symmetric matrix which can be destabilising)
  % Cint - damping that is due to internal material effects in the shaft,
  % rather than external fluid interaction. 

if nargin<6
    Kom=zeros(2);
    Cint=zeros(2);
end

ombw=zeros(size(Omvals));
omfw=ombw;

for ii=1:length(Omvals)
    Omega=Omvals(ii);
    
    %stationary whirl calc
    AA = [ zeros(2)    eye(2)     ;
        -M\( K + Omega*Kom )     -M\( C + Cint + Omega*G) ];
    [ V, D]=eig(AA);
    %pick the pair of modes where sign of om indicates whirling dir
    flt = imag( V(2,:)./V(1,:) )<0; %|pick [1,-j] mode & judge dir by d sign f corrspJ eigval: BW if - , FW if +, and greater one is FW_mode 
    s=diag(D);
    oms= imag(s(flt));
    ombw(ii)=min(oms);
    omfw(ii)=max(oms);
end