function [ nlf ] = ZilliContactStiffness( q , rc, ks  )
%ZilliContactStiffness: Evaluates nonlinear force (if any) due to contact
%stiffness -i.e. the conservative part
% q is the u,v displacement of the rotor (column vector)
% rc the clearance
% ks the stator stiffness
% works on a time histoyry of col vectors

r= sqrt(sum(q.*q)); %displacement magnitude

flt=r>=rc;

nlf=zeros(size(q));

if any(flt)
    nlf(:,flt)= -( ks*(r(flt)-rc)./r(flt) ).*q(:,flt);
end


% if r<rc;
%     nlf=[ 0 ;
%           0 ];
% else
%     nlf =  -( ks*(r-rc)/r )*q;
% end


end

