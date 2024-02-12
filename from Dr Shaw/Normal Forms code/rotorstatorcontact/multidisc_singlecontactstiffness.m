function [ nlf ] = multidisc_singlecontactstiffness( q , rc, ks , statordofs )
%multidisc_singlecontactstiffness: Evaluates nonlinear force (if any) due to contact
%stiffness -i.e. the conservative part
% q is the u,v displacement of the rotor (column vector)
% rc the clearance
% ks the stator stiffness
% statordofs are the x and y dofs at the stator
% works on a time histoyry of col vectors

r= sqrt(sum(q(statordofs,:).*q(statordofs,:))); %displacement magnitude - note this doesn't work for multiple stators!

flt=r>=rc;

nlf=zeros(size(q));

if any(flt)
    nlf(statordofs,flt)= -( ks*(r(flt)-rc)./r(flt) ).*q(statordofs,flt);
end


% if r<rc;
%     nlf=[ 0 ;
%           0 ];
% else
%     nlf =  -( ks*(r-rc)/r )*q;
% end


end

