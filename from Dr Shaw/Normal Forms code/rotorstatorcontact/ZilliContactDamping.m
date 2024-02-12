function [ nlf ] = ZilliContactDamping( q , qdot , rc, cs  )
%ZilliContactDamping: Evaluates nonlinear force (if any) due to damping while in contact
%    -i.e. the nonconservative part
% q is the u,v displacement of the rotor (column vector)
% qdot the velocities of u , v
% rc the clearance
% cs the stator damping
% can work on on a time history of column vectors

rsq= sum(q.*q); %displacement magnitude squared
rdot_over_r = sum(q.*qdot)./rsq; %radial displacement vel over r

flt=rsq>=rc^2;

nlf=zeros(size(q));

if any(flt)
    nlf(:,flt) = -( cs*rdot_over_r(flt) ).*q(:,flt);
end
% if rsq<rc^2;
%     nlf=[ 0 ;
%           0 ];
% else
%     nlf =  -( cs*rdot_over_r )*q;
% end


end

