function [rModel] = reduce_rotor_model(model,nr,assumeexternaldamping)
%Defines reduced modal model for rotor based on static modes.

% NOTE - NO FLUID BEARINGS FOR THIS SCRIPT
% NOTE - NO SHAFT DAMPING

if nargin<3
    assumeexternaldamping=false;
end

rModel=model; % aim to tag new reduced matrices into model
rModel.nr=nr;

Node_Def = model.node;
Force_Def = model.force;
Bearing_Def = model.bearing;

% check to make sure there are no speed dependent bearings
nbearing = size(Bearing_Def,1);
const_bearing = 1;
for ib=1:nbearing
    if Bearing_Def(ib,1) > 6.5
        const_bearing = 0;
    end
end
if ~const_bearing
    disp(['>>>> Error in runup.m - only onstant stiffness bearings implemented'])
    time = [];
    response = [];
    speed = [];
    return
end


% obtain model of rotor
[M0,C0,C1,K0,K1] = rotormtx(model);
if assumeexternaldamping==true
    K1=zeros(size(K1));
end
[nnode,~] = size(Node_Def);
ndof = 4*nnode;
rModel.ndof=ndof;

% sort out zeroed DoF and determine bearing model
[Mb,Cb,Kb,zero_dof] = bearmtx(model,0.0);
dof = 1:ndof;
dof(zero_dof) = [];

% calculate machine model
M = M0 + Mb;
K = K0 + Kb;
C = C0 + Cb;
%AS fix - reduce size of matrices to exclude fixed DOF
M=M(dof,dof);
K=K(dof,dof);
K1=K1(dof,dof);
C=C(dof,dof);
C1=C1(dof, dof);

% sort out forcing
[nforce,~] = size(Force_Def);
ubforce = zeros(ndof,1);
jot = sqrt(-1);
for iforce = 1:nforce
    if Force_Def(iforce,1) == 1  % unbalance force
        node = Force_Def(iforce,2);
        unbal_mag = Force_Def(iforce,3);
        unbal_phase = Force_Def(iforce,4);
        force_dof = [4*node-3; 4*node-2];
        ubforce(force_dof) = ubforce(force_dof) + unbal_mag*exp(jot*unbal_phase)*[1; -1i];
    end
    if Force_Def(iforce,1) == 2  % unbalance moment
        node = Force_Def(iforce,2);
        unbal_mag = Force_Def(iforce,3);
        unbal_phase = Force_Def(iforce,4);
        force_dof = [4*node-1; 4*node];
        ubforce(force_dof) = ubforce(force_dof) + unbal_mag*exp(jot*unbal_phase)*[1i; 1];
    end
end
ubforce(zero_dof) = [];  % remove force from zeroed dof

% reduce the model
[ndofz,~] = size(M);
nr = round(nr);
% model reduction based on undamped modes
[eigvec,eigval] = eig(K,M);
[eigval,isort] = sort(diag(eigval));
eigvec = eigvec(:,isort);
Tr = eigvec(:,1:nr);
eigmaxr = sqrt(eigval(nr))/(2*pi);
disp(['Maximum frequency of reduced system is ' num2str(eigmaxr) ' Hz'])
% disp(sqrt(eigval(1:nr))/(2*pi))

rModel.Tr=Tr;
rModel.Mr = Tr.'*M*Tr;
rModel.Cr = Tr.'*C*Tr;
rModel.C1r = Tr.'*C1*Tr;
rModel.Kr = Tr.'*K*Tr;
rModel.K1r = Tr.'*K1*Tr;
rModel.fr = Tr.'*ubforce;

Mr = rModel.Mr;
Cr =  rModel.Cr;
C1r =  rModel.C1r;  % part of C matrix determined by gyroscopics
Kr =  rModel.Kr;
K1r =  rModel.K1r;
fr =  rModel.fr;

rModel.ndofz=ndofz;

rModel.A = [zeros(nr,nr) eye(nr,nr); -Mr\Kr -Mr\Cr];
rModel.A1 = [zeros(nr,nr) zeros(nr,nr); -Mr\K1r -Mr\C1r];
rModel.B = [zeros(nr,1); Mr\fr];

rModel.zero_dof= zero_dof;

   
end



    
        
        
