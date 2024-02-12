function [time,response,speed,impacts,endstate,qhist] = statorcontactsim(model,alpha,tspan,q0)
%
%  function  timesim.m
%
%   [time, response, speed,impacts,endstate] = statorcontactsim(model, alpha, tspan,nr,q0)
%
% Calculates the response  using time integration
%
%    time         is time vector
%
%    response     is the response output and is a 2 dimensional array
%                 The indices are DoFs and time
%
%    speed        gives the instantaneous rotor speed during runup
%
%    impacts      boolean values cooressponding to time series,
%                   is true whenever the state of rotor contact changes
%
%    model        is the rotor model. It must have reduced modal properties
%    calculated with reduce_rotor_model first. 
%    
%    alpha        [a2 a1 a0], rotor angle = a2*t^2 + a1*t + a0
%
%    tspan        the time span vector for ode45
%
%   q0            starting state variable in form  [ q qdot ] in reduced
%   modal coords
%
%   qhist         entire history in modal coords, with velocites as well
%   (how I should have stored data to begin with!)


% NOTE - NO FLUID BEARINGS FOR THIS SCRIPT
% NOTE - NO SHAFT DAMPING

nr = model.nr;
if nargin < 4 % no supplied start state, assume at rest
    q0 = zeros(1, (model.nr)*2 );
end

Node_Def = model.node;
Force_Def = model.force;
Bearing_Def = model.bearing;
Stator_Def = model.stator;  % assume that stators have been defined
% - no point using contact if not!

% check to make sure there are no speed dependent bearings
nbearing = size(Bearing_Def,1);
const_bearing = 1;
for ib=1:nbearing
    if Bearing_Def(ib,1) > 6.5
        const_bearing = 0;
    end
end
if ~const_bearing
    disp(['>>>> Error in statorcontactsim.m - only onstant stiffness bearings implemented'])
    time = [];
    response = [];
    speed = [];
    return
end

Tr=model.Tr;
Mr = model.Mr;
Cr =  model.Cr;
C1r =  model.C1r;  % part of C matrix determined by gyroscopics
Kr =  model.Kr;
if isfield(model,'K1r')
K1r =  model.K1r;
end
fr =  model.fr;

A = model.A;
A1 = model.A1;
B = model.B;

zero_dof = model.zero_dof;
dof = 1:(model.ndof);
dof(zero_dof) = [];

%variables for contact detection- AS 6/10/2015
nstators = size(Stator_Def,1);
contact_dofnos = zeros(2, nstators); contact_dofhack=contact_dofnos;
for ii=1:nstators
    % compile a list of the Dofs involved in contact
    contact_dofnos(:,ii) = (Stator_Def(ii,1)-1)*4+ [1;2];
end
for zdof =zero_dof  %a somewhat inelegant way of adjusting contact dofs for zeroed dofs
    contact_dofhack(contact_dofnos>zdof)= contact_dofhack(contact_dofnos>zdof)+1;
end
contact_dofnos=contact_dofnos-contact_dofhack;
%create a subset of the transformation matrix that gives potential contact
%dofs only
Tr_contactx=Tr(contact_dofnos(1,:),:);
Tr_contacty=Tr(contact_dofnos(2,:),:);
%precalculate as much as possible
clrsq = ( Stator_Def(:,2) .* Stator_Def(:,2) )'; 
clrnc=sqrt(clrsq);
ist=ones(1,nstators);
dr=zeros(1,nstators);
ks = Stator_Def(:,4);
cs = Stator_Def(:,5);

%time stepping is aware of which stators are in contact and which aren't
% this data is stored in a logical array of 0s and 1s
% these 0s and 1s give a binary number to each permutation of contacts
% this number indexes a list of the  step speed for each configuration, so
% that the intial step size is matched to the contact config and hopefully
% always optimal, with minimal reiterations. 
opts = odeset('Events',@contactchange,'RelTol',1e-8,'AbsTol',1e-8);
t0=tspan(1);
maxstep=0.02*(tspan(2)-tspan(1));
ncontactconfigs = 2^nstators;
contactstate=false(nstators,1);
twopowers = 2.^(0:(nstators-1));
meanspeed=alpha(2)+alpha(1)*tspan(2);
initstep=ones(1,ncontactconfigs)*(2*pi/meanspeed)/40; %guess a reasonable step size
if nargin< 4
    q0=zeros(2*nr,1);
    %note if q0 is specified it must have the correct dimensions
end
%determine initial contact state
clrn=contactchange(0,q0 );
contactstate=clrn<=0;

%buffer for output
BlockSize = 1000000;
timehist = zeros(BlockSize,1);
qhist = zeros(BlockSize,nr*2);
impacts=false(BlockSize,1);
datastartpos=1;
CurrBuffersize=BlockSize;
nreallocations=0;
ChangedStatorNo=0;
retrys=0;
while t0 < tspan(2);
    cstateindex=twopowers*contactstate+1;
    Tr_contactxactive=Tr(contact_dofnos(1,contactstate),:);
    Tr_contactyactive=Tr(contact_dofnos(2,contactstate),:);
    ksactive=ks(contactstate);
    clractive=clrnc(contactstate)';
    contactdofsactive = contact_dofnos(:,contactstate);
    
    dphi = 2*alpha(1)*t0 + alpha(2);
    if dphi>0;
        maxstep=( 2*pi/dphi )/4;  %limit maximum time step to one quarter forcing period
    end
    opts=odeset(opts,'MaxStep',maxstep,'InitialStep',initstep(cstateindex));
    [time,q,TE,YE,IE] = ode45(@deriv,[t0 tspan(2)],q0,opts);
    %make a check that the contact state was appropriate for the step just
    %executed - if not run it again with the opposite contact state
    if ChangedStatorNo ~=0
        snodeno=Stator_Def(ChangedStatorNo,1);
        sx = Tr_contactx(ChangedStatorNo,:)*q(:,1:nr)';
        sy= Tr_contacty(ChangedStatorNo,:)*q(:,1:nr)';
        srsq=sx.^2 + sy.^2;
        if sum(srsq) >= clrsq*length(srsq)
            %it appears contact was active
            contact_state_correct= contactstate(ChangedStatorNo)==true;
         
        else
            %it appears contact was  not active
            contact_state_correct= contactstate(ChangedStatorNo)==false;
       end
        if contact_state_correct==false
            warning('An incorrect contact state encountered, initial RPM: %f, contact state flag  was: %i, radial vel: %e,  retrying.. ',...
                alpha(2)*60/(2*pi),contactstate(ChangedStatorNo),rdottimesr(ChangedStatorNo))
            retrys=retrys+1;
            if retrys> 1
                %should never get here- neither options works so assume
                %contact active, ie. flip and try again if it isn't
                warning('Could not resolve contact state! Assuming contact active. ')
                if contactstate(ChangedStatorNo)==true
                    retrys =0;
                else
                    contactstate(ChangedStatorNo) = true;
                    continue
                end
            else
                %flip it try again
                contactstate(ChangedStatorNo) = ~contactstate(ChangedStatorNo);
                continue
            end
        else
            retrys=0; %everything seems ok, reset the retry counter and carry on
        end
    end
    if isempty(IE)
        ChangedStatorNo=0; %will occur at end of simulation
    else
%         if length(IE)>1
%             warning('Length IE > 1.')
%         end
        ChangedStatorNo=IE(end);%sometimes 2 events are reported, one at the start
        % and one at the end, so this gets the right one
        
        % update flags for stators in contact:      
        %obtain radial velocity at affected stator, and use to determine
        %whether contact is occuring in the next segment 
        rdottimesr = contactradialvel(   YE(end,:)' )  ;
        contactstate(ChangedStatorNo)= ...
                                rdottimesr(ChangedStatorNo) >= 0 ;
                            
                            
                            
      %%%%%%% BELOW ARE ALL THE NUMERICALLY DODGY WAYS THIS TEST USED TO BE DONE                      
% % %                             %         contactstate(ChangedStatorNo)= xor(contactstate(ChangedStatorNo) ,true);
% % %         contactstate(ChangedStatorNo)= contactchange(TE(end),YE(end,:)' )<=0; % AS - this is more robust than 'flipping' as 
% % %         % it did before (see above) and if in doubt (i.e. clearance==0)
% % %         % assumes contact active so the solution cannot diverge
% % %         
% % % %         if contactstate(ChangedStatorNo)==0 && contactchange(TE(end),YE(end,:)' )<0
% % % %             warning('Contact state inactive but contactchange function negative')
% % % %         end
% % % %         if contactstate(ChangedStatorNo)==1 && contactchange(TE(end),YE(end,:)' )>0
% % % %             warning('Contact state active but contactchange function positive')
% % % %         end

        
    end
    initstep(cstateindex)=time(end)-time(end-1); %save time step for current contact state
    
    %set inital conditions for next batch based on impact model
    t0=time(end);
    q0=q(end,:)'; %default options
%     qstate=q0(1:nr);
%     qvel=q0((nr+1):end);
    if ChangedStatorNo>0
        switch Stator_Def(ChangedStatorNo,3)
            case 1
                %do nothing 
            otherwise
                warning('Unidentified contact type')
        end
    end
    
    %add chunk to history buffer
    if t0<tspan(2)
        chunklength=size(time,1)-1;%avoid duplicating time point at beginning and end of chunks
    else
        chunklength=size(time,1);
    end
    dataendpos=datastartpos+chunklength-1;
    if dataendpos>CurrBuffersize
        %increase buffer size
        extraspaceneeded = dataendpos-CurrBuffersize;
        extraspaceroundedup =ceil(extraspaceneeded/BlockSize)*BlockSize;
        timehist = [ timehist; zeros(extraspaceroundedup,1) ];
        qhist = [ qhist; zeros(extraspaceroundedup,nr*2) ];
        impacts = [ impacts; false(extraspaceroundedup,1) ];
        CurrBuffersize=CurrBuffersize+extraspaceroundedup;
        nreallocations=nreallocations+1;
        if nreallocations>10
            warning('10 reallocations occurred - increase block size.')
            nreallocations=0;
        end
    end
    timehist(datastartpos:dataendpos)=time(1:chunklength);
    qhist(datastartpos:dataendpos,:)=q(1:chunklength,:);
    impacts(dataendpos+1)=true;
    datastartpos=dataendpos+1;
    
end
time = timehist(1:dataendpos).';
impacts = impacts(1:dataendpos).';
npts = dataendpos;
response = zeros(model.ndof,npts);
response(dof,:) = Tr*qhist(1:npts,1:nr).';
endstate = q0; %for possible restarts
qhist = qhist(1:npts,:);

speed = 2*alpha(1)*time + alpha(2)*ones(size(time));

    function [qdot] = deriv(t,q )
        phi = alpha(1)*t*t + alpha(2)*t + alpha(3);
        d_phi = 2*alpha(1)*t + alpha(2);
        dd_phi = 2*alpha(1);
        qdot = A*q + A1*q*d_phi + real( B*(d_phi*d_phi-1i*dd_phi)*exp(1i*phi) );
        
        %add contact forces
        if ~isempty(contactdofsactive)
            %stiffness
            xyc  =[ Tr_contactxactive*q(1:nr),    Tr_contactyactive*q(1:nr)  ]'; %location of contact nodes
            rcsq=sum(xyc.*xyc);
            rc=sqrt(rcsq)'; %distance r from centre
            fk = ( ksactive.*(clractive./rc-1)  )'; %works out as the magnitude of force divided by rc
            Fxyk=[fk.*xyc(1,:) ; fk.*xyc(2,:) ];  % directs 'a' in the right direction and multiplies by rc to get force
            
            %damping
            xyvel = [ Tr_contactxactive*q((nr+1):end),    Tr_contactyactive*q((nr+1):end)  ]'; 
            rvel =sum(   xyc .* xyvel )./(rcsq'); %extra division by rc saves nomrmalising position vectors later
            fc = -cs.*rvel;
            Fxyc=[fc.*xyc(1,:) ; fc.*xyc(2,:) ];
            
            %expand forces applied at DOFs and convert back to modal sys
            Fxyall=zeros(model.ndofz,1);
            Fxyall( contactdofsactive ) =Fxyk+Fxyc;
            qdot( (nr+1):end ) = qdot( (nr+1):end )  + Tr.'*Fxyall;
            
        end
    end

    function [clrnc,isterminal,direction] = contactchange(~,q )
        xy =[ Tr_contactx*q(1:nr),    Tr_contacty*q(1:nr)  ]';
        clrnc = clrsq-sum(xy.*xy); %negative if in contact
        isterminal=ist;
        direction=dr;
    end

    function rdottimesr = contactradialvel(q )
        x= Tr_contactx*q(1:nr);
        y= Tr_contacty*q(1:nr);
        xdot= Tr_contactx*q((nr+1):end);
        ydot= Tr_contacty*q((nr+1):end);
        rdottimesr = x.*xdot + y.*ydot ; 
     
    end
end



    
        
        
