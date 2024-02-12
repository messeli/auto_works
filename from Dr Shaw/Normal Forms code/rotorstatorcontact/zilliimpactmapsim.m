function [time,response,impacts,endstate] = zilliimpactmapsim(model,Omega,tspan,y0,dt0)
%
%  function  
%
%   [time,response,speed,impacts,endstate] = zilliimpactmapsim(model,alpha,tspan,y0)
%
% Calculates the response of an overhung rotor with a snubber ring  using 
% time integration and rigid impact events to model contact. ie. contacts
% are instantaneous, with post impact velocity calculated with  a
% coeffecicient of restitution applied to the normal velocity. 
% Problems with chattering are handled via Nordmark and Piiroinen's chatter
% mapping technique, see: 
% A.B Nordmark and P.T Piiroinen.
% Simulation and stability analysis of impacting systems with complete chattering.
% Nonlinear dynamics, 58(1-2):85–106, 2009.
%
%    time         is time vector
%
%    response     is the response output and is a 2 dimensional array
%                 The indices are DoFs and frequency
%
%
%    impacts      boolean values cooressponding to time series,
%                   is true whenever the state of rotor contact changes
%
%    model        is the rotor model created with CreateZilliMdl 
%                    (note, it must also have and additional member e specified 
%                    for the coefficient of restitution )
%
%    Omega        rotor speed - note simulation is rotating coords, where
%    synchronous force direction is always assumed to act along postive
%    x-axis
%
%    tspan        the time span vector for ode45
%
%   y0            starting state variable
%
%    dt0          Initial time step

rndid=cputime;  % uses this to uniquely mark a single simulation in outputs, 
                    % useful when debugging a brute force simulation. 

if nargin < 4 % no supplied start state, assume at rest
    y0 = zeros(1, 4 );
end

if nargin < 5
    dt0=0.00001*2*pi/(Omega*8); %initial time step
end

% rotating coordinate simulation
% oobtain rotating system matrices
M=model.M;
C=model.C;
[G_,K_,Kc_] = zillRotatingMatrices( M,model.G,C,model.K,Omega );
if isfield(model,'Cint')
    Cint=model.Cint;
else
    Cint=0*C;
end

% create linear state space matrix
AA_ = [ zeros(2)    eye(2)     ;
    -M\(K_+Kc_)     -M\( C + Cint +  Omega*G_) ];

%synchronous forcing vector
b_ = Omega^2*model.m*model.epsln*[1;0];

%simulation functions
% basic linear system
f_ls = @(t,y) AA_*y +  [ 0; 0; M\b_ ];
% functions needed for chatter mapping and the 'stuck system'
delta = model.rc; % clearance
rx=model.d;%coef of restitution
ax_ = @(t,y) [  -(-y(1)*y(3) -y(2)*y(4))*y(1)/delta^3-y(3)/delta;
                -(-y(1)*y(3) -y(2)*y(4))*y(2)/delta^3-y(4)/delta;
                -y(1)/delta;
                -y(2)/delta; ];   %(eq 3 from Nordmark Piiroinen 2009)
Wx_ = @(t,y) ( (1+rx) / delta ) * [ 0, 0, y(1), y(2) ]'; %(eq 7 (I think) from Nordmark Piiroinen 2009)
% alternative function to use when stuck in sliding contact (eq 10 from Nordmark Piiroinen 2009):
f_stck = @(t,y) f_ls(t,y)  + ( (ax_(t,y)'*f_ls(t,y) )/ (1 +rx )) * Wx_(t,y);

%set up time span and steps
t0=tspan(1);
tend=tspan(2);
nlstep=0.01*2*pi/(Omega*8);
lstep=0.01*2*pi/(Omega*8);

%buffer for output - this avoids multiple reallocations of memory when
%sticking lots of short segments of simulation together (for performance)
BlockSize = 1000000;
timehist = zeros(BlockSize,1);
yhist = zeros(BlockSize,4);
impacts=false(BlockSize,1);
datastartpos=1;
CurrBuffersize=BlockSize;
nreallocations=0;

% run sim
ops=odeset( 'Events', @(t,y) zillEvent(t,y(1:2),model.rc), ...
    'RelTol', 1e-8,'AbsTol',1e-8, 'MaxStep' , 2*pi/(Omega*8) );
opsl=odeset( ops, 'InitialStep',dt0);
opslstk =odeset( ops,'Events', @(t,y) unstickevent( y , AA_ , b_ , M , rx, delta ) );
bStuck=false;  % flag is true if rotor is in constant sliding contact along the stator
while t0<tend
    %run simulation until contact occurs or simulation ends:
    if bStuck==false
        %no contact
        [t,y ] = ode45(f_ls,[t0 tend],y0,opsl);
        opsl=odeset( opsl, 'InitialStep',t(end)-t(end-1));
    else
        %sliding contact
        [t,y ] = ode45(f_stck,[t0 tend],y0,opslstk);
        opslstk =odeset( opslstk, 'InitialStep',t(end)-t(end-1));
        fprintf('Chatter ended at time %f on id %f !\n' , t(end) , rndid );
    end
    
    %set up next stage
    t0=t(end); 
    
    %add new chunk of simulation to history buffer
    if t0<tspan(2)
        chunklength=size(t,1)-1;%avoid duplicating time point at beginning and end of chunks
    else
        chunklength=size(t,1);
    end
    dataendpos=datastartpos+chunklength-1;
    if dataendpos>CurrBuffersize
        %increase buffer size if it is too small
        extraspaceneeded = dataendpos-CurrBuffersize;
        extraspaceroundedup =ceil(extraspaceneeded/BlockSize)*BlockSize;
        timehist = [ timehist; zeros(extraspaceroundedup,1) ];
        yhist = [ yhist; zeros(extraspaceroundedup,4) ];
        impacts = [ impacts; false(extraspaceroundedup,1) ];
        CurrBuffersize=CurrBuffersize+extraspaceroundedup;
        nreallocations=nreallocations+1;
        if nreallocations>10
            warning('10 reallocations occurred - increase block size.')
            nreallocations=0;
        end
    end
    timehist(datastartpos:dataendpos)=t(1:chunklength);
    yhist(datastartpos:dataendpos,:)=y(1:chunklength,:);
    impacts(dataendpos+1)=true;
    datastartpos=dataendpos+1;
    
    %evaluate new initial conditions
    if bStuck==true   %if exiting the sliding condition, simply start a new block
                            % no need for impact map
        y0=y(end,:)';
        bStuck=false;
        continue
    end
    
    %apply impact map
    t0=t(end);
    x=y(end,:)';
    vx = ( - x(1)*x(3) - x(2)*x(4) ) /delta;
    Fx=f_ls(t0,x);
    a1= -x(1)*x(3) -x(2)*x(4) ;
    ax = [ -a1*x(1)/delta^3-x(3)/delta;
        -a1*x(2)/delta^3-x(4)/delta;
        -x(1)/delta;
        -x(2)/delta; ]' * Fx;
    if vx>0  && t(end) < tend
        fprintf('positive vx %e at time %f on id %f !\n', vx , t0, rndid );
    else
        %                 disp('negative vx')
    end
    % test for a chatter starting..
    if abs(vx)< 1e-5 && ax < 0 && t(end) < tend
        % apply chatter mapping
        Wx = ( (1+rx) / delta ) * [ 0, 0, x(1), x(2) ]';
        xstar = x +  ( 1/(1-rx) ) * ( 2 * Fx * rx/ax + Wx ) *vx;%(eq 30 from Nordmark Piiroinen 2009)
        deltstar =  ( 1/(1-rx) ) * ( 2*rx/ax ) * vx;%(eq 31 from Nordmark Piiroinen 2009)
        
        fprintf('Chattering detected at vx %e time %f on id %f !\n', vx , t0 ,rndid);
        
        %append to time sig and set new initial conditions
        t0=t0+deltstar;
        t=[ t ; t0 ];
        y0 = xstar';
        y = [ y ; y0 ];
        
        fprintf('Chattering complete at time %f on id %f !\n' , t0 , rndid );
        
        bStuck=true;
        
    else
        %normal impact
        ypos=y(end,1:2)';
        yvel=y(end,3:4)';
        %project velocity in radial/tangential direction
        rmat = [ ypos(1) ypos(2)
            -ypos(2) ypos(1) ]/delta;
        rtvel = rmat*yvel;
        
        %apply reflection/reduction of radial vel
        rtvel(1) = -rx*rtvel(1);
        y0 = [ ypos ; rmat'*rtvel ]; %transform back to normal coords
        
        bStuck =false;
    end
    
    
    
end
% truncate empty space from the history buffers before returning data.
time = timehist(1:dataendpos).';
impacts = impacts(1:dataendpos).';
response = yhist(1:dataendpos,:).';
endstate = y0; %for possible restarts

end






