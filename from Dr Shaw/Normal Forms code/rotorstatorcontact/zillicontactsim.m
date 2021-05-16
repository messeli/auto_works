function [time,response,impacts,endstate] = zillicontactsim(model,Omega,tspan,y0)
%
%  function  timesim.m
%
%   [time,response,speed,impacts,endstate] = zillicontactsim(model,alpha,tspan,y0)
%
% Calculates the response  using time integration
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
%    model        is the rotor model.
%
%    Omega        rotor speed - note simulation is rotating coords, where
%    synchronous force direction is always assumed to act along postive
%    x-axis
%
%    tspan        the time span vector for ode45
%
%   y0            starting state variable
%

if nargin < 4 % no supplied start state, assume at rest
    y0 = zeros(1, 4 );
end


% rotating coordinate simulation
%rotating sys matrices
M=model.M;
C=model.C;
[G_,K_,Kc_] = zillRotatingMatrices( M,model.G,C,model.K,Omega );

%state space matrix
AA_ = [ zeros(2)    eye(2)     ;
    -M\(K_+Kc_)     -M\( C +  Omega*G_) ];

%synchronous forcing vector
b_ = -Omega^2*model.m*model.epsln*[1;0];

%simulation functions (linear and nonlinear/contacting)
fls = @(t,y) AA_*y +  [ 0; 0; M\b_ ];
fnls = @(t,y) fls(t,y) + ...
    [0;  0 ; M\zillNx_tsim(  y,   model.beta , model.cs  ) ];

%determine initial contact state
incontact = sum(y0(1:2).^2) >= (model.rc)^2;

%set up time span and steps
t0=tspan(1);
tend=tspan(2);
nlstep=0.00001*2*pi/(Omega*8);
lstep=0.00001*2*pi/(Omega*8);

%buffer for output
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
while t0<tend
    if incontact==false
        opsl=odeset( ops, 'InitialStep',lstep);
        
        [t,y ] = ode45(fls,[t0 tend],y0,opsl);
        
        %check that we have genuinely been out of contact
        r=sqrt( sum( y( : , 1:2 )'.^2 ) );
        if mean(r) >= model.rc
            incontact=true;
            continue  %repeat the whole chunk with same ICs but contact assumed
        else
            lstep=t(end)-t(end-1); %update time step, rest of loop advances simulation
        end
        
    else
        opsnl=odeset( ops, 'InitialStep',nlstep);
        
        [t,y] = ode45(fnls,[t0 tend],y0,opsnl);
        
        nlstep=t(end)-t(end-1);
    end
    
    %set up next stage
    t0=t(end);
    y0=y(end,:)';
    
    %add chunk to history buffer
    if t0<tspan(2)
        chunklength=size(t,1)-1;%avoid duplicating time point at beginning and end of chunks
    else
        chunklength=size(t,1);
    end
    dataendpos=datastartpos+chunklength-1;
    if dataendpos>CurrBuffersize
        %increase buffer size
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
    
    %evaluate new contact state based on velocity
    if y0(3:4)'*y0(1:2) >= 0
        incontact = true;
    else
        incontact = false;
    end
    
    
end
time = timehist(1:dataendpos).';
impacts = impacts(1:dataendpos).';
npts = dataendpos;
response = yhist(1:dataendpos,:).';
endstate = y0; %for possible restarts

end






