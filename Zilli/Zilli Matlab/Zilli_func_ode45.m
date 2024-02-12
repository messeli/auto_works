function [T,Q] = Zilli_func_ode45(thetaP,qn,tend,tol,coordinateSystem) 
% [t,q] = Zilli_func_ode45(6, [0.5;0;0;0.5],1000,1e-7,'sta')
% thetaP = 5 ;
% qn = [0.5;0;0;0.5] ;
% tend = 1000 ;
% tol = 1e-7 ;
% coordinateSystem = 'sta' ;

if nargin <= 4 % number of arguments input
    f = @f_sta ;
else
    if     coordinateSystem == "sta", f = @f_sta ;
    elseif coordinateSystem == "rot", f = @f_rot ; 
    end
end
 
theta = 0 ;
isContact = sqrt( (qn(1))^2 + (qn(3))^2 ) >= 1 ; 
hTry = 1e-8 ;% For tol (relTol). an initial h of 8e-8 is used.
tn = 0 ; 
tnm1 = 0 ;
t = [tn ] ;
q = [qn' ] ;
T = [t ] ; 
Q = [q ] ;
isInitialStep = true ; 

tic
disp(['STARTED  thetaP = ', num2str(thetaP)])

while t(end) < tend
  if isInitialStep
    options = odeset('RelTol',tol,'AbsTol',tol*1e-2,'initialStep',hTry,'Events',@myEventFcn);%,'outputFcn',@myoutputfunc) ;
  else
    options = odeset('RelTol',tol,'AbsTol',tol*1e-2,'Events',@myEventFcn);%,'outputFcn',@myoutputfunc) ;
  end
  [t,q,te,qe,ie] = ode45(f,[tn tend],qn,options) ;% Stops upon event.

  if sqrt( (q(end-1,1))^2 + (q(end-1,3))^2 ) < 1 % Lastly no contact.
      isContact = 1 ;% Upon event there is contact. 
  elseif sqrt( (q(end-1,1))^2 + (q(end-1,3))^2 ) >= 1 % Lastly in contact
      isContact = 0 ;% Upon event no contact. 
  end
  qn = q(end,:) ;
  tn = t(end) ;
  Q = [Q; q] ; 
  T = [T; t] ;
    
end
Q = Q.' ; T = T.' ;% Take transpose to make it 4 long rows. 

TOC = toc ;
disp(['FINISHED thetaP = ', num2str(thetaP), ' in ',num2str(TOC),' sec.'])

% Zilli_individualplot(t,q,thetaP,[0.98 1], true) ; 


%% FUNCTIONS %%

function dq = f_sta(tn,qn)
    %fprintf("%6.2f of %4.2f\n",tn/tend*100,thetaP) % Track progress.
    
    epsH = 3.53e-1 ; 
    mH = 0.9 ;% 0.367 ; 
    beta = 1.32 ; 
    JpH = 0.143 ; 
    a_b = 1.40 ; 
    zeta = 0.01 ;%0.01 ; 
    g = 0 ; 
    c = 0.01 ;% Not needed as g=0. I couldn't find from the Zilli's paper.
    cS = c*a_b ;
    thetaPP = 0 ;
    
    q1 = qn(1) ;% phi_x_Hat
    q2 = qn(2) ;% phi_x_Hat_prime
    q3 = qn(3) ;% phi_y_Hat
    q4 = qn(4) ;% phi_y_Hat_prime
    
%   isContact = sqrt( q1^2 + q3^2 ) >= 1 ;% slows down as calcd each time.
    
    % Prepare thetaP and tnm1
    if ~isInitialStep % If not the initial step. 
      hnm1 = (tn-tnm1) ;% h for prev step used to calc thetaP for this step.
      theta = theta + hnm1*thetaP ;% tnm1: time atW prev step was calcd.
      tnm1 = tn ;%next step preparation. 
    else % If it is the initial step, update the angle with specified hTry. 
      theta = theta + (hTry/4)*thetaP ; %The specified first step size is always devided by four - see in hMat ??
      isInitialStep = false ;
    end
    % After the simulation is terminated with event, the tnm1 is kept for
    % the next increment. So the next increment is going to continue from
    % the "if ~isInitialStep" part. 

    dq1 = q2;
    
    dq2 = - JpH * thetaP * q4 ...
          - 2 * zeta * q2 ...
          - q1 ...
          + mH * epsH * ( thetaPP * cos(theta) - thetaP^2 * sin(theta) )...
          + mH * g / cS ...
          + isContact * beta * q1 * ( 1/sqrt( q1^2 + q3^2 ) - 1 ) ;
    
    dq3 = q4;
    
    dq4 = + JpH * thetaP * q2 ...
          - 2 * zeta * q4 ...
          - q3 ...
          + mH * epsH * ( thetaPP * sin(theta) + thetaP^2 * cos(theta) )...
          + isContact * beta * q3 * ( 1/sqrt( q1^2 + q3^2 ) - 1 ) ;

    dq = [dq1;dq2;dq3;dq4];
    
end
%%
function dq = f_rot(tn,qn)
    epsH = 3.53e-1 ;
    mH = 0.9 ;% 0.367 ; 
    beta = 1.32 ; 
    JpH = 0.143 ;
    zeta = 0.01 ;
    thetaPP = 0 ; 
    
    q1 = qn(1) ;% u : the x-axis of rotating frame. 
    q2 = qn(2) ;% u_prime
    q3 = qn(3) ;% v
    q4 = qn(4) ;% v_prime
    
%   isContact = sqrt( q1^2 + q3^2 ) >= 1 ;% slows down as calcd each time.
    
    dq1 = q2 ; 
    
    dq2 = - thetaP * (JpH-2) * q4 ...
          - 2 * zeta * q2 ...
          + ( thetaP^2 * (1-JpH)  -  1 ) * q1 ...
          + 2 * zeta * thetaP * q3 ... % this is - @2019 Shaw, but then shoots to inf.
          + mH*epsH*thetaP^2 ...
          + isContact * ( 1/sqrt(q1^2+q3^2) - 1 ) * beta * q1 ;
    
    dq3 = q4 ; 
    
    dq4 = + thetaP * (JpH-2) * q2 ...
          - 2 * zeta * q4 ...
          + ( thetaP^2 * (1-JpH)  -  1 ) * q3 ...
          - 2 * zeta * thetaP * q1 ... % this is + @2019 Shaw, but then shoots to inf. 
          - mH*epsH*thetaPP ...
          + isContact * ( 1/sqrt(q1^2+q3^2) - 1 ) * beta * q3 ;
      
    dq = [dq1;dq2;dq3;dq4];

    % In matrix form: 
%     dq = [0 1 0 0
%          thetaP^2*(1-JpH)-1+isContact*(1/sqrt(q1^2+q3^2)-1) , -2*zeta,  2*zeta*thetaP ,  thetaP*(2-JpH)
%          0 0 0 1 
%         -2*zeta*thetaP , -thetaP*(2-JpH) , thetaP^2*(1-JpH)-1+isContact*(1/sqrt(q1^2+q3^2)-1) , -2*zeta]*qn ; 
    
end


%% contactChange event

function [value,isterminal,direction] = myEventFcn(tn,qn)
    value = sqrt( qn(1)^2 + qn(3)^2 ) - 1 ;
%   isterminal = 0 ;% Not terminated as isContact is updated inside f.
    isterminal = 1 ;% Terminated upon a contactChange to update isContact. 
    direction = 0 ;
end
%% myoutputfunc 
% 
% function status = myoutputfunc(t,q,flag)
%     persistent hMat
%     switch flag
%         case 'init'
%             hMat = [hMat h];% Initially the h must be hDefault, otherwise it is the last h before contactChange
%         case ''
%             hMat = [hMat, t(end)-t(end(-1) ];
%         case 'done' % when it's done
%             assignin('base','hMat',hMat); % get the data to the workspace.
%     end
%     status = 0;
% end

end


