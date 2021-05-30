% function zilli_noncons_normalforms_rj2

% Analyse the Zilli system with nonconservative effects using Normal forms and
% an AFT HB method, plus energy transfer equation.
% Alex S
clear
close all
clc

%Define system matrices - Stationary system
mass=1;  % nondimensionalised -remember this is an effective mass, = Id+a^2 mdisc
Jp=0.14;
zeta=0.01;
k=1;
[M,G,C,K] = zillStationaryMatrices( mass, Jp , zeta, k );%mass, Jp, zeta, k
m=3.67e-1;   % out of balance mass
epsln=  0.353/2;  % out of balance distance
f=m*epsln;

J=[ 0 -1;
    1  0 ];  %skew symmetric matrix

% contact stiffness and clearance
beta =  1.32e1;
rc=1;%This cannot change, we assume a nondimensionalised system above
nq0 = @(q) ZilliContactStiffness(q,rc,beta);

% contact damping = assume a piecewise damping constant that is equivalent
% to a coefficient of restitution d
d = 1;
zeta_s =sqrt(    ( log(d) /pi)^2   /   (  1 + ( log(d) /pi)^2  )     );
cs = zeta_s * 2 * sqrt(beta);
nqzeta = @(q,qdot) ZilliContactDamping(q,qdot,rc,cs); %note- this doesn't include nonconservative forcing, only dampinng

%overall nonlin function
nq = @(q,qdot) ZilliContactStiffness(q,rc,beta)+ZilliContactDamping(q,qdot,rc,cs);

%model case structure used in time sim
zmdl = CreateZilliMdl( mass, Jp, zeta, k , m, epsln, beta, rc ,cs);

%% linear system critical speed and syncronous response- stationary and rotating coords
%campbell diag
Omvals=6.0:0.025:6.4;%5.8:0.02:10;%3.8%:0.05:5.0;%3.37:0.0005:3.39;%2:0.01:2.5;%2.2:0.01:2.5;% 3.2:0.01:4.5;%3.7:0.01:3.9;   %0.2:0.2:6.4;%
nOmvals=length(Omvals);
[ omfw , ombw ] = zillWhirlSpeeds_Stationary( M,G,C,K,Omvals );
[ omfw_ , ombw_ ] = zillWhirlSpeeds_Rotating( M,G,C,K,Omvals );
gmfigure
subplot(1,2,2)
mrksize=3;
plot( Omvals ,   omfw_ , 'k-+' ,'markersize',mrksize)
hold on
plot( Omvals ,   ombw_ , 'k-^' ,'markersize',mrksize)
plot( Omvals ,   2*omfw_ , 'k--' )
plot( Omvals ,   3*omfw_ , '--' )
plot( Omvals ,   4*omfw_ , '--' )
plot( Omvals ,   2*ombw_ , 'k:' )
% plot( Omvals ,   3*ombw_ , ':' )
grid minor
%grid on
fntsz=12;
xlabel('$\hat{\Omega}$' , 'interpreter' , 'latex', 'fontsize',fntsz  )
ylabel('Rotating system angular velocity' , 'interpreter' , 'latex', 'fontsize',fntsz )
xlim([ 0 Omvals(end)]);
title('(b)',  'interpreter' , 'latex' , 'fontsize',fntsz)
%interpolate crossing points
Om21 = GetCrossings( Omvals , ombw_ - 2*  omfw_  )
om21 = interp1(   Omvals , ombw_ , Om21) 
plot(Om21,om21,'ro','markersize',15)
Om32 = GetCrossings( Omvals , 2*ombw_ - 3*  omfw_  )
om32 = interp1(   Omvals , 2* ombw_ , Om32) 
plot(Om32,om32,'ro','markersize',15)
legend({'$\tilde{\omega}_{n2}$', '$\tilde{\omega}_{n1}$', ...
    '$2\tilde{\omega}_{n2}$',  '$3\tilde{\omega}_{n2}$' ,  '$4\tilde{\omega}_{n2}$' , ...
    '$2\tilde{\omega}_{n1}$'  } , 'interpreter' , 'latex' , 'fontsize',fntsz ...
    ,'location','southwest','autoupdate','off')
subplot(1,2,1) %'traditional' campbell
plot( Omvals ,   omfw , 'k-+' ,'markersize',mrksize)
hold on
plot( Omvals ,  abs( ombw), 'k-^','markersize',mrksize )
plot( Omvals ,   Omvals, ':' )
xlim([ 0 Omvals(end)]);
% ylim([0 3]);
%grid on
legend({'$|{\omega}_{n2}|$', '$|{\omega}_{n1}|$', 'Drive speed'  } , 'interpreter' , 'latex' , 'fontsize',fntsz ...
    ,'location','northwest','autoupdate','off')
xlabel('$\hat{\Omega}$' , 'interpreter' , 'latex', 'fontsize',fntsz  )
ylabel('Stationary system whirl speed' , 'interpreter' , 'latex', 'fontsize',fntsz )
title('(a)',  'interpreter' , 'latex' , 'fontsize',fntsz)

%% oob response
X_=zeros(2,length(Omvals));
for ii=1:nOmvals
    Omega=Omvals(ii);
    
    %rotating sys matrices - note underscores indicate a rotating system
    %property
    [G_,K_,Kc_] = zillRotatingMatrices( M,G,C,K,Omega );
    
    %calc response to out of balance
    %rotating
    b_ = Omega^2*m*epsln*[1;0];
    X_(:,ii) = ( K_+Kc_) \ b_  ; %note this is a real calculation in rotating, which is  nice
    
end
omcrit=GetCrossings(Omvals,omfw-Omvals)
gmfigure
% plot(Omvals, abs(X(1,:))  );
hold on
plot(Omvals, sqrt(X_(1,:).^2 + X_(2,:).^2) ,'-');
%grid on

%% Solve periodic (in rotating frame) responses - forced/damped
Ndof=2;
%AFT parameters
Nft = 2^10; %number of fourier terms/time steps in time evaluation
nguess=4;
Us = zeros(2,nguess*nOmvals);
solvecodes=zeros(1,nguess*nOmvals);
omrs=zeros(1,nguess*nOmvals);
Omvals=[Omvals Omvals Omvals Omvals];
% pick resonant terms - hard coded for now
Sk = [ 1 2 ];  %the kth modes
Sl = [ -3 -2 ];      % the lth harmonics of omega_r
NResTerms=length(Sk);  %number of resonant terms

forb=gmfigure;
cptns={'(a) Upper branch solution','(b) Lower branch solution'};

for ii =1:(nguess*nOmvals)
    Omega=Omvals(ii);
    
    % obtain rotating frame matrices
    [G_,K_,Kc_] = zillRotatingMatrices( M,G,C,K,Omega );
    b_ = [ f*Omega^2;
        0  ];
    
    % (conservative linear) state space matrix and nonlin functions
    A_ = [   zeros(2)    ,  eye(2)  ;
        -M\K_        ,  -Omega*(M\G_)];
    % nonlinear function
    ny = @(y) [   zeros(Ndof,size(y,2));  M\nq0(y(1:Ndof,:)) ]  ...
        + [ zeros(Ndof,2*Ndof) ; -M\Kc_,    -M\C    ] * y ...
        + [ zeros(Ndof,size(y,2));  M\nqzeta(y([1,2],:),y([3,4],:)) ] ...
        + [ zeros(Ndof,size(y,2));  -M\b_*ones(1,size(y,2)) ];
    
    %modal transformation matrices
    [Phi_2n,Lambda_2n]=eig(A_);
    %pick out eigenvalues of form [ a, -ja , ..]'
    flt = imag(   Phi_2n(2,:) ./ ( Phi_2n(1,:) )   ) < 0;
    Phi = Phi_2n(:,flt);
    Lambda = Lambda_2n( flt , flt) ;
    %sort eigenvalues/vectors
    [~,ix]=sort(abs(imag(diag(Lambda))+Omega));  %sorted by magnitude of stationary systen whirl speed
    Phi = Phi(:,ix);
    Lambda = Lambda(ix,ix);
    %reassemble full size matrices
    Phi_2n = [ Phi conj(Phi) ];
    Lambda_2n = [ Lambda    ,   zeros(2);
        zeros(2)  ,   conj(Lambda) ];
    
    %create (modally transformed) nonlinear functions
    Phinv_2n = inv(Phi_2n);
    np = @(p) Phinv_2n(1:Ndof,:) * ny( 2*real(Phi*p) );
    
    
    %create NF transformed/reduced system variables
    % note - I am sure some precalculation could speed up this bit but this
    % method will hopefully get something working!
    Lambda_u = Lambda( Sk , Sk );
    uflt = false(NResTerms,Ndof); % a matrix to pull out just the resonant modes from modally transformed funcs
    for jj=1:NResTerms
        uflt(jj,Sk(jj))=true;
    end
    
    % set up HO
    Y0=-[X_(:,mod(ii-1,nOmvals)+1); 0;0 ]; %the linear out of balance response in rotating coords
    P0=Phinv_2n(1:Ndof,:)*Y0;
    H0 = zeros(Ndof,Nft);
    H0(:,1) = P0;  %note - we use the Matlab convention for FFT harmonics i.e.
    % 1 = DC component ,  2:(nft/2) - positive freq components,
    % (nft/2+1):nft - negative freq components with -1 last.
    % This is followed in calculations (but fixed for graphs etc).
    % shorthand: matlabharmonic = mod(trueharmonic,nft)+1;
    for jj=1:NResTerms %zero the resonant components -they end up in U
        H0(Sk(jj),mod(Sl(jj),Nft)+1)=0;
    end
    
    fctr=2.5;
    switch floor((ii-1)/nOmvals)
        case 0
            Uin0 = fctr*[ 1 -1 0  ]';
        case 1
            Uin0 = fctr*[ 1 0 -1  ]';
        case 2
            Uin0 = fctr*[ 1 1 0   ]';
        case 3
            Uin0 = fctr*[ 1 0 1   ]';
        otherwise
            Uin0 = [ 0 0 0   ]'; %shouldn't get here!
    end
    for kk=1:3
        opt= optimoptions('fsolve','Display','none','useparallel',true,'FunctionTolerance',1e-8);
        NFHB = @(Uin) NF_AFT_target( Uin, Lambda_u, np , H0, Sk, Sl , 2 ,2 );
        [U_,~,solvecodes(ii)]=fsolve( NFHB , Uin0 ,opt);
        %     solvecodes(ii)=sc;
        U = [ U_(1);
            U_(2)+1j*U_(3); ];
        Us(:,ii)=U;
        [~,omrs(ii),Np0]=NFHB(U_);
        
        %evaluate H (note H0 is just a guess at H or a previous iteration)
        Psidiag  = [ 0 1:(Nft/2-1) (-Nft/2:-1) ]*1j*omrs(ii)  ;  %i.e. a row matrix of the diagonal of Psi
        Beta = zeros(size(Np0));
        for jj=1:Ndof
            Beta(jj,:) = Psidiag-Lambda(jj,jj);
        end
        H = Np0 ./ Beta; %note this will include some large values at the resonant terms - but these are about to be overwritten
        
        H0=H;  %iterate H to get a better sol
        Uin0=U_;
    end
    % evaluate P and then p
    P=H;
    for jj=1:NResTerms
        P(Sk(jj),mod(Sl(jj),Nft)+1)=U(jj); %copy the resonant terms into P
    end
    p = ifft(P,Nft,2)*Nft;
    
    % translate back to spatial coordinates
    y = 2*real(Phi *p);
    x = y(1:Ndof,:);
    
    %plot orbit - if it is interesting!
    if abs(U(2))> 0.001 && abs(U(1))> 0.001 && solvecodes(ii)>0
        
        
       fg=gmfigure;
        subplot(1,2,1)
        ho=plot(x(1,:),x(2,:),'--','linewidth',2);
        title(['Orbit $\Omega = ' sprintf('%f',Omega) '$'] , 'interpreter' , 'latex' , 'fontsize',fntsz);
      %  grid on
      xlim(1.1*[ -1 1])
      ylim(1.1*[ -1 1])
        axis square
        hold on
        rectangle('position',2*[-0.5 -0.5 1 1],'curvature',[1 1])
        ax=subplot(1,2,2);
        plot(-Nft/2:(Nft/2-1),abs(P(:,[(Nft/2+1:end)  1:(Nft/2)])'))
        %grid on
        xlabel('Harmonic $\ell$',  'interpreter' , 'latex' , 'fontsize',fntsz)
        ylabel('Complex modal amplitude',  'interpreter' , 'latex' , 'fontsize',fntsz);
        xlim([-Nft/2 Nft/2])
        ax.XTick=[-Nft/2:floor(Nft/16):Nft/2];
        legend({'$|P_{1,\ell}|$','$|P_{2,\ell}|$'}, 'location','Northeast',  'interpreter' , 'latex' , 'fontsize',fntsz,'autoupdate','off')
        title(sprintf('(b) $\\omega_r = %.3f$',omrs(ii)),  'interpreter' , 'latex' , 'fontsize',fntsz);
        xlim([-12 12]);
        %time simulation test!
        [~,resp ] = zillicontactsim(zmdl,Omega,[0 30*2*pi/omrs(ii)],y(:,1));
        subplot(1,2,1)
        plot(resp(1,:),resp(2,:),'-','linewidth',0.5);
        legend({'Analytical','Numerical'}, ...
            'location','southoutside',  'interpreter' , 'latex' , ...
            'fontsize',10,'orientation','horizontal','autoupdate','off')
        %         tflt=tt<2*pi/omrs(ii);
        %            plot(resp(1,tflt),resp(2,tflt),'--','linewidth',2);
        uistack(ho,'top');
        plot(x(1,1),x(2,1),'k*','markersize',10);
        fprintf('ii=%i, Omega=%f, omega_r=%f , U1=%f  , U2=%f+j%f , fig %i \n' , ii, Omega, omrs(ii) , U_(1), U_(2), U_(3) , fg.Number);
        
        
    end
    
    
    
end
flt=solvecodes>0;
gmfigure
plot(Omvals(flt),abs(Us(1,flt)),'+')
hold on
plot(Omvals(flt),abs(Us(2,flt)),'o')
legend({sprintf('$\\left| U_{%i,%i} \\right|$',Sk(1),Sl(1)), ...
    sprintf('$\\left| U_{%i,%i} \\right|$',Sk(2),Sl(2))}, ...
    'interpreter' , 'latex' , 'fontsize',fntsz ,  'location', 'northwest' ...
    ,'autoupdate','off')
xlabel('$\hat{\Omega}$','interpreter' , 'latex' , 'fontsize',fntsz)
ylabel('Complex modal amplitude','interpreter' , 'latex' , 'fontsize',fntsz)
%grid on
xlim([Omvals(1) Omvals(end)])
[solvecodes; Omvals]
gmfigure
plot(real(Us(1,flt)),imag(Us(1,flt)),'+')
hold on
plot(real(Us(2,flt)),imag(Us(2,flt)),'o')
%grid on
gmfigure
plot(Omvals(flt),omrs(flt),'+-')
xlabel('$\hat{\Omega}$','interpreter' , 'latex' , 'fontsize',fntsz)
ylabel('$\omega_r$','interpreter' , 'latex' , 'fontsize',fntsz)
% grid on
gmfigure
plot(omrs(flt),abs(Us(1,flt)),'+')
hold on
plot(omrs(flt),abs(Us(2,flt)),'o')
legend({sprintf('$\\left| U_{%i,%i} \\right|$',Sk(1),Sl(1)), ...
    sprintf('$\\left| U_{%i,%i} \\right|$',Sk(2),Sl(2))}, ...
    'interpreter' , 'latex' , 'fontsize',fntsz ,  'location', 'northwest' ...
    ,'autoupdate','off')
xlabel('$\omega_r$','interpreter' , 'latex' , 'fontsize',fntsz)
ylabel('Complex modal amplitude','interpreter' , 'latex' , 'fontsize',fntsz)
% grid on

% flt=solvecodes==1;
% if any(flt)
%     plot(Omvals(flt),0,'*')
% end
% flt=solvecodes<1;
% if any(flt)
%     plot(Omvals(flt),0,'o')
% end
% flt=solvecodes==2;
% if any(flt)
% plot(Omvals(flt),0,'x')
% end
% flt=solvecodes==3;
% if any(flt)
% plot(Omvals(flt),0,'+')
% end
% CreateLatexReport('nf_aft')