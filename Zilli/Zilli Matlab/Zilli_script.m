clc, clear,% close all

%% Solve for a range of thetaP (= theta_prime = omega_hat)
solvers = ["mySolver", "ode45"] ;
coords = ["sta" "rot"] ;
% 'sta' - stationary in {phixH; phixHP; phiyH; phiyHP}
% 'rot' - rotating X in {    u;     uP;    vH;     vP}
solver = solvers(2) ;
for coord = coords
startClock = tic ;
coord
tend = 1000 ;
tol = 1e-8 ; 
thetaP_range = 0.02:0.02:7 ;%6.02:0.02:6.42 % 0.02:0.02:4 % 4.02:0.02:7  % 
T = {};Q = {};
p = 1 ; remainderIndexCount = length(thetaP_range) ;
while  remainderIndexCount > 0
  indexRange = p:p+min(10,remainderIndexCount)-1 ; 
  disp(['NEW INDEX BATCH ' num2str(indexRange)])
  switch solver 
    case "mySolver"
      for i = indexRange
        %! DECORATION AROUND PARFOR: IT DOES NOT DISTRIBUTE THE LOAD AFTER 
        %! THE FIRST DISTRIBUTION. IT SWITCHES TO NORMAL FOR !!!
        %! THIS DECORATION HELPED DECREASE ONE FROM 950 TO 110 SEC. 
        thetaP = thetaP_range(i) ;
        qn = [0.5;0;0;0.5] ;% defined in phixH phixHP phiyH phiyHP coords.
        if coord == "rot" 
          qn = [1 0 0 0;0 1 -thetaP 0;0 0 1 0;thetaP 0 0 1]^(-1)*[0 0 1 0;0 0 0 1;-1 0 0 0; 0 -1 0 0]*qn ;
        end
        [T{i},Q{i},isContactMat{i},hMat{i}] = Zilli_func(thetaP,qn,tend,tol,coord);
      end    

    case "ode45"
      parfor i = indexRange
        thetaP = thetaP_range(i) ; 
        qn = [0.5;0;0;0.5] ;% defined in phixH phixHP phiyH phiyHP coords.
        if coord == "rot" 
          qn = [1 0 0 0;0 1 -thetaP 0;0 0 1 0;thetaP 0 0 1]^(-1)*[0 0 1 0;0 0 0 1;-1 0 0 0; 0 -1 0 0]*qn ;
        end
        [T{i},Q{i}] = Zilli_func_ode45(thetaP,qn,tend,tol,coord) ;
        hMat{i} = [];
        for j = 1:length(T{i})-1
          hMat{i}(j) = T{i}(j+1)-T{i}(j) ; 
        end
      end
  end

  remainderIndexCount = length(thetaP_range) - indexRange(end); 
  p = indexRange(end) + 1 ; 
end
finishClock = toc(startClock) 
disp(['COMPLETED in ', num2str( finishClock / 60 ), ' mins.'])

%% Visualise a trajectory: Individual checks. 
% for thetaP_check = 1.1%3.38 %[1.01 3.03 3.21 3.42 3.49 3.62 3.67 3.68 3.82 3.85] 
%   I = find( abs( thetaP_range-thetaP_check ) <= 0.001 );
%   Zilli_individualplot(T{I},Q{I},thetaP_check,[0.90 1], false) ; 
%   % TF: Animation
% end 

%% Visualise solutions: A broader view with  3D-Plot & Poincare section

%-OOO Prepare the figure
figure
ax1 = subplot(121) ; hold(ax1,'on') ; grid(ax1,'on')
ax2 = subplot(122) ; hold(ax2,'on') ; grid(ax2,'on')

if coord == "sta"
  xlabel(ax1,"q1 (Phi_x or -Y^H^a^t)")
  ylabel(ax1,"q2 (Phi_y or X^H^a^t)")
  zlabel(ax1,"Omega^H^a^t (or thetaP )")
  title (ax1,['Trajectory stack for ' solver ' ' coord]) 
elseif coord == "rot"
  xlabel(ax1,"q1 (u)")
  ylabel(ax1,"q2 (v)")
  zlabel(ax1,"Omega^H^a^t (or thetaP )")
  title (ax1,['Trajectory stack for ' solver ' ' coord]) 
end
xlim(ax1,[-3 3]), ylim(ax1,[-3 3])

xlabel(ax2,"Omega^H^a^t (or thetaP)")
ylabel(ax2,"Normalised Ampitude")
title (ax2,"Amplitude")
ylim(ax2,[0 5])
%-OOO END

for i=1:length(thetaP_range)
  from = round(0.98*length(Q{i}),0) ;  
  to = length(Q{i}) ;
  thetaP = thetaP_range(i) ; 

  %-111 3D Visualise 
  plot3( ax1, ...
         Q{i}(1,from:end) , ...
         Q{i}(3,from:end) , ...
         thetaP*ones(1,length(from:to)) ,...
         "k." )
  %-111 END

  %-222 Amplitude plot through from:to indexes
  rH = sqrt( Q{i}(1,:).^2 + Q{i}(3,:).^2 ) ;
  plot( ax2, thetaP*ones(1,length(from:to)) , rH(from:to) , "k.")
  %-222 END



  %-AAA Stroboscopic Poincare section 
  %-UUU Pack all period-close indexes into IPeriod.
  TPeriod = 2*pi/thetaP ;
  k = 1 ; r = 1 ; IPeriod = []; 
  while k*TPeriod <= T{i}(end)
    if k*TPeriod > T{i}(from) && k*TPeriod <= T{i}(to)
      [~,IPeriod(r)] = min( abs( T{i} - k*TPeriod ) ) ;
      r = r+1 ;
    end
    k = k+1 ;
  end
  %-UUU END
  plot3( ax1, ...
         Q{i}(1,IPeriod) ,...
         Q{i}(3,IPeriod) , ...
         thetaP*ones(1,length(IPeriod)),...
         "bo" )

  plot( ax2 , thetaP*ones(1,length(IPeriod)) , rH(IPeriod) , "bo" )
  %-AAA END

end



end