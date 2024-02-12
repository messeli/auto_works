function [hbalance, omega_r, Np0 ]=NF_AFT_target( Uin, Lambdau, np, H0, Sk, Sl , zin  , zeqn )
%|:aFunc&itsInputs r gvn as input to odr func 

%add the 'zero coeefficient' and then convert Uin to complex form
Uin2 = [ Uin( 1:(zin-1) ); 0; Uin(zin:end) ];
U = Uin2(1:2:end) + 1j*Uin2(2:2:end);
NResTerms=length(U); % this is just the number of terms considered resonant -typically no more than 2

%get some parameters from H0
[~,Nft] = size(H0);

%add U to H0 to get P %|eq 16 of Shaw2019 NF article
P0=H0;   % remember this only apporixmates P - it P assuming H=H0
for ii=1:NResTerms
    P0(Sk(ii),mod(Sl(ii),Nft)+1) = U(ii); %note - H0 at this term should always be 0 
    %|:P=U+H, but H for resonant terms is zero
    %|:Also note that the P=U+H is origY P*uStar=U*uStar+H*uStar in d analytic NF theory.
    %|...Bt here instead f uStar, we hav u=U*t (eq12 of Shaw2019); t(l) = exp(1j*l*omr*time) 
    %|...The difficulty is due to the omeg_r being unknown. 
    %|...Also, d method doesnt refer t d time ever bt t d U (similarly, t H0 P0) 
    %|:H0 and P0 will have a dim f  dof*n_fourier (ie here designated by N*nf)
    %|...U'd hav d sam dim (eq12) but we already state dt drr only 2 modes, so U is reducD t 2*1,
    %|...whr 1st el is the complex amp f Sl(1)*wr freq and the second is the 2nd simY.
end

%time series of p
p = ifft(P0,Nft,2)*Nft; %2 makes it go over rows not columns %|eq 17 of Shaw2019 NF article
%|:Sizes [p]=Ndof*N_timeData=2*1024  and  [np]=Ndof*N_fourierTransform=2*1024

%evaluate np and Np %|eq 18*19 of Shaw2019 NF article
Np0  = fft(np(p),Nft,2)/Nft;

%reduce Np to Nu  (just the resonant mode/frequency combinations) 
Nu=zeros(NResTerms,1);
for ii =1:NResTerms
    Nu(ii) = Np0(Sk(ii),mod(Sl(ii),Nft)+1);
end
%|:eq 22 (with 21 skipped), of Shaw2019 NF article.
%|:Indeed dss similar to the lines 13-15 here, bt as each component cn b isolatD into a diff eq 
%|...and becoz we only need d resonant parts of the whole Nu (dim: 2*n_ft), 
%|...we just get the res part of Nu, whc s already Np due to eq 22 of the 2019Shaw article.

% find omega_r using one row of eqns
% calculate the 'reduced' Psi matrix need for the HB calculation
Psi1 = diag(1j * Sl ); 
PsiU1 = Psi1*U;
LambdauU_plus_Nu = Lambdau*U +Nu;
A1 = ExpandComplex2RealVec(PsiU1);
A2 =  ExpandComplex2RealVec(LambdauU_plus_Nu); %losing all ability to concoct sensible array names
%use row zeqn of this to calculate omega %|zeqn=2,so like! article says,dss Im f first! line f eq23.
omega_r = A2(zeqn) / A1(zeqn);
% %use least squares for omega_r
% omega_r = sum(A2.*A1)/sum(A1.*A1);
Psi = Psi1*omega_r;

%evaluate the reduced harmonic balance and break into real and imaginary
%components
hb=ExpandComplex2RealVec(      (Psi-Lambdau)*U - Nu            );

%return result, ignoring the solved equation
hbalance= hb(  [ 1:(zeqn-1)  (zeqn+1):end ] );%|Dss minimised in the nonliner solver. 


