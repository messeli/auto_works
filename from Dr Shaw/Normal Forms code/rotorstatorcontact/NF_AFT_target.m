function [hbalance, omega_r, Np0 ]=NF_AFT_target( Uin, Lambdau, np, H0, Sk, Sl , zin  , zeqn )

%add the 'zero coeefficient' and then convert Uin to complex form
Uin2 = [ Uin( 1:(zin-1) ); 0; Uin(zin:end) ];
U = Uin2(1:2:end) + 1j*Uin2(2:2:end);
NResTerms=length(U); % this is just the number of terms considered resonant -typically no more than 2

%get some parameters from H0
[~,Nft] = size(H0);

%add U to H0 to get P
P0=H0;   % remember this only apporixmates P - it P assuming H=H0
for ii=1:NResTerms
    P0(Sk(ii),mod(Sl(ii),Nft)+1) = U(ii); %note - H0 at this term should always be 0 %|P=U+H, but H for resonant terms is zero
end

%time series of p
p = ifft(P0,Nft,2)*Nft; %2 makes it go over rows not columns

%evaluate np and Np
Np0  = fft(np(p),Nft,2)/Nft;

%reduce Np to Nu  (just the resonant mode/frequency combinations)
Nu=zeros(NResTerms,1);
for ii =1:NResTerms
    Nu(ii) = Np0(Sk(ii),mod(Sl(ii),Nft)+1);
end

% find omega_r using one row of eqns
% calculate the 'reduced' Psi matrix need for the HB calculation
Psi1 = diag(1j * Sl ); 
PsiU1 = Psi1*U;
LambdauU_plus_Nu = Lambdau*U +Nu;
A1 = ExpandComplex2RealVec(PsiU1);
A2 =  ExpandComplex2RealVec(LambdauU_plus_Nu); %losing all ability to concoct sensible array names
%use row zeqn of this to calculate omega
omega_r = A2(zeqn) / A1(zeqn);
% %use least squares for omega_r
% omega_r = sum(A2.*A1)/sum(A1.*A1);
Psi = Psi1*omega_r;

%evaluate the reduced harmonic balance and break into real and imaginary
%components
hb=ExpandComplex2RealVec(      (Psi-Lambdau)*U - Nu            );

%return result, ignoring the solved equation
hbalance= hb(  [ 1:(zeqn-1)  (zeqn+1):end ] );


