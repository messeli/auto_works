function r = ExpandComplex2RealVec( z )

r= zeros(2*length(z),1);

r(2:2:end)= imag(z);
r(1:2:(end-1))=real(z);