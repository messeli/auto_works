function crossing = GetCrossings(xvals,yvals)
my_sign = yvals(1);
for i=2:length(xvals)
   if my_sign~=sign(yvals(i))
      I = [i-1,i] ;
      break
   end
end
crossing = xvals(I(1)) - abs( yvals(I(1))/yvals(I(2)) ) * ( xvals(I(1)) - xvals(I(2)) ) ;

end

