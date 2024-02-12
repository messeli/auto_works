% like figure but gives nice classical proportions..
function fho=gmfigure(h,rows)
if nargin==0
fh=figure;
else
    fh=figure(h);
end
if nargin<2
rows=1;
else
    %
end


op=fh.OuterPosition;%[left bottom width height]
op(3) =op(4)*(1+sqrt(5))/2;
op(2)=op(2)-op(4)*(rows-1);
op(4)=op(4)*rows;
fh.OuterPosition=op;
if nargout==1
    fho=fh;
end