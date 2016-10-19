function [n i j]=mintersect(y, ylag1)

pos0=1; pos1=1;

i=zeros(length(y),1); j=i; n=i;
c=0;
while (pos0<numel(y)+1) && (pos1<numel(ylag1)+1)
    
    if y(pos0)>ylag1(pos1)
        pos1=pos1+1;
    elseif y(pos0)<ylag1(pos1)
        pos0=pos0+1;
    else
        c=c+1;
        n(c)=ylag1(pos1);
        i(c)=pos0;
        j(c)=pos1;
        pos0=pos0+1;
        pos1=pos1+1;
    end
    
end

i=i(1:c); j=j(1:c); n=n(1:c);