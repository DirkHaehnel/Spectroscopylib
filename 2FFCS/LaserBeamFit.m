function [err, c, x] = LaserBeamFit(p,lam,z,y,flag,pic,weight)

% Fitting the Gaussian radius of a Gauss beam
% (c) Joerg Enderlein, http://www.joerg-enderlein.de (2008)

if length(p)==1
    p = [p 0];
end
if length(p)>2
    x = p(1)*sqrt(1+(lam*(z-p(2))/pi/p(3)^2).^2);
else
    x = p(1)*sqrt(1+(lam*(z-p(2))/pi/p(1)^2).^2);
end
if nargin>3 
    if nargin<5 || isempty(flag)
        c = x(:)\y(:);
    else
        c = 1;
    end
    x = x*c;
    if nargin>5 && ~isempty(pic)
        plot(z,y,'o',z,x); drawnow
    end
    if nargin>6 && ~isempty(weight)
        err = sum((y(:)-x(:)).^2./x(:).*weight(:));
    else
        err = sum((y(:)-x(:)).^2./x(:));
    end
else
    err = x;
end
