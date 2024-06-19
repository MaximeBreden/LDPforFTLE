function y = icumsum(x,dim,dir)

% A version of the Matlab built-in cumsum function which works for intvals

if isempty(x)
    y = x;
else
    if nargin < 3
        dir = 'forward';
    end

    if not(exist('intval','file')) || not(isintval(x(1)))
            y = cumsum(x,dim,dir);
    else
        rndold = getround;
        setround(-1)
        yinf = cumsum(x.inf,dim,dir);
        setround(1)
        ysup = cumsum(x.sup,dim,dir);
        y = infsup(yinf,ysup);
        setround(rndold)    % set rounding to previous value
    end
end