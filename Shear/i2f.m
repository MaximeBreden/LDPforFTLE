function y = i2f(x,opt)

%Routine for dealing with data without having to worry whether the input is
%a usual float or an intval. If the input is an intval, the default output
%is the supremum

if exist('intval','file') && isintval(x(1))
    if nargin == 2 
        if strcmp(opt,'inf')
            y = inf(x);
        elseif strcmp(opt,'mid')
            y = mid(x);
        elseif strcmp(opt,'sup')
            y = sup(x);
        else
            error("Invalid input for 'opt'")
        end   
    else
        y = sup(x);
    end
else
    y = x;
end