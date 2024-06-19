function y = f2i(x,test)

%Routine for converting floats to intvals when needed (typically an array 
%of zeros), without generating an error when Intlab is not used.

if exist('intval','file') && isintval(test)
    y = intval(x);
else
    y = x;
end