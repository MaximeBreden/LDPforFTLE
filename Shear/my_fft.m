function Y = my_fft(X,n,dim)

if exist('intval','file') && isintval(X(1))
    if nargin == 3
        sz = size(X);
        if dim==3 && length(sz)==3
            Y = intval(zeros(sz));
            for i = 1:sz(1)
                Y(i,:,:) = transpose(verifyfft(transpose(squeeze(X(i,:,:)))));
            end
        elseif dim==2 && length(sz)==2
            Y = transpose(verifyfft(transpose(X)));           
        else
            error("this case is not implemented yet")
        end
    else
        Y = verifyfft(X);
    end     
else
    if nargin == 3
        Y = fft(X,n,dim);
    elseif nargin == 2
        Y = fft(X,n);
    else
        Y = fft(X);
    end
end