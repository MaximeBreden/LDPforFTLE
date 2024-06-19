function F = F_eigenproblem(X,problem,ext,type)

% The map F(u,lambda) = (L*u - lambda*u, normalization condition) 

% The optional parameter ext can be used to suppress truncation errors in
% Fourier.

sz = size(X);
N = (sz(1)-2)/2;
K = sz(2)-1;
if K > 0
    if isfield(problem,'data_type')
        type = problem.data_type;
    else
        type = 'grid';
        warning("By default, we assume x to be of grid type")
    end
end
u = X(1:2*N+1,:);
lambda = X(2*N+2,:);
u_normalization = problem.normalization;
[alpha,b,sigma,p] = obtain_parameters(problem);

% Parts not depending on k
Der = transpose(1i*(-N:N));
over_12 = b/2*Der; 
dia_12 = 2*sigma^2*Der.^2 + b*Der;
under_12 = b/2*Der;
over_0 = 1/2*(1i*b)/2;
dia_0 = 1/2*(-2*alpha);
under_0 = 1/2*(-1i*b)/2;
if nargin >= 3 && ext
    L_12 = spdiags([over_12,dia_12,under_12], [0,-1,-2], 2*N+3, 2*N+1 );
    L_0 = spdiags(ones(2*N+1,1)*[over_0,dia_0,under_0], [0,-1,-2], 2*N+3, 2*N+1 );
    I = spdiags(ones(2*N+1,1), -1, 2*N+3, 2*N+1 );
else
    L_12 = spdiags([over_12,dia_12,under_12], [1,0,-1], 2*N+1, 2*N+1 );
    L_0 = spdiags(ones(2*N+1,1)*[over_0,dia_0,under_0], [1,0,-1], 2*N+1, 2*N+1 );
    I = spdiags(ones(2*N+1,1), 0, 2*N+1, 2*N+1 );
end

% Parts depending on p
if K>0 && strcmp(type,'coeffs')
    Mp = spdiags( ones(K+1,1)*[p(2),p(1),p(2)], [1,0,-1], K+1, K+1 );
    Mp(1,2) = 2*p(2);
    Mlambda = convomat_Cheb(transpose(lambda));
    pu = u * transpose(Mp);
    lambdau = u * transpose(Mlambda);
else
    pu = u .* repmat(p,[2*N+1,1]);
    lambdau = u .* repmat(lambda,[2*N+1,1]);
end

% Normalization condition
if u == u_normalization
    F_normalization = zeros(1,K+1);
else
    if K>0 && strcmp(type,'coeffs')
        error("Not implemented yet")
    else
        F_normalization = sum(conj(u_normalization).*u,1) - ones(1,K+1);
    end
end
F = [L_12*u + L_0*pu - I*lambdau; F_normalization];





