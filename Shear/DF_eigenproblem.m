function DF = DF_eigenproblem(X,problem,ext)

% The Frechet derivative of the map 
% F(u,lambda) = (L*u - lambda*u, normalization condition) 

% The optional parameter ext can be used to suppress truncation errors in
% Fourier

sz = size(X);
N = (sz(1)-2)/2;
K = sz(2)-1;
if K > 0
    if not(isfield(problem,'data_type')) || ( isfield(problem,'data_type') && not(strcmp(problem.data_type,'grid')) )
        warning("x is assumed to be of grid type")
    end
end

u = X(1:2*N+1,:);
lambda = X(2*N+2,:);
u_normalization = problem.normalization;
[alpha,b,sigma,p] = obtain_parameters(problem);

if nargin >= 3 && ext
    Der = transpose(1i*(-(N+1):N+1));
else
    Der = transpose(1i*(-N:N));
end
over_12 = b/2*Der; 
dia_12 = 2*sigma^2*Der.^2 + b*Der;
under_12 = b/2*Der;
over_0 = 1/2*(1i*b)/2;
dia_0 = 1/2*(-2*alpha);
under_0 = 1/2*(-1i*b)/2;
if nargin >= 3 && ext
    L_12 = spdiags([over_12,dia_12,under_12], [0,-1,-2], 2*N+5, 2*N+3 );
    L_0 = spdiags(ones(2*N+3,1)*[over_0,dia_0,under_0], [0,-1,-2], 2*N+5, 2*N+3 );
    I_N = spdiags(ones(2*N+3,1), -1, 2*N+5, 2*N+3 );
    u = [zeros(2,K+1);u;zeros(2,K+1)];
    u_normalization = [zeros(1,K+1);u_normalization;zeros(1,K+1)];
else
    L_12 = spdiags([over_12,dia_12,under_12], [1,0,-1], 2*N+1, 2*N+1 );
    L_0 = spdiags(ones(2*N+1,1)*[over_0,dia_0,under_0], [1,0,-1], 2*N+1, 2*N+1 );
    I_N = spdiags(ones(2*N+1,1), 0, 2*N+1, 2*N+1 );
end
p = reshape(p, [1,1,K+1]);
lambda = reshape(lambda, [1,1,K+1]);
if K == 0
    Lml = L_12 + p*L_0 - lambda*I_N;
else %also works for K=0 but only after converting sparse matrices to full ones
    Lml = superkron(L_12,ones(1,1,K+1)) + superkron(L_0,p) - superkron(I_N,lambda);
end

DF = [Lml, -reshape(u,[size(u,1),1,K+1]);
      reshape(conj(u_normalization),[1,size(u_normalization,1),K+1]), zeros(1,1,K+1)];
