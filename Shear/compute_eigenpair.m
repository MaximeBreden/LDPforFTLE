function X_grid = compute_eigenpair(K,N,problem,show)

% Approximate computation of the first eigenvalue of L_p, first using
% Matlab's eigensolver, and then refined using Newton's method.

if nargin < 4
    show = false;
end

[alpha,b,sigma,~,~,~,p_grid] = obtain_parameters(problem);

[m0,m1] = coeffs_operator(alpha,b,N);
Der = spdiags( transpose(1i*(-N:N)), 0, 2*N+1, 2*N+1 );
% The Toeplitz matrix associated to m0 (discrete convolution product)
M0 = spdiags( m0, [1,0,-1], 2*N+1, 2*N+1 );
% The Toeplitz matrix associated to m1 (discrete convolution product)
M1 = spdiags( m1, [1,0,-1], 2*N+1, 2*N+1 );

X_grid = zeros(2*N+2,K+1);
tol = 1e-12;
it_max = 100;
for j = 1:K+1
    p = p_grid(j);
    L = 2*sigma^2*Der^2 + M1*Der + p/2*M0;% The differential operator
    [eigenvect, eigenval] = eig(full(L));
    eigenval = diag(eigenval);
    [eigenval, ind] = sort(eigenval,'ComparisonMethod','real');
    lambda = eigenval(end);
    u = eigenvect(:,ind(end));
    X = symmetrize([u;lambda]);
    problem.normalization = u;

    % Newton
    problem.p = p;
    X_grid(:,j) = Newton_eigenproblem(X,problem,tol,it_max,show);
end