function [m0,m1] = coeffs_operator(alpha,b,N)

% b*(1+cos(phi)) in the exponential basis
m1 = ones(2*N+1,1) * [b/2, b, b/2];

% b*sin(phi) - 2*alpha in the exponential basis
m0 = ones(2*N+1,1) * [1i*b/2, -2*alpha, -1i*b/2];
