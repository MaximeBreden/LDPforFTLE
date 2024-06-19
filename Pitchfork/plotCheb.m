function u = plotCheb(u,a,b,nb)
% Plots the function represented by the Chebyshev coefficients in u, on
% [a,b] and discretized with step size nb.
% Normalization: 
% u(x)=u_0 + 2\sum_{k\geq 1} u_k T_k(2*x/(b-a)+(b+a)/(b-a)).

if nargin < 4
    nb = 1e2;
end
if nargin < 2
    a = -1;
    b = 1;
end

K = length(u)-1;

u(2:end)=2*u(2:end);

X = (a:(b-a)/nb(1):b)';
MK = cos(acos(2*X/(b-a)-(b+a)/(b-a))*(0:K));%Very bad way of coding this
eval = MK*u;

% if min(eval) < -1e-5
%     eval = -eval;
%     u = -u;
% end

plot(X,eval,'Linewidth',2)
    