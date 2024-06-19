function B = findB(prblm,nb)

% Computes the (close to) largest B such that V0 <= V

alpha = prblm.ialpha;
sigma = prblm.isigma;
p = prblm.ip;

A = prblm.iA;

x0 = sqrt( max( alpha+sigma*sqrt( max( 2*(A-3*(p-1/2)), 0 ) ) , 0 ) );
x0 = max(sup(x0),1);

phi = @(x) x.^2 .* ( (x.^2-alpha).^2/(2*sigma^2) + 3*(p-1/2)-A ) + alpha*(1/2-p);

B_asymptotic = 0;
B_local = Inf;

it = 0;
while B_asymptotic <= B_local && it < 10
    B_asymptotic = phi(x0);
    x = linspace(-x0,x0,nb);
    x(1) = -x0;
    x(end) = x0;
    x = infsup(x(1:end-1),x(2:end));
    B_local = min(phi(x));
    it = it + 1;
    x0 = 10*x0;
end

B = inf(min(B_asymptotic,B_local));




