%%% Non rigorous computation of the asymptotic Lyapunov exponent 
%%% lambda(alpha) for different values of alpha. 

clear variables
close all
clc

sigma = 1;
L = 100; % We truncate the real line to [-L,L]
tab_alpha = linspace(-2.6,2,1000);
tab_lambda = 0*tab_alpha;

j = 1;
for alpha = tab_alpha
    p = @(x) exp(1/sigma^2*(alpha*x.^2-1/2*x.^4)); % stationary density
    c = integral(p,-L,L); % normalization
    tab_lambda(j) = alpha - 1/c*integral(@(x) 3*x.^2.*p(x),-L,L);
    j = j+1;
end

figure
plot(tab_alpha,tab_lambda,'Linewidth',2)
xlabel('$\alpha$','Interpreter','Latex')
ylabel('$\lambda(\alpha)$','Interpreter','Latex')
title('Non rigorous computation of the asymptotic Lyapunov exponent','Interpreter','Latex')
set(gca,'FontSize',15) 
axis tight
