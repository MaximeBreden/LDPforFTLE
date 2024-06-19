%%% Non rigorous computation of I(0) = - inf_p Lambda(p) for different 
%%% values of alpha. When alpha_min is close to 0, the output gets 
%%% unreliable unless K is large enough (or unless the domain [a,b] is 
%%% small enough). For rigorous computations, see script_pitchfork_single.m
%%% and script_pitchfork_all.m.

clear variables
close all
clc

format short

%% Initialization
% Parameters of the original SDE (alpha, sigma)
alpha_min = 0.001;
alpha_max = 10;
% tab_alpha = linspace(alpha_min,alpha_max,200);
tab_alpha = logspace(log10(alpha_min),log10(alpha_max),50);
sigma = 1;

prblm.sigma = sigma;

% Bounded domain [a,b] on which we restrict ourselves for the
% eigenfunctions
a = -5;
b = 5;

prblm.a = a;
prblm.b = b;

% Number of Chebyshev modes used to compute approximate eigenfunctions
K = 200;

% Parameter for the eigenvalue computation. Using Matlab's eigs with the
% option for directly getting a few eigenvalues with largest real part does
% not seem to work well for this problem. We instead compute eigenvalues
% which are the closest to 0 (nb_eigs_SM of them), and then take the one
% with largest real part among those. Of course this could in principle
% miss the eigenvalue with largest real part, but this never happened for
% our examples, and this procedure is only used in the nonrigorous part of 
% the code.
global nb_eigs_SM
nb_eigs_SM = 10;

%% Approximation of I(0) as a function of alpha
pmin = 0;
pmax = 1e3;
tol = 1e-6;
it_max = 1e3;
shift = 10;

[tab_I0,tab_pstar] = compute_rate_function_numerics(prblm,K,tab_alpha,pmin,pmax,tol,it_max,shift);

%% Figures
figure
plot(tab_alpha,tab_pstar,'Linewidth',2)
xlabel('$\alpha$','Interpreter','Latex')
ylabel('$p^*_\alpha$','Interpreter','Latex')
title('Non rigorous computation of the minimizer of $\Lambda_{\alpha}$','Interpreter','Latex')
set(gca,'FontSize',15)   
axis tight

figure
plot(tab_alpha,tab_I0,'Linewidth',2)
xlabel('$\alpha$','Interpreter','Latex')
ylabel('$\mathcal{I}_{\alpha}(0)$','Interpreter','Latex')
title('Non rigorous computation of the rate function','Interpreter','Latex')
set(gca,'FontSize',15)     
axis tight

if max(tab_pstar) > 100
    figure
    semilogy(tab_alpha,tab_pstar,'Linewidth',2)
    xlabel('$\alpha$','Interpreter','Latex')
    ylabel('$p^*_\alpha$','Interpreter','Latex')
    title('Non rigorous computation of the minimizer of $\Lambda_{\alpha}$','Interpreter','Latex')
    set(gca,'FontSize',15)  
    axis tight
end

if max(tab_I0) > 100
    figure
    semilogy(tab_alpha,tab_I0,'Linewidth',2)
    xlabel('$\alpha$','Interpreter','Latex')
    ylabel('$\mathcal{I}_{\alpha}(0)$','Interpreter','Latex')
    title('Non rigorous computation of the rate function','Interpreter','Latex')
    set(gca,'FontSize',15)   
    axis tight
end

% save('data_plot_nonrig.mat','tab_alpha','tab_I0')