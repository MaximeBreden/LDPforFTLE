clear variables
close all
clc

format short

%%% /!\ This script will not run without Intlab.

%% Initialization

% Parameter sigma of the original SDE
isigma = intval('1'); 

% Values of alpha at which we are going to rigorously compute the rate
% function at 0
alpha_min = 0.001;
alpha_max = 10;
nb = 50;
tab_alpha = logspace(log10(alpha_min),log10(alpha_max),nb);

% For each of these values of alpha, we then select in a somewhat ad-hoc 
% manner the interval [a,b] and the number K of Chebyshev modes used for 
% computing approximate eigenvectors, as well as the interval [pmin,pmax]
% that is considered for values of p when first computing Lambda(p)
% approximately, and the values of A used for defining the quadratic
% potential in the base problem of the homotopy method.

mask1 = (tab_alpha >= 0.1);
mask2 = min( tab_alpha < 0.1, tab_alpha >= 0.01);
mask3 = (tab_alpha < 0.01);

% Values of a and b
tab_b = zeros(1,nb);
tab_b(mask1) = 2*log10(tab_alpha(mask1))+4;
tab_b(mask2) = 0.5 + 1.75*(log10(tab_alpha(mask2))+2);
tab_b(mask3) = 0.2 + 0.3*(log10(tab_alpha(mask3))+3);
tab_a = -tab_b;

% Values of K
tab_K = 60*ones(1,nb);
tab_K(mask1) = 50 + 75*(log10(tab_alpha(mask1))+1);
tab_K = round(tab_K);

% Values of pmin
tab_pmin = 0.01*ones(1,nb);
tab_pmin(mask2) = 0.3./tab_alpha(mask2).^2;
tab_pmin(mask3) = 0.3./tab_alpha(mask3).^2;

% Values of pmax
tab_pmax = zeros(1,nb);
tab_pmax(mask1) = 10./tab_alpha(mask1);
tab_pmax(mask2) = 0.5./tab_alpha(mask2).^2;
tab_pmax(mask3) = 0.5./tab_alpha(mask3).^2;

% Values of A
tab_A = ones(1,nb);
mask_A = min( tab_alpha < 1, tab_alpha >= 0.01);
tab_A(mask_A) = 1./tab_alpha(mask_A).^2;
tab_A(mask3) = 1./tab_alpha(mask3).^2;

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

tab_alpha = intval(tab_alpha);
tab_ipstar = intval(zeros(1,nb));
tab_I0 = intval(zeros(1,nb));
for n = 1:nb
    % Selecting parameters
    ialpha = tab_alpha(n)
    a = tab_a(n);
    b = tab_b(n);
    K = tab_K(n);
    alpha = mid(ialpha);
    sigma = mid(isigma);
    prblm.ialpha = ialpha;
    prblm.isigma = isigma;
    prblm.alpha = alpha;
    prblm.sigma = sigma;
    prblm.ia = intval(a);
    prblm.ib = intval(b);
    prblm.a = a;
    prblm.b = b;
    
    %% Numerically finding an approximate minimizer pstar for Lambda
    tol = 1e-6;
    it_max = 1e3;
    nb_p_num = 0;
    pmin = tab_pmin(n);
    pmax = tab_pmax(n);
    pstar = find_argmin_Lambda(prblm,K,pmin,pmax,nb_p_num,tol,it_max);
    
    %% Rigorously enclosing the minimizer pstar of Lambda 
    prblm.p = pstar;
    prblm.ip = intval(pstar);
    A = tab_A(n);
    show = 0;
    iLambdapstar = compute_Lambda(prblm, K, A, show); % rigorous enclosure of Lambda at pstar
    
    nb_dec_dp = 5 - round(log10(pstar));
    delta_p = 10^(-nb_dec_dp);
    delta_p_max = 10^round(log10(pstar));
    success = false;
    
    while not(success) && delta_p <= delta_p_max
        pleft = pstar - delta_p;
        ipleft = intval(pleft);
        prblm.p = pleft;
        prblm.ip = ipleft;
        iLambdapleft = compute_Lambda(prblm, K, A); % rigorous enclosure of Lambda at some pleft < pstar
        
        pright = pstar + delta_p;
        ipright = intval(pright);
        prblm.p = pright;
        prblm.ip = ipright;
        iLambdapright = compute_Lambda(prblm, K, A); % rigorous enclosure of Lambda at some pright > pstar
        
        % If Lambda(pstar) <= Lambda(pleft) and Lambda(pstar) <= Lambda(pright)
        % by convexity of Lambda the minimum of Lambda must lie between
        % pleft and pright.
        success = (iLambdapstar <= iLambdapleft) && (iLambdapstar <= iLambdapright);
        delta_p = 10*delta_p;
    end
    
    if success
        ipstar = infsup(inf(ipleft),sup(ipright));
        tab_ipstar(n) = ipstar;
    else
        fprintf("\n\nUnable to rigorously enclose the minimizer of Lambda\n")
        return
    end
    
    
    %% Rigorously enclosing the minimal value of Lambda
    
    % By taking np>2 we can subdivide the interval ipstar containing the 
    % minimizer pstar into smaller subintervals, in order to get a tighter 
    % enclosure of Lambda(ipstar)
    
    nb_p = 2;
    tab_p = linspace(inf(ipstar),sup(ipstar),nb_p);
    tab_p(1) = inf(ipstar);
    tab_p(end) = sup(ipstar);
    tab_p = infsup(tab_p(1:end-1),tab_p(2:end));
    tab_Lambda = intval(zeros(1,nb_p-1));
    
    for j = 1:nb_p-1
        prblm.p = mid(tab_p(j));
        prblm.ip = tab_p(j);
        tab_Lambda(j) = compute_Lambda(prblm, K, A); % rigorous enclosure of Lambda on the j-th subinterval
    end
    
    iLambda_min = min(tab_Lambda);
    tab_I0(n) = -iLambda_min;
end

%% Plots
figure
semilogy(tab_alpha.mid,tab_I0.inf,'.b','markersize',10)
hold on
semilogy(tab_alpha.mid,tab_I0.sup,'.r','markersize',10)
xlabel('$\alpha$','Interpreter','Latex')
ylabel('$\mathcal{I}_{\alpha}(0)$','Interpreter','Latex')
title('Rigorous computation of the rate function','Interpreter','Latex')
set(gca,'FontSize',15) 
axis tight







