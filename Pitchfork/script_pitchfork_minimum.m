clear variables
close all
clc

format short

%%% /!\ This script will not run without Intlab.

%% Initialization

% Parameter sigma of the original SDE
isigma = intval('1');
tab_alpha = [intval('1.225'),intval('1.227'),intval('1.229')];

a = -4;
b = 4;
K = 150;
pmin = 0.01;
pmax = 10;
A = 1;

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

tab_ipstar = intval(zeros(1,3));
tab_I0 = intval(zeros(1,3));
for n = 1:3
    ialpha = tab_alpha(n)
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
    
    %% Numerically finding an approximate minimizer for Lambda
    tol = 1e-6;
    it_max = 1e3;
    nb_p_num = 0;
    pstar = find_argmin_Lambda(prblm,K,pmin,pmax,nb_p_num,tol,it_max);
    
    %% Rigorously enclosing the minimizer of Lambda
    prblm.p = pstar;
    prblm.ip = intval(pstar);
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
    
    % Since the 3 values of alpha are quite close to one another, and close
    % to the minimum where alpha -> I_alpha(0) is flat, we need rather
    % tight enclosures of I_alpha(0) to garantuee that we indeed have a
    % minimum. In each case, we therefore subdivide the interval ipstar 
    % into smaller subintervals, in order to get a tighter enclosure of 
    % Lambda(ipstar).
    nb_p = 2000; % number of subintervals
    tab_p = linspace(inf(ipstar),sup(ipstar),nb_p);
    tab_p(1) = inf(ipstar);
    tab_p(end) = sup(ipstar);
    tab_p = infsup(tab_p(1:end-1),tab_p(2:end));
    tab_Lambda = intval(zeros(1,nb_p-1));
    
    for j = 1:nb_p-1
        fprintf("\nSubinterval %i out of %i\n",j,nb_p-1)
        prblm.p = mid(tab_p(j));
        prblm.ip = tab_p(j);
        tab_Lambda(j) = compute_Lambda(prblm, K, A); % rigorous enclosure of Lambda on the j-th subinterval
    end
    
    iLambda_min = min(tab_Lambda);
    tab_I0(n) = -iLambda_min;
end


figure
semilogy(tab_alpha.mid,tab_I0.sup,'.r','markersize',10)
hold on
semilogy(tab_alpha.mid,tab_I0.inf,'.b','markersize',10)
xlabel('$\alpha$','Interpreter','Latex')
ylabel('$\mathcal{I}_{\alpha}(0)$','Interpreter','Latex')
legend('Upper bounds','Lower bounds')
title('Rigorous computation of the rate function','Interpreter','Latex')
set(gca,'FontSize',15) 
axis tight

%% Proof of the existence of a minimum of alpha -> I_alpha(0)
test = (tab_I0(2) < tab_I0(1)) && (tab_I0(2) < tab_I0(3));
if test
    fprintf("\n\nWe have proven there exists a minimum of I_alpha(0) for alpha in:\n")
    fprintf("[%f,%f]\n",inf(tab_alpha(1)),sup(tab_alpha(3)))
else
    fprintf("\n\nUnable to prove the existence of a minimum of I_alpha(0).\n")
end
    