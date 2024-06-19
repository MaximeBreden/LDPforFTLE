clear variables
close all
clc

format long

%%% /!\ This script will not run without Intlab.

%% Initialization

% Parameters of the original SDE
isigma = intval('1');
ialpha = intval('1');

% Bounded domain [a,b] on which we restrict ourselves for the numerics
a = -4;
b = 4;

% Number of Chebyshev modes used to compute approximate eigenfunctions
K = 100;

% Interval of p in which we look for the minimizer of Lambda
pmin = 0.01;
pmax = 2;

% The interval containing the minimizer will be subdivided into nb_p-1 pieces of (approximately) equal length.
% Must be >= 2. 
nb_p = 20; 

% Coefficient A in the potential V_0 = Ax^2+B used for the base
% problem (B is then computed automatically later so that V_0 <= V)
A = 1; 

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

%% Numerically finding an approximate minimizer for Lambda
tol = 1e-6;
it_max = 1e3;
nb_p_num = 101;
pstar = find_argmin_Lambda(prblm,K,pmin,pmax,nb_p_num,tol,it_max);

%% Rigorously enclosing the minimizer of Lambda
prec = 5;
nb_dec = round(log10(pstar));
ipstar = intval(num2str(pstar,prec)); %only used to select a pstar which is easier to display in the paper
pstar = mid(ipstar);
prblm.p = pstar;
prblm.ip = ipstar;
show = 1;
iLambdapstar = compute_Lambda(prblm, K, A, show); % rigorous enclosure of Lambda at pstar

nb_dec_dp = 5 - round(log10(pstar));
delta_p = 10^(-nb_dec_dp);
delta_p_max = 10^round(log10(pstar));
success = false;

while not(success) && delta_p <= delta_p_max
    ipleft = intval(num2str(pstar-delta_p,prec));
    prblm.p = mid(ipleft);
    prblm.ip = ipleft;
    iLambdapleft = compute_Lambda(prblm, K, A); % rigorous enclosure of Lambda at some pleft < pstar
    
    ipright = intval(num2str(pstar+delta_p,prec));
    prblm.p = mid(ipright);
    prblm.ip = ipright;
    iLambdapright = compute_Lambda(prblm, K, A); % rigorous enclosure of Lambda at some pright > pstar
    
    % If Lambda(pstar) <= Lambda(pleft) and Lambda(pstar) <= Lambda(pright)
    % by convexity of Lambda the minimum of Lambda must lie between
    % pleft and pright.
    success = (iLambdapstar <= iLambdapleft) && (iLambdapstar <= iLambdapright);
    delta_p = 10*delta_p;
end

if success
    ipleft
    infsup(iLambdapleft)
    ipstar
    infsup(iLambdapstar)
    ipright
    infsup(iLambdapright)
    ipstar = infsup(inf(ipleft),sup(ipright));
    fprintf("\n\nThe minimizer of Lambda is included in:\n")
    infsup(ipstar)
else
    fprintf("\n\nUnable to rigorously enclose the minimizer of Lambda\n")
    return
end


%% Rigorously enclosing the minimal value of Lambda

% By taking np>2 we can subdivide the interval ipstar containing the 
% minimizer pstar into smaller subintervals, in order to get a tighter 
% enclosure of Lambda(ipstar)
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
fprintf("\n\nI(0) is included in:\n")
infsup(-iLambda_min)




