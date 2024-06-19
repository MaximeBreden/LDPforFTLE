clear variables
close all
clc

format short

%%% /!\ This script will run without Intlab, but of course the obtained
%%% results are then not guaranteed to be immune to rounding errors. 
%%% In order to get guaranteed results, start Intlab in your current
%%% directory, and put the use_Intlab variable below to true.

use_Intlab = false;

%% Initialization

% Parameters of the model, initialized as string to be safe when using interval arithmetic
alpha = '1';
b = '5';
sigma = '1';

% We study the MLE Lambda(p) for p in [p_mid-p_rad,p_mid+p_rad]
p_mid = 1;
p_rad = 5;

N = 20; % number of Fourier modes used
K = 50; % number of Chebyshev modes used

center = false; % put to true if you want to automatically re-center p around the minimum of Lambda
p_rad_center = 0.1; % not used unless center = true

% Weights for the proof, initialized as string to be safe when using interval arithmetic
nu = '1.0';
eta = '1.01';

alpha = str2double(alpha);
b = str2double(b);
sigma = str2double(sigma);
problem.alpha = alpha;
problem.b = b;
problem.sigma = sigma;
problem.continuation.p_mid = p_mid;
problem.continuation.p_rad = p_rad;

theta_grid = ((K:-1:0)*pi)/K;
cheb_grid = cos(theta_grid); % Chebyshev grid on [-1,1]
p_grid = p_mid + cheb_grid*p_rad; % Chebyshev grid on [p_mid-p_rad,p_mid+p_rad]
problem.continuation.p_grid = p_grid;


%% Numerical computation of an eigenpair (which should correspond to Lambda(k))

show = false;
X_grid = compute_eigenpair(K,N,problem,show); % Approximate eigenpairs, for each k in p_grid

nb_Cheb = 101; % Nb of points used in Chebyshev for the plot
nb_Four = 101; % Nb of points used in Fourier for the plot
plot_sols(X_grid,problem,nb_Cheb,nb_Four)
fig_counter = 0;

%% Recentering? 
% If we only care about where Lambda takes its minimum, which is all we
% need in order to compute I(0), we can use center = true. Then, the code
% automatically narrows the interval [pmin,pmax] around the minimum, which
% then facilitates the validation.

if center    
    lambda_grid = X_grid(end,:);
    lambda_coeffs = grid2coeffs(lambda_grid);
    pts_Cheb = linspace(-1,1,1e3);
    [~,index] = min(eval_Cheb(lambda_coeffs,pts_Cheb)); % Finding where Lambda takes its minimum
    if index == 1 || index == length(pts_Cheb)
        fprintf("The minimum of Lambda (if it exists) may not be within the selected range of p's. No recentering done.\n")
    else
        figure(1)
        title("Not validated over that range of p's")
        figure(2)
        title("Not validated over that range of p's")
        % Redefining p_mid and p_rad
        p_mid = p_mid+pts_Cheb(index)*p_rad;
        fprintf("Recentering the values of p, the new p_mid is %g\n",p_mid)
        p_rad = p_rad_center;
        p_grid = p_mid + cheb_grid*p_rad;
        problem.continuation.p_mid = p_mid;
        problem.continuation.p_rad = p_rad;
        problem.continuation.p_grid = p_grid;
        X_grid = compute_eigenpair(K,N,problem,show); % Approximate eigenpairs, for each p in the new p_grid
        plot_sols(X_grid,problem,nb_Cheb,nb_Four)
        fig_counter = 2;
    end
end


%% Validation of the eigenpair

% Computation of (the finite part of) the approximate inverse A, sampled on p_grid
A_grid = zeros(2*N+2,2*N+2,K+1);
for j = 1:K+1
    X = X_grid(:,j);
    problem.normalization = X(1:end-1);
    problem.p = i2f(p_grid(j));
    A_grid(:,:,j) = inv(full(DF_eigenproblem(X,problem)));
end


%%% First ``validation'' without interval arithmetics
problem.proof.nu = str2double(nu);
problem.proof.eta = str2double(eta);
X_coeffs = grid2coeffs(X_grid);
X_coeffs(end,:) = real(X_coeffs(end,:));
A_coeffs = grid2coeffs(A_grid); % The A we use in the proof is defined in term of Chebyshev coefficients (rather than in terms of values at p_grid)
fprintf('\nValidation of the eigenpairs (without interval arithmetic)\n')
[rmin,rmax] = proof_eigenproblem(X_coeffs,A_coeffs,problem);
if isnan(rmin)
    return
end

%%% True validation with interval arithmetics
true_proof = false;
if exist('intval','file') && use_Intlab
    inu = intval(nu);
    ieta = intval(eta);
    iX_coeffs = intval(X_coeffs);
    iA_coeffs = intval(A_coeffs);
    iproblem.alpha = intval(alpha);
    iproblem.b = intval(b);
    iproblem.sigma = intval(sigma);
    iproblem.proof.nu = inu;
    iproblem.proof.eta = ieta;
    iproblem.continuation.p_mid = intval(p_mid);
    iproblem.continuation.p_rad = intval(p_rad);
    fprintf('True proof for the eigenpairs (with interval arithmetic)\n')
    [irmin,irmax]=proof_eigenproblem(iX_coeffs,iA_coeffs,iproblem);
    if not(isnan(irmin))
        true_proof = true;
    end
else
    fprintf("You need Intlab in order to run the rigorous proof\n\n")
end


%% Proving that we have the correct eigenpairs (by checking that the eigenfunctions are positive)

u_coeffs = X_coeffs(1:end-1,:);
u_min = min( eval_FourCheb( u_coeffs, linspace(-1,1,1e3), linspace(0,2*pi,1e3), true ), [], 'all') ;
fprintf("Approximate minimal value of the eigenfunctions: %g\n",u_min)
if u_min>0
    fprintf("It seems we computed the correct eigenpairs.\n")
else
    fprintf("Careful, we may not have the correct eigenpairs\n")
end

if true_proof
    fprintf("Currently trying to check positivity rigorously ...\n") 
    N_FFT = 2^7;
    K_FFT = 2^7;
    N_FFT_max = 2^12;
    K_FFT_max = 2^12;
    iu_coeffs = intval(u_coeffs);
    true_proof = prove_positivityFFT(iu_coeffs,irmin,N_FFT,K_FFT,N_FFT_max,K_FFT_max);
end

%% Computing I(0)

format long

lambda_coeffs = real(X_coeffs(end,:));
if true_proof
    [I0,iI0] = compute_I0(lambda_coeffs,iproblem,irmin);
else
    I0 = compute_I0(lambda_coeffs,problem);
end

%% Computing the asymptotic LE and the asymptotic variance

dlambda_coeffs = 1/p_rad * derCheb(lambda_coeffs);
d2lambda_coeffs = 1/p_rad * derCheb(dlambda_coeffs);

if not( p_mid-p_rad < 0 && p_mid+p_rad>0 ) 
    fprintf("\n0 is not contained within the selected range of p's, unable to compute the asymptotic LE and the asymptotic variance\n")
else
    LE = eval_Cheb(dlambda_coeffs,(0-p_mid)/p_rad);
    fprintf("\nApproximate value of the asymptotic LE: %g\n",LE)
    if true_proof
        ip_rad = intval(p_rad);
        fprintf("LE rigorously enclosed:\n")
        ilambda_coeffs = intval(lambda_coeffs);
        idlambda_coeffs = 1/ip_rad * derCheb(ilambda_coeffs);
        i0_resc = (0-p_mid)/ip_rad;
        iLE = eval_Cheb(idlambda_coeffs,i0_resc) + 1/ip_rad * constant_derivative(ieta,1,inf(i0_resc),sup(i0_resc)) * midrad(0,irmin);
        infsup(iLE)
    end
            
    AV = eval_Cheb(d2lambda_coeffs,(0-p_mid)/p_rad);
    fprintf("\nApproximate value of the asymptotic variance: %g\n",AV)
    
    if true_proof
        fprintf("AV rigorously enclosed:\n")
        id2lambda_coeffs = 1/ip_rad * derCheb(idlambda_coeffs);
        iAV = eval_Cheb(id2lambda_coeffs,i0_resc) + 1/ip_rad^2 * constant_derivative(ieta,2,inf(i0_resc),sup(i0_resc)) * midrad(0,irmin);
        infsup(iAV)
    end
end

%% More plots
pts_Cheb = linspace(-1,1,1e3);
figure
plot(p_mid+pts_Cheb*p_rad,eval_Cheb(dlambda_coeffs,pts_Cheb), 'linewidth', 2)
xlabel('$p$', 'Interpreter', 'latex')
ylabel("$\Lambda'(p)$", 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight


figure
plot(p_mid+pts_Cheb*p_rad,eval_Cheb(d2lambda_coeffs,pts_Cheb), 'linewidth', 2)
xlabel('$p$', 'Interpreter', 'latex')
ylabel("$\Lambda''(p)$", 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight

%% Full rate function (NOT RIGOROUS)

temp = eval_Cheb(dlambda_coeffs,[-1;1]);
a_min = temp(1);
a_max = temp(2);

a_mid = (a_min+a_max)/2;
a_rad = 0.99*(a_max-a_min)/2;
problem.rate.a_mid = a_mid;
problem.rate.a_rad = a_rad;
S = 64;

nb_a = 1e3;
q_grid = approx_inverse(dlambda_coeffs,S,problem,nb_a); % Initial guess
q_coeffs = grid2coeffs(q_grid);
it_max = 1e2;
tol = 1e-10;
show = false;
q_coeffs = Newton_rate(q_coeffs,dlambda_coeffs,d2lambda_coeffs,problem,tol,it_max,show);

%%% Plots 
pts_a = linspace(-1,1,1e3);
a_grid = a_mid+pts_a*a_rad;

figure
qa_grid = real(eval_Cheb(q_coeffs,(a_grid-a_mid)/a_rad));
rate_grid = a_grid .* qa_grid - eval_Cheb(lambda_coeffs,(qa_grid-p_mid)/p_rad);
plot(a_grid,rate_grid, 'linewidth', 2)
xlabel('$r$', 'Interpreter', 'latex')
ylabel("$\mathcal{I}(r)$", 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight

min_I = min(rate_grid);
if min_I ~= rate_grid(1) && min_I ~= rate_grid(end)
    T = 10;
    figure
    plot(a_grid, exp(-T*rate_grid), 'linewidth', 2)
    xlabel('$r$', 'Interpreter', 'latex')
    ylabel('$e^{-T\mathcal{I}(r)}$', 'Interpreter', 'latex')
    title(['T = ',num2str(T)], 'Interpreter', 'latex')
    set(gca,'FontSize',15)
    axis tight
end

figure(1+fig_counter)
figure(2+fig_counter)
