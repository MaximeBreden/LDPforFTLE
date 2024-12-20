clear variables
close all
clc

format short

%%% /!\ This script will run without Intlab, but of course the obtained
%%% results are then not guaranteed to be immune to rounding errors. 
%%% In order to get guaranteed results, start Intlab in your current
%%% directory, and put the use_Intlab variable below to true.

use_Intlab = true;

%% Initialization

option = 'rate_function'; 
% option = 'asymptotic_variance'; 
%%% Change option between 'rate_function' and 'asymptotic_variance' in
%%% order to get all the results stated in the paper. We do not do both at
%%% once (although in principle we could) because each requires rather 
%%% different ranges of p.

% Parameters of the model, initialized as string to be safe when using interval arithmetic
alpha = '1';
sigma = '1';

N = 20; % number of Fourier modes used
K = 10; % number of Chebyshev modes used

% Weights for the proof, initialized as string to be safe when using interval arithmetic
nu = '1.0';
eta = '1.01';

if strcmp(option,'rate_function')
    nb = 50; % number of points in b
    tab_b = logspace(log10(2.1-2), log10(20-2), nb) + 2;
    tab_p_mid = 1./(tab_b-2).^2;
    p_rad = 5;
    p_rad_center = 0.1;
    % We study the MLE for p in [p_mid-p_rad,p_mid+p_rad]
    % The value of p_rad is used to numerically find the minimum, and then
    % p_rad_center is used for the computation of I(0)
    tab_I0 = zeros(1,nb);
    if exist('intval','file') && use_Intlab
        tab_iI0 = intval(tab_I0);
    end
elseif strcmp(option,'asymptotic_variance')
    nb = 50; % number of points in b
    tab_b = linspace(0, 20, nb);
    tab_p_mid = zeros(1,nb);
    p_rad = 2^-4;
    % We study the MLE for p in [-p_rad,+p_rad]
    tab_LE = zeros(1,nb);
    tab_AV = zeros(1,nb);
    if exist('intval','file') && use_Intlab
        tab_iLE = intval(tab_LE);
        tab_iAV = intval(tab_AV);
    end
else
    error("invalid choice for option")
end

alpha = str2double(alpha);
sigma = str2double(sigma);
problem.alpha = alpha;
problem.sigma = sigma;
problem.continuation.p_rad = p_rad;

theta_grid = ((K:-1:0)*pi)/K;
cheb_grid = cos(theta_grid); % Chebyshev grid on [-1,1]

for n = 1:nb
    fprintf("\n")
    b = tab_b(n)
    p_mid = tab_p_mid(n);
    problem.b = b;
    problem.continuation.p_mid = p_mid;
    problem.continuation.p_rad = p_rad;
    p_grid = p_mid + cheb_grid*p_rad; % Chebyshev grid on [p_mid-p_rad,p_mid+p_rad]
    problem.continuation.p_grid = p_grid;

    %% Numerical computation of an eigenpair (which should correspond to Lambda(k))
    
    show = false;
    X_grid = compute_eigenpair(K,N,problem,show); % Approximate eigenpairs, for each k in p_grid
   
    if strcmp(option,'rate_function')      
        %% Recentering around the numerically detected minimizer within [p_mid-p_rad,p_mid+p_rad]
        lambda_grid = X_grid(end,:);
        lambda_coeffs = grid2coeffs(lambda_grid);
        pts_Cheb = linspace(-1,1,1e4);
        [~,index] = min(eval_Cheb(lambda_coeffs,pts_Cheb)); % Finding where Lambda takes its minimum
        if index == 1 || index == length(pts_Cheb)
            fprintf("The minimum of Lambda (if it exists) may not be within the selected range of p's. No recentering done.\n")
        else
            % Redefining p_mid and p_rad
            p_mid = p_mid+pts_Cheb(index)*p_rad;
            % fprintf("Recentering the values of p, the new p_mid is %g\n",p_mid)
            p_grid = p_mid + cheb_grid*p_rad_center;
            problem.continuation.p_mid = p_mid;
            problem.continuation.p_rad = p_rad_center;
            problem.continuation.p_grid = p_grid;
            X_grid = compute_eigenpair(K,N,problem,show); % Approximate eigenpairs, for each p in the new p_grid
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
        iproblem.continuation.p_mid = intval(problem.continuation.p_mid);
        iproblem.continuation.p_rad = intval(problem.continuation.p_rad);
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
        N_FFT = 2^12;
        K_FFT = 2^2;
        N_FFT_max = 2^15;
        K_FFT_max = 2^12;
        iu_coeffs = intval(u_coeffs);
        true_proof = prove_positivityFFT(iu_coeffs,irmin,N_FFT,K_FFT,N_FFT_max,K_FFT_max);
    end
    
    %% Further computations using the moment Lyapunov function
    lambda_coeffs = real(X_coeffs(end,:));
    if strcmp(option,'rate_function')
        %% Computing I(0)        
        if true_proof
            [I0,iI0] = compute_I0(lambda_coeffs,iproblem,irmin);
            tab_iI0(n) = iI0;
        else
            I0 = compute_I0(lambda_coeffs,problem);
            tab_I0(n) = I0;
        end
    elseif strcmp(option,'asymptotic_variance')
        %% Computing the asymptotic LE and the asymptotic variance        
        dlambda_coeffs = 1/p_rad * derCheb(lambda_coeffs);
        d2lambda_coeffs = 1/p_rad * derCheb(dlambda_coeffs);
        if not( p_mid-p_rad < 0 && p_mid+p_rad>0 ) 
            fprintf("\n0 is not contained within the selected range of p's, unable to compute the asymptotic LE and the asymptotic variance\n")
        else
            LE = eval_Cheb(dlambda_coeffs,(0-p_mid)/p_rad);
            tab_LE(n) = LE;
            fprintf("\nApproximate value of the asymptotic LE: %g\n",LE)
            if true_proof
                ip_rad = intval(p_rad);
                fprintf("LE rigorously enclosed:\n")
                ilambda_coeffs = intval(lambda_coeffs);
                idlambda_coeffs = 1/ip_rad * derCheb(ilambda_coeffs);
                i0_resc = (0-p_mid)/ip_rad;
                iLE = eval_Cheb(idlambda_coeffs,i0_resc) + 1/ip_rad * constant_derivative(ieta,1,inf(i0_resc),sup(i0_resc)) * midrad(0,irmin);
                tab_iLE(n) = iLE;
                infsup(iLE)
            end
                    
            AV = eval_Cheb(d2lambda_coeffs,(0-p_mid)/p_rad);
            tab_AV(n) = AV;
            fprintf("\nApproximate value of the asymptotic variance: %g\n",AV)
            if true_proof
                fprintf("AV rigorously enclosed:\n")
                id2lambda_coeffs = 1/ip_rad * derCheb(idlambda_coeffs);
                iAV = eval_Cheb(id2lambda_coeffs,i0_resc) + 1/ip_rad^2 * constant_derivative(ieta,2,inf(i0_resc),sup(i0_resc)) * midrad(0,irmin);
                tab_iAV(n) = iAV;
                infsup(iAV)
            end
        end
    end
end

%% Plots
if strcmp(option,'rate_function')
    figure
    if true_proof
        plot(tab_b,tab_iI0,'.k','Markersize',15)
        title('Rigorous computation of the rate function','Interpreter','Latex')
        fprintf("\nThe error bound for each value of I(0) is a most:\n%g\n",max(rad(tab_iI0)))
    else
        plot(tab_b,tab_I0,'k','Linewidth',2)
        title('Non-rigorous computation of the rate function','Interpreter','Latex')
    end
    xlabel('$b$','Interpreter','Latex')
    ylabel('$\mathcal{I}_{b}(0)$','Interpreter','Latex')
    xticks([2 5 10 15 20])
    xticklabels({'2','5','10','15','20'})
    set(gca,'FontSize',15) 
    axis tight
    xlim([2 20])

elseif strcmp(option,'asymptotic_variance')
    figure
    if true_proof
        plot(tab_b,mid(tab_iLE),'.k','Markersize',15)
        title('Rigorous computation of the asymptotic Lyapunov exponent','Interpreter','Latex')
        fprintf("\nThe error bound for each asymptotic Lyapunov exponent is a most:\n%g\n",max(rad(tab_iLE)))
    else
        plot(tab_b,tab_LE,'k','Linewidth',2)
        title('Non-rigorous computation of the asymptotic Lyapunov exponent','Interpreter','Latex')
    end
    xlabel('$b$','Interpreter','Latex')
    ylabel('$\Lambda_b''(0)$','Interpreter','Latex')
    set(gca,'FontSize',15) 
    axis tight

    figure
    if true_proof
        plot(tab_b,mid(tab_iAV),'.k','Markersize',15)
        title('Rigorous computation of the asymptotic variance ','Interpreter','Latex')
        fprintf("\nThe error bound for each asymptotic variance is a most:\n%g\n",max(rad(tab_iAV)))
    else
        plot(tab_b,tab_AV,'k','Linewidth',2)
        title('Non-rigorous computation of the asymptotic variance','Interpreter','Latex')
    end
    xlabel('$b$','Interpreter','Latex')
    ylabel('$\Lambda_b''''(0)$','Interpreter','Latex')
    set(gca,'FontSize',15) 
    axis tight
end