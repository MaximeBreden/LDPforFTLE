function [rmin,rmax] = proof_eigenproblem(X_coeffs,A_coeffs,problem)

% Computation of the bounds needed for the Newton-Kantorovich argument

%% Init
sz = size(X_coeffs);
N = (sz(1)-2)/2; 
K = sz(2)-1;
[alpha,b,sigma] = obtain_parameters(problem);
nu = problem.proof.nu;
eta = problem.proof.eta;
p_mid = problem.continuation.p_mid;
p_rad = problem.continuation.p_rad;
    
if exist('intval','file') && isintval(X_coeffs(1))
    ipi = intval('pi');
    K2 = 2^(ceil(log2(2*K)))+1; % K2-1 is the smallest power of two >= 2*K
    K3 = 2^(ceil(log2(3*K)))+1; % K3-1 is the smallest power of two >= 3*K
else
    ipi = pi;
    K2 = 2*K+1;
    K3 = 3*K+1;
end


%% Y

AFX_coeffs_K3 = f2i(zeros(2*N+4,K3),X_coeffs(1));
Ind_N_finite = [2:2*N+2, 2*N+4];
Ind_N_spill = [1, 2*N+3];

%%% F(X) %%%
X_coeffs_K2 = [X_coeffs,zeros(2*N+2,K2-(K+1))];
problem.data_type = 'coeffs';
problem.p = [p_mid,p_rad/2,zeros(1,K2-2)];
problem.normalization = X_coeffs_K2(1:end-1,:);
FX_coeffs_K2 = F_eigenproblem(X_coeffs_K2,problem,true); %F(X) in coefficients
FX_grid_K3 = coeffs2grid([FX_coeffs_K2,zeros(2*N+4,K3-K2)]); %F(X) on the extended grid


%%% Finite part (in Fourier) %%%
A_coeffs_K3 = cat(3, A_coeffs, zeros(2*N+2,2*N+2,K3-(K+1)));
A_grid_K3 = coeffs2grid(A_coeffs_K3);
AFX_grid_K3 = f2i(zeros(2*N+2,K3),X_coeffs(1));
for j = 1:K3
    AFX_grid_K3(:,j) = A_grid_K3(:,:,j) * FX_grid_K3(Ind_N_finite,j);
end
AFX_coeffs_K3(Ind_N_finite,:) = grid2coeffs(AFX_grid_K3);

%%% Spillover terms (in Fourier) %%%
if mod(K3-1,K2-1)==0
    r = round((K3-1)/(K2-1));
    FX_grid_K2 = FX_grid_K3(:,1:r:end);
else
    FX_grid_K2 = coeffs2grid(FX_coeffs_K2);
end
Aspill = diag(-1./(2*sigma^2*[-(N+1),N+1].^2));
AFX_grid_K2_tail = Aspill * FX_grid_K2(Ind_N_spill,:);
AFX_coeffs_K3(Ind_N_spill,1:K2) = grid2coeffs(AFX_grid_K2_tail);

%%% norm computation %%%
Y = norm_nueta(AFX_coeffs_K3,nu,eta);
disp(['Y = ',num2str(i2f(Y))])


%% Z1

ADFX_coeffs_K2 = f2i(zeros(2*N+6,2*N+4,K2),X_coeffs(1));
Ind_N_finite = [3:2*N+3, 2*N+6];
Ind_N_spill = [1:2, 2*N+4:2*N+5];

%%% DF(X) %%%
X_grid_K2 = coeffs2grid(X_coeffs_K2);
cheb_grid_K2 = cos( ((K2-1:-1:0)*ipi)/(K2-1) );
p_grid_K2 = p_mid + cheb_grid_K2*p_rad;
problem.data_type = 'grid';
problem.p = p_grid_K2;
problem.normalization = X_grid_K2(1:end-1,:);
DFX_grid_K2 = DF_eigenproblem(X_grid_K2,problem,true);


%%% Finite part (in Fourier) %%%
if mod(K3-1,K2-1)==0
    r = round((K3-1)/(K2-1));
    A_grid_K2 = A_grid_K3(:,:,1:r:end);
else
    A_grid_K2 = coeffs2grid(A_coeffs_K3(:,:,1:K2));
end
ADFX_grid_K2 = f2i(zeros(2*N+2,2*N+4,K2),X_coeffs(1));
for j = 1:K2
    ADFX_grid_K2(:,:,j) = A_grid_K2(:,:,j) * DFX_grid_K2(Ind_N_finite,:,j);
end
ADFX_coeffs_K2(Ind_N_finite,:,:) = grid2coeffs(ADFX_grid_K2);

%%% Spillover part (in Fourier) %%%
Aspill = -1./(2*sigma^2*[-(N+2):-(N+1),N+1:N+2].^2)';
Aspill = repmat(Aspill,[1,2*N+4,K2]);
ADFX_grid_K2_tail = Aspill .* DFX_grid_K2(Ind_N_spill,:,:);
ADFX_coeffs_K2(Ind_N_spill,:,:) = grid2coeffs(ADFX_grid_K2_tail);

%%% norm computation %%%
Iext = spdiags(ones(2*N+3,1),-1,2*N+5,2*N+3);
Iext = [Iext, zeros(2*N+5,1);
        zeros(1,2*N+3), 1];
B = cat(3,full(Iext),zeros(2*N+6,2*N+4,K2-1)) - ADFX_coeffs_K2;
Z1_finite = norm_nueta(B, nu, eta);
disp(['Z1_finite = ',num2str(i2f(Z1_finite))])

%%% Tail part (in Fourier) %%%
p = f2i(zeros(1,K+1),X_coeffs(1));
p(1) = p_mid;
p(2) = 1/2*p_rad;
lambda_coeffs = X_coeffs(2*N+2,:);
lpap = lambda_coeffs + alpha*p;
Z1_tail = 1/(2*sigma^2) * ( 1/(N+2)^2 * (abs(b)*(N+2)+norm_nueta(lpap,nu,eta)) + ...
                            1/(2*(N+1)^2) * abs(b) * (2*(N+2)+norm_nueta(p,nu,eta)) * nu);

disp(['Z1_tail = ',num2str(i2f(Z1_tail))])
      
Z1 = max(Z1_finite,Z1_tail);
disp(['Z1 = ',num2str(i2f(Z1))])
   
%% Z2
Z2 = max( norm_nueta(A_coeffs,nu,eta), 1/(2*(sigma*(N+1))^2) );
disp(['Z2 = ',num2str(i2f(Z2))])


%% Cheking that we have a contraction

Delta = (1-Z1)^2-2*Y*Z2;
if i2f(Z1)<1 && i2f(Delta,'inf')>0
    rmin = i2f( (1-Z1-sqrt(Delta)) / Z2, 'sup' );
    rmax = i2f( (1-Z1) / Z2, 'inf' );
    disp(['[r_min, r_max] = [',num2str(rmin),', ',num2str(rmax),']'])
    fprintf('Validation successful :)\n\n')
else
    rmin = NaN;
    rmax = NaN;
    fprintf('The validation failed :(\n\n')
end