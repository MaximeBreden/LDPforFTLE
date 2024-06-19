function [iLambdap,prblm] = compute_Lambda(prblm, K, A, show, para_homotopy, nb_sub_for_B)

% Rigorous computation of Lambda(p) the first eigenvalue of H_p, using the
% Rayleigh-Ritz and Lehmann-Maehly bounds together with the homotopy 
% method.
% The field prblm contains the paramters alpha and sigma of the SDE, as
% well as the spatial domain [a,b] used to approximate the problem. The
% approximate eigenvectors are computed as Chebyshev polynomials of order
% K. Subsequent inputs are optionals.

% A is used to define the base problem (for which we use a quadratic
% potential V0(x) = Ax^2 + B). Given A, we then find B such that V0 <= V (V
% being the original order 6 potential appearing in H_p). This B is 
% selected by optimizing over a sample of size nb_sub_for_B.
% show can be used to display intermediate information during the homotopy
% method. 
% para_homotopy is a field containing there parameters allowing to
% fine-tune some of the homotopy steps.


global nb_eigs_SM

if nargin < 6
    nb_sub_for_B = 1e3;
    if nargin < 5 
        para_homotopy.tol_simple_eig = 10^-7;%For the validation of simple VS clustered eigenvalues
        para_homotopy.show_s = 0;%Put to 1 to display what the value of s is during the whole process
        para_homotopy.show_enclosure = 0;%Put to 1 to display the obtained enclosures at each homotopy breakpoint
        if nargin < 4
            show = false;
            if nargin < 3
                A = 1;
            end
        end
    end
end

%% Computation of a suitable base problem

iA = intval(A);
prblm.iA = iA;
prblm.A = A;

B = findB(prblm,nb_sub_for_B); % Rigorously finds a B (which is close to the largest possible) such that V0 <= V
iB = intval(B);
prblm.iB = iB;
prblm.B = B;

if show > 0
    % Plot displaying the original potential together with the one of the base
    % problem
    alpha = prblm.alpha;
    sigma = prblm.sigma;
    p = prblm.p;
    a = prblm.a;
    b = prblm.b;
    
    gamma = @(x) (alpha*x-x.^3).^2/(2*sigma^2) + (1/2-p)*(alpha-3*x.^2);
    g = @(x) A*x.^2+B;
    X = linspace(a,b,1e4);
    figure
    plot(X,gamma(X),'b',X,g(X),'r','Linewidth',2)
    xlabel('x')
    legend('$V$','$V^0$','Interpreter','Latex')
    title('Potentials for $H$ and for the base problem $H^0$','Interpreter','Latex')
    set(gca,'FontSize',15) 
    axis tight
    drawnow
end

%% Initialization of all the matrices needed for the homopty
op = initialize_H(prblm,K); % Construct the discretized version of H_p

%% First (nonrigorous) computation of the principal eigenvalue of H, together with a rigorous computation of the first M eigenvalues of H0
[EVect_H,EVal_H] = eigs(op.Trunc*op.H*op.Mat_BC,op.Trunc*op.Mat_BC,nb_eigs_SM,'SM');
[EVal_H,ind] = sort(real(diag(EVal_H)));
PEVal_H = EVal_H(1);

if show > 1
    PEVect_H = op.Mat_BC*EVect_H(:,ind(1));
    if PEVect_H(1)<0
        PEVect_H = -PEVect_H;
    end
    figure
    plotCheb(PEVect_H,a,b,1e3);
    xlabel('x')
    title(["Approximate principal eigenfunction of H, corresponding eigenvalue: ",num2str(PEVal_H)])
    fprintf("\n\nApproximate principal eigenvalue of H: PEVal_H = %g\n",PEVal_H)
end

% Determining the minimal number of required eigenvalues of the base
% problem
sigma = prblm.sigma;
isigma = prblm.isigma;
margin_nb_eig = 1;
M = ceil((PEVal_H-op.shift-B+margin_nb_eig)/(sigma*sqrt(2*A))+1/2);
EVal_H0 = iB + ((0:M-1)'+1/2)*isigma*sqrt(2*iA) + op.ishift;
if show > 10
    fprintf("\nFirst eigenvalues of H0:\n")
    EVal_H0
    fprintf("We kept the first %i eigenvalues of the base problem\n",M)
end

% Computations of the associated eigenvectors of the base problem (does not
% have to be rigorous)
[EVect_H0_num,EVal_H0_num] = eigs(op.Trunc*op.H0*op.Mat_BC,op.Trunc*op.Mat_BC,min(5*M,K-3),'SM');
[EVal_H0_num,ind0] = sort(real(diag(EVal_H0_num)));

EVal_H0_num = EVal_H0_num(1:M);
EVect_H0_num = op.Mat_BC*EVect_H0_num(:,ind0(1:M));

% Sanity check
rel_err = abs((EVal_H0-EVal_H0_num)./EVal_H0);
if max(rel_err)>1e-3
    fprintf('\nFirst eigenvalues of the base problem, theoretical and numerical:\n')
    [EVal_H0.mid EVal_H0_num]
    figure
    plotCheb(EVect_H0_num(:,end),prblm.a,prblm.b,1e3);
    title("Last eigenfunction of H0 kept")
    error("Something is wrong with the eigenvalues/eigenfunctions of the base problem, you might need a larger domain and/or more Chebyshev modes. You could also try to increase A.")
end

%% Homotopy to get the first eigenvalue of H
if show > 10
    fprintf('\n\nHomotopy between H0 and H...\n')
end
eigenval_0 = mid(EVal_H0);
ieigenval_0 = EVal_H0;
eigenvect_0 = EVect_H0_num;
s = 0;
steps_s = 0;
margin_H = margin_for_homotopy(eigenval_0);%This parameter controls how the breakpoints s_k for the homotopy are chosen 
if show > 10
    fprintf('\nmargin_H  =  %f\n',margin_H)
end
% tol_simple_eig = 10^-7;%For the validation of simple VS clustered eigenvalues
% show_s = 1;%Put to 1 to display what the value of s is during the whole process
% show_enclosure = 1;%Put to 1 to display the obtained enclosures at each homotopy breakpoint

%The homotopy method 
% tic;
[ilower_eig_H,iupper_eig_H,~,~] = homotopy(s,steps_s,op,eigenval_0,eigenvect_0,ieigenval_0.inf,margin_H,para_homotopy);
% toc;

iLambdap = -(infsup(ilower_eig_H(1),iupper_eig_H(1))-op.ishift); %rigorous enclosure of the first eigenvalue of H
if show > 10
    fprintf("\nHomotopy succesful, here is the enclosure for Lambda(p):\n")
    iLambdap
    fprintf('radius of the enclosure: %g\n',rad(iLambdap))
end