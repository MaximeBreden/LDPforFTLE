function [tab_I0,tab_pstar] = compute_rate_function_numerics(prblm,K,tab_alpha,pmin,pmax,tol,it_max,shift)

%% Init
if nargin < 4 || isempty(pmin)
    pmin = -10;
end
if nargin < 5 || isempty(pmax)
    pmin = 1000;
end
if nargin < 6 || isempty(tol)
    tol = 1e-6;
end
if nargin < 7 || isempty(it_max)
    it_max = 1e3;
end
if nargin < 8 || isempty(shift)
    shift = 10;
end

global nb_eigs_SM

% Parameters of the SDE
sigma = prblm.sigma;
% Domain
a = prblm.a;
b = prblm.b;

I = speye(K+1);

% derivatives (on [-1,1], need to add a factor 2/(b-a) for each derivative
% when on [a,b] !!!)
MDK = derCheb_mat(K+1)*2/(b-a);
MD2K = MDK^2;

% vectors 1 and x in Chebyshev (on [-1,1])
e0 = zeros(K+1,1);
e0(1) = 1;
e1 = zeros(K+1,1);
e1(2) = 1/2;

% rescaled x and powers of x (on [a,b])
x = (b-a)/2*e1+(b+a)/2*e0;
Mx = convomat_coscos(x);
x2 = Mx*x;
Mx2 = Mx*Mx;
x3 = Mx*x2;

% Boundary conditions
Mat_BC =[zeros(4,K-3);speye(K-3)];
BC = [1, 2*ones(1,K);
       1, 2*(-1).^(1:K);
       (0:K).^2;
       (-1).^(0:K).*(0:K).^2];    
BC_to_inv = BC(:,1:4);
BC_remainder = BC(:,5:end);         
Mat_BC(1:4,:) = -BC_to_inv\BC_remainder;            
Trunc = [speye(K-3),zeros(K-3,4)];


%%
tab_I0 = 0*tab_alpha;
tab_pstar = 0*tab_alpha;
nb = length(tab_alpha);
coeff = 2^-5;
for j = nb:-1:1
    alpha = tab_alpha(j)
    pstar_alpha = trisection(pmin,pmax,alpha);
    tab_pstar(j) = pstar_alpha;
    tab_I0(j) = - Lambdap(pstar_alpha,alpha);
    if j>1
        pmin = pstar_alpha - 2*tol;
        deltapmax_approx = pstar_alpha*alpha/tab_alpha(j-1)^2 * (alpha-tab_alpha(j-1));
        pmax = pstar_alpha + coeff*deltapmax_approx + 2*tol;
    end
end   
    

function pstar = trisection(pmin,pmax,alpha)  
pmax_init = pmax;
Lpmin = Lambdap(pmin,alpha);
Lpmax = Lambdap(pmax,alpha);
it = 0;
while pmax-pmin>tol && it<it_max
    p1 = (2*pmin + pmax) / 3;
    p2 = (pmin + 2*pmax) / 3;
    Lp1 = Lambdap(p1,alpha);
    Lp2 = Lambdap(p2,alpha);
    if Lp1<Lp2
        if Lpmin<Lp1
            pmax = p1;
            Lpmax = Lp1;
        else
            pmax = p2;
            Lpmax = Lp2;
        end
    else
        if Lpmax<Lp2
            pmin = p2;
            Lpmin = Lp2;
        else
            pmin = p1;
            Lpmin = Lp1;
        end
    end
    it = it+1;
end
pstar = (pmin + pmax) / 2; %The approximate minimizer
if pmax_init-pstar<=tol
    coeff = 2*coeff;
    pmax = pstar_alpha + coeff*deltapmax_approx + 2*tol;
    pstar = trisection(pmin,pmax,alpha);
end
end

function Lp = Lambdap(p,alpha)
    temp = x3 - alpha*x;
    h = convomat_coscos(temp)*temp;
    Mh = 1/(2*sigma^2) * convomat_coscos(h);
    Hp = -sigma^2/2*MD2K + (1/2-p)*(alpha*I-3*Mx2) + Mh + shift*I;
    [~,Lambdap] = eigs(Trunc*Hp*Mat_BC,Trunc*Mat_BC,nb_eigs_SM,'SM');
    Lp = -min(real(diag(Lambdap))) + shift;
end

end

