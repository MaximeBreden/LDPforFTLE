function pstar = find_argmin_Lambda(prblm,K,pmin,pmax,opt,tol,it_max)

% Non-rigorous computation of the minimizer p of Lambda, in [pmin,pmax].
% The field prblm contains the paramters alpha and sigma of the SDE, as
% well as the spatial domain [a,b] used to approximate the problem. This
% approximation is based on truncated Chebyshev series of order K. 

% If the optional parameter opt is > 0, a plot of Lambda using a 
% subdivision of [pmin,pmax] of size opt will be produced.
% tol and it_max are optional parameters using for the (basic)
% minimum-finding algorithm.

% The nested function Lambdap uses a discretized version of H_p together
% with Matlab's eigenvalue solver in order to approximate Lambda(p).

global nb_eigs_SM

if nargin < 7
    it_max = 1e3;
    if nargin < 6
        tol = 1e-6;
        if nargin < 5
            opt = 0;
        end
    end
end

%% Init

% Parameters of the SDE
alpha = prblm.alpha;
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
temp = x3 - alpha*x;
h = convomat_coscos(temp)*temp;
Mh = 1/(2*sigma^2) * convomat_coscos(h);

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

%% Non rigorous plot of Lambda (optional)
if opt>0
    tab_p = linspace(pmin,pmax,opt);
    tab_Lambda = 0*tab_p;
    for j = 1:opt
        tab_Lambda(j) = Lambdap(tab_p(j));
    end
    figure
    plot(tab_p,tab_Lambda,'k','Linewidth',2)
    xlabel('p')
    ylabel('$\Lambda(p)$','Interpreter','Latex')
    title('Non rigorous computation of the moment Lyapunov function','Interpreter','Latex')
    set(gca,'FontSize',15) 
    axis tight
    drawnow
    [~,ind] = sort(tab_Lambda);
    ind = ind(1);
    if ind == 1
        pmax = tab_p(2);
        % fprintf("\nThe minimum of Lambda may be to the left of the selected interval [pmin,pmax]\n")
    elseif ind == opt
        pmin = tab_p(end-1);
        % fprintf("\nThe minimum of Lambda may be to the right of the selected interval [pmin,pmax]\n")
    else
        pmin = tab_p(ind-1);
        pmax = tab_p(ind+1);
    end
end

%% Trisection to find the minimum

Lpmin = Lambdap(pmin);
Lpmax = Lambdap(pmax);

it = 0;
while pmax-pmin>tol && it<it_max
    p1 = (2*pmin + pmax) / 3;
    p2 = (pmin + 2*pmax) / 3;
    Lp1 = Lambdap(p1);
    Lp2 = Lambdap(p2);
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

if pmax-pmin>tol
    fprintf("\nThe minimization may not have converged: pmin = %g, pmax = %g\n",pmin,pmax)
end
pstar = (pmin + pmax) / 2; %The approximate minimizer

Hpstar = -sigma^2/2*MD2K + (1/2-pstar)*(alpha*I-3*Mx2) + Mh;
[EVect,EVal] = eigs(Trunc*Hpstar*Mat_BC,Trunc*Mat_BC,nb_eigs_SM,'SM');
[~,ind] = sort(real(diag(EVal)));
PEVect = Mat_BC*EVect(:,ind(1));
if PEVect(1)<0
    PEVect = -PEVect;
end

if opt>0
    figure
    plotCheb(PEVect,prblm.a,prblm.b,1e3);
    xlabel('x')
    title("Approximate principal eigenfunction for the approximate minimizer p*")
    drawnow
end

function Lp = Lambdap(p)
    Hp = -sigma^2/2*MD2K + (1/2-p)*(alpha*I-3*Mx2) + Mh;
    [~,Lambdap] = eigs(Trunc*Hp*Mat_BC,Trunc*Mat_BC,nb_eigs_SM,'SM');
    Lp = -min(real(diag(Lambdap)));
end

end