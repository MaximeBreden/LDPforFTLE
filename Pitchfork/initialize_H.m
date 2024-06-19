function operator_H = initialize_H(prblm,K)

% All the matrix representations of the operators involved in H

% Parameters of the SDE
alpha = prblm.alpha;
sigma = prblm.sigma;
p = prblm.p;
ialpha = prblm.ialpha;
isigma = prblm.isigma;
ip = prblm.ip;

% Domain
a = prblm.a;
b = prblm.b;
ia = prblm.ia;
ib = prblm.ib;

% Coefficients of the base problem
if isfield(prblm,'A')
    A = prblm.A;
    B = prblm.B;
    iA = prblm.iA;
    iB = prblm.iB;
end

iK = K+6;
I = speye(K+1);
iI = speye(iK+1);

% derivatives (on [-1,1], need to add a factor 2/(b-a) for each derivative
% when on [a,b] !!!)
MDK = derCheb_mat(K+1)*2/(b-a);
MD2K = MDK^2;
iMDK = derCheb_mat(iK+1)*2/(ib-ia);
iMD2K = iMDK^2;

% vectors 1 and x in Chebyshev (on [-1,1])
e0 = zeros(K+1,1);
e0(1) = 1;
e1 = zeros(K+1,1);
e1(2) = 1/2;
ie0 = zeros(iK+1,1);
ie0(1) = 1;
ie1 = zeros(iK+1,1);
ie1(2) = 1/2;

% rescaled x and powers of x (on [a,b])
x = (b-a)/2*e1+(b+a)/2*e0;
Mx = convomat_coscos(x);
x2 = Mx*x;
Mx2 = Mx*Mx;
x3 = Mx*x2;
temp = x3 - alpha*x;
h = convomat_coscos(temp)*temp;
Mh = 1/(2*sigma^2) * convomat_coscos(h);
ix = (ib-ia)/2*ie1+(ib+ia)/2*ie0;
iMx = convomat_coscos(ix);
ix2 = iMx*ix;
iMx2 = iMx*iMx;
ix3 = iMx*ix2;
itemp = ix3 - ialpha*ix;
ih = convomat_coscos(itemp)*itemp;
iMh = 1/(2*isigma^2) * convomat_coscos(ih);

% Shifting the spectrum to be safely away from zero
if isfield(prblm,'A')
    PEVal_H0 = 1/2*sigma*sqrt(2*A)+B;
    shift = max(-PEVal_H0,0) + 1;
    ishift = intval(shift);
else
    shift = 10;
    ishift = intval(10);
end

% the operators H and the associated base problem H^{(0)}
H = -sigma^2/2*MD2K + (1/2-p)*(alpha*I-3*Mx2) + Mh + shift*I;
iH = -isigma^2/2*iMD2K + (1/2-ip)*(ialpha*iI-3*iMx2) + iMh + ishift*iI;
if isfield(prblm,'A')
    H_0 = -sigma^2/2*MD2K + A*Mx2 + B*I + shift*I;
    iH_0 = -isigma^2/2*iMD2K + iA*iMx2 + iB*iI + ishift*iI;
end

% Matrices of the scalar product (by default, the scalar product is 
% 1/2*int_{-1}^1 f*g, need to multiply by b-a to recover the scalar product
% on [a,b] with weight 1.
Mat_scal = (b-a)*scalCheb_mat(K);
iMat_scal = (ib-ia)*iscalCheb_mat(iK);

% Boundary conditions
Mat_BC =[zeros(4,K-3);speye(K-3)];
BC = [1, 2*ones(1,K);
       1, 2*(-1).^(1:K);
       (0:K).^2;
       (-1).^(0:K).*(0:K).^2];    
BC_to_inv = BC(:,1:4);
BC_remainder = BC(:,5:end);         
Mat_BC(1:4,:) = -BC_to_inv\BC_remainder;
iMat_BC =[zeros(4,K-3);speye(K-3)];
iBC = [1, 2*ones(1,K);
        1, 2*(-1).^(1:K);
        (0:K).^2;
        (-1).^(0:K).*(0:K).^2];    
iBC_to_inv = iBC(:,1:4);
iBC_remainder = iBC(:,5:end);         
iMat_BC(1:4,:) = -iBC_to_inv\iBC_remainder;

Trunc = [speye(K-3),zeros(K-3,4)];
Reset_for_BC = [zeros(K-3,4),speye(K-3)];
Pad = [speye(K+1);zeros(6,K+1)];


%% Structure for H

operator_H.Id = I;
operator_H.H = H;
if isfield(prblm,'A')
    operator_H.H0 = H_0;
end
operator_H.Mat_scal = Mat_scal;
operator_H.Mat_BC = Mat_BC;
operator_H.Reset_for_BC = Reset_for_BC;
operator_H.Trunc = Trunc;
operator_H.Pad = Pad;
operator_H.shift = shift;

operator_H.iH = iH;
if isfield(prblm,'A')
    operator_H.iH0 = iH_0;
end
operator_H.iMat_scal = iMat_scal;
operator_H.iMat_BC = iMat_BC;
operator_H.ishift = ishift;



