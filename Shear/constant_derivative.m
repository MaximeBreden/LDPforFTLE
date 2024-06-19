function val = constant_derivative(eta,nb,a,b)

% The constant C = C(nb,eta,a,b) allowing to control derivatives from  the
% weighted l^1 norms: sup_{[a,b]} |u^{(nb)}| <= C * || u ||_{\ell^1_\eta}.

if not( eta>1 )
    val = Inf;
    return
end

if a<-1 && b>1
    val = Inf;
    fprintf("This case in not implemented")
    return
else
    x = max(abs([a,b]));
    if exist('intval','file') && isintval(eta)
        x2 = intval(x)^2;
    else
        x2 = x^2;
    end
end

if x > 2/(eta+1/eta)
    dist = (eta+1/eta)/2 - x;
else
    dist = (eta-1/eta)/2 * sqrt(1-x2);
end
val1 = factorial(nb) / dist^nb;

if nb == 1
    nlim = 2/log(eta);
    n = floor(i2f(nlim,'inf')) : ceil(i2f(nlim,'sup'));
    val2 = max(n.^2./eta.^n);
elseif nb == 2
    nlim = 4/log(eta);
    n = 2 : max(ceil(i2f(nlim,'sup')),2);
    val2 = max(n.^2.*(n.^2-1)./eta.^n) / 3;
else
    %Estimate based on the Chebyshev series not implemented for nb > 2
    val2 = Inf;
end

val = min(val1,val2);

    
