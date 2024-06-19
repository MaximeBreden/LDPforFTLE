function v = derCheb(u)

% The output v contains the Chebyshev coefficients of u'. The convention is
% that u = u_0 + 2*\sum u_k T_k. 

K = size(u,2)-1;
if K == 0
    v = zeros(size(u,1),1);
else
    v = zeros(size(u));
    if exist('intval','file') && isintval(u(1))
        v = intval(v);
    end
    even_u = (2:2:K);
    odd_u = (1:2:K);
    even_v = (0:2:K-1);
    odd_v = (1:2:K-1);
    v(:,even_v+1) = 2*icumsum( u(:,odd_u+1)*diag(odd_u), 2, 'reverse' );
    v(:,odd_v+1) = 2*icumsum( u(:,even_u+1)*diag(even_u), 2, 'reverse' );
end
