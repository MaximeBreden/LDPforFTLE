function v_coeffs = grid2coeffs(v_grid)

% Goes from values at the Chebyshev grid to Chebyshev coefficients, along
% the last dimension of v_grid. For instance, if v_grid is a 3d-tensor,
% each 1d vector v_grid(i,j,:) is interpreted as the values of a function
% at the different points x_k of the grid.
% Conventions: For Chebyshev coefficients v_0,...,v_K, the associated 
% function v is given by v(x) = v_0 + 2*sum_{k=1}^K v_k T_k(x). The
% Chebyshev grid is on [-1,1], and made of Chebsyhev points of the second
% kind: x_k = cos ((K-k)*pi/K), k=0,...,K.
% The computation is done using the inverse FFT.

sz = size(v_grid);
K = sz(end)-1;
n = length(sz);
C = repmat({':'},1,n-1); % used to do nothing on the n-1 first dimensions

v_temp = cat( n, v_grid(C{:},K+1:-1:2), v_grid(C{:},1:K) );
v_coeffs = my_ifft(v_temp,2*K,n);
if exist('intval','file') && isintval(v_grid(1))
    v_coeffs(C{:},2:K) = intersect( v_coeffs(C{:},2:K), v_coeffs(C{:},2*K:-1:K+2) );
else
    v_coeffs(C{:},2:K) = ( v_coeffs(C{:},2:K) + v_coeffs(C{:},2*K:-1:K+2) ) / 2;
end
v_coeffs = v_coeffs(C{:},1:K+1);
v_coeffs(C{:},K+1) = v_coeffs(C{:},K+1)/2;

if isreal(v_grid)
    v_coeffs = real(v_coeffs);
end

% % The 1d version
% K = length(v_grid)-1;
% v_temp = [ v_grid(K+1:-1:2); v_grid(1:K) ];
% v_coeffs = ifft(v_temp);
% v_coeffs = v_coeffs(1:K+1);
% % v_coeffs = [ v_coeffs(1); (v_coeffs(2:K)+v_coeffs(2*K:-1:K+2))/2; v_coeffs(K+1) ];