function v_grid = coeffs2grid(v_coeffs)

% Goes from Chebyshev coefficients to values at the Chebyshev grid, along
% the last dimension of v_coeffs. For instance, if v_coeffs is a 3d-tensor,
% each 1d vector v_coeffs(i,j,:) is interpreted as a sequence of Chebyshev
% coefficients.
% Conventions: For Chebyshev coefficients v_0,...,v_K, the associated 
% function v is given by v(x) = v_0 + 2*sum_{k=1}^K v_k T_k(x). The
% Chebyshev grid is on [-1,1], and made of Chebsyhev points of the second
% kind: x_k = cos ((K-k)*pi/K), k=0,...,K.
% The computation is done using the FFT.

sz = size(v_coeffs);
K = sz(end)-1;
n = length(sz);
C = repmat({':'},1,n-1); % used to do nothing on the n-1 first dimensions

v_temp = cat( n, v_coeffs, v_coeffs(C{:},K:-1:2) );
v_temp(C{:},K+1) = 2*v_temp(C{:},K+1);
v_grid = my_fft(v_temp,2*K,n);
if exist('intval','file') && isintval(v_coeffs(1))
    % v_grid(C{:},2:K) = intersect( v_grid(C{:},2:K), v_grid(C{:},2*K:-1:K+2) );
    % %The intersection can output NaN for very small inputs
else
    v_grid(C{:},2:K) = ( v_grid(C{:},2:K) + v_grid(C{:},2*K:-1:K+2) ) / 2;
end
v_grid = v_grid(C{:},K+1:-1:1);

if isreal(v_coeffs)
    v_grid = real(v_grid);
end

% % The 1d version
% K = length(v_coeffs)-1;
% v_temp = [ v_coeffs; v_coeffs(K:-1:2) ];
% v_grid = fft(v_temp);
% v_grid = v_grid(K+1:-1:1);
% % v_grid = [ v_grid(K+1); (v_grid(K+2:2*K)+v_grid(K:-1:2))/2; v_grid(1) ];