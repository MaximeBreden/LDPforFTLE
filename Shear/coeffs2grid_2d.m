function v_grid = coeffs2grid_2d(v_coeffs)

% Goes from Fourier x Chebyshev coefficients to values at a grid (uniform
% in the angle, and Chebyshev in space)
% Conventions: For Chebyshev coefficients v_0,...,v_K, the associated 
% function v is given by v(x) = v_0 + 2*sum_{k=1}^K v_k T_k(x). The
% Chebyshev grid is on [-1,1], and made of Chebsyhev points of the second
% kind: x_k = cos ((K-k)*pi/K), k=0,...,K.
% The computation is done using the FFT.

sz = size(v_coeffs);
K = sz(2)-1;
v_F = flip(my_fft(ifftshift(v_coeffs,1)),1);

v_temp = [v_F, v_F(:,K:-1:2)];
v_temp(:,K+1) = 2*v_temp(:,K+1);
v_grid = my_fft(v_temp,2*K,2);
if exist('intval','file') && isintval(v_coeffs(1))
    % v_grid(:,2:K) = intersect( v_grid(:,2:K), v_grid(:,2*K:-1:K+2) );%    
    % %The intersection can output NaN for very small inputs
else
    v_grid(:,2:K) = ( v_grid(:,2:K) + v_grid(:,2*K:-1:K+2) ) / 2;
end
v_grid = v_grid(:,K+1:-1:1);

