function val = eval_Cheb(u_coeffs,pts)

% v of size (2N+1)x(K+1)

K = size(u_coeffs,2)-1;
if size(pts,1) > 1
    pts = transpose(pts);
end

Mat = cos( (0:K)' * acos(pts) );
Mat(2:K+1,:) = 2*Mat(2:K+1,:);
val = u_coeffs*Mat;

