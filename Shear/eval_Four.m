function val = eval_Four(u_coeffs,pts,real_output)

% v of size (2N+1)x(K+1)

N = (size(u_coeffs,1)-1)/2;
if size(pts,2) > 1
    pts = transpose(pts);
end

Mat = exp(1i*pts*(-N:N));
if exist('intval','file') && isintval(pts(1))
    Mat = intersect(Mat,cintval(0,1));
end
val = Mat*u_coeffs;

if nargin > 2 && real_output
    val = real(val);
end