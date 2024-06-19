function val = eval_FourCheb(u_coeffs,pts_Cheb,pts_Four,real_output)

% u of size (2N+1)x(K+1)

sz = size(u_coeffs);
N = (sz(1)-1)/2;
K = sz(2)-1;
if size(pts_Four,2) > 1
    pts_Four = transpose(pts_Four);
end
if size(pts_Cheb,1) > 1
    pts_Cheb = transpose(pts_Cheb);
end

Mat_Four = exp(1i*pts_Four*(-N:N));
if exist('intval','file') && isintval(pts_Four(1))
    Mat_Four = intersect(Mat_Four,cintval(0,1));
end
Mat_Cheb = cos( (0:K)' * acos(pts_Cheb) );
Mat_Cheb(2:K+1,:) = 2*Mat_Cheb(2:K+1,:);

val = Mat_Four*u_coeffs*Mat_Cheb;

if nargin > 2 && real_output
    val = real(val);
end