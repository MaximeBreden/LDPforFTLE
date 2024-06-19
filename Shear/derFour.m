function v = derFour(u)

N = (size(u,1)-1)/2;
Mat_Der = diag(1i*(-N:N));
if exist('intval','file') && isintval(u(1))
    Mat_Der = intval(Mat_Der);
end
v = Mat_Der * u;