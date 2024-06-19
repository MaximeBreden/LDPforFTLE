function Mat=iscalCheb_mat(K)

%Computation of the matrix (bilinear form) corresponding to the L^2 scalar 
%product in the Cheybyshev basis : <f,g> = 1/2 * int_{-1}^1 f(x)*g(x) dx

k = 0:K;
Mat = -2*(1./(intval(bsxfun(@plus,k',k)).^2-1)+1./(intval(bsxfun(@minus,k',k)).^2-1));
Mat(1:2:K+1,2:2:K+1) = 0;
Mat(2:2:K+1,1:2:K+1) = 0;
Mat(1,:) = Mat(1,:)/2;
Mat(:,1) = Mat(:,1)/2;


