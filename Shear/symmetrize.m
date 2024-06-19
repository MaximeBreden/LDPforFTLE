function Y = symmetrize(X)

N = (length(X)-2)/2;
u = X(1:2*N+1);
lambda = X(2*N+2);

Y = [ (u+flip(conj(u)))/2; real(lambda)];