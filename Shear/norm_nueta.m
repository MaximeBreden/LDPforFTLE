function val = norm_nueta(X,nu,eta)

sz = size(X);
l = length(sz);
if l == 2
    N = (sz(1)-2)/2;
    K = sz(2)-1;
    w_N = [nu.^abs(-N:N), 1];
    w_K = [1 2*eta.^(1:K)];
    val = w_N * abs(X) * w_K';
elseif l == 3
    Nr = (sz(1)-2)/2;
    Nc = (sz(2)-2)/2;
    K = sz(3)-1;
    w_Nr = [nu.^abs(-Nr:Nr), 1];
    w_Nc = [nu.^abs(-Nc:Nc), 1];
    w_K = [1 2*eta.^(1:K)];
    w_K = repmat( reshape(w_K,[1,1,K+1]), [2*Nr+2,2*Nc+2,1] );
    XK = sum( abs(X).*w_K, 3);
    val = max( (w_Nr * XK) ./ w_Nc );
end