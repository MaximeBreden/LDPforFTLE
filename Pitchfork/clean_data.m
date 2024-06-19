function v = clean_data(v,tol)

%Putting to 0 the higher order Chebyshev modes which seem to be too
%large due to numerical errors.

if nargin < 2
    tol = 1e-12;
end

vv = sum(abs(v),2);
indices = 1:length(vv);
mask = min( vv > tol, indices' >= length(vv)/2 );
[~,ind_min] = min(vv(mask));
indices_mask = indices(mask);
v(indices_mask(ind_min)+1:end,:) = 0;