function q_grid = approx_inverse(dlambda_coeffs,S,problem,nb)

cheb_grid_S = cos( ((S:-1:0)*pi)/S );
a_grid_S = problem.rate.a_mid + cheb_grid_S*problem.rate.a_rad;
pts = linspace(-1,1,nb)';
dlambda_val = eval_Cheb(dlambda_coeffs,pts);
if dlambda_val(1) >= a_grid_S(1)
    error("This should not happen, make sure the values of a_mid and a_rad chosen are compatible with the range of \Lambda' on the selected k-interval")
end

q_grid = zeros(1,S+1);
ind = 1;
pts = problem.continuation.p_mid + pts*problem.continuation.p_rad;
for s = 1:S+1
    while dlambda_val(ind) < a_grid_S(s) && ind < nb+1
        ind = ind+1;
    end
    if ind == nb+1
        error("This should not happen, make sure the values of a_mid and a_rad chosen are compatible with the range of \Lambda' on the selected k-interval")
    end
    l = (dlambda_val(ind)-a_grid_S(s)) / (dlambda_val(ind)-dlambda_val(ind-1));
    q_grid(s) = (1-l)*pts(ind) + l*pts(ind-1);
end
