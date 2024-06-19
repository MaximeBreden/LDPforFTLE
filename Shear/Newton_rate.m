function q_coeffs = Newton_rate(q_coeffs,dlambda_coeffs,d2lambda_coeffs,problem,tol,it_max,show)

if nargin < 5
    show = false;
    if nargin < 4
        it_max = 20;
        if nargin < 3
            tol = 10^(-12);
        end
    end
end

if show
    fprintf('Newton iterations, residual errors:\n')
end

p_mid = problem.continuation.p_mid;
p_rad = problem.continuation.p_rad;

a_mid = problem.rate.a_mid;
a_rad = problem.rate.a_rad;

S = length(q_coeffs)-1;
a_coeffs = zeros(1,S+1);
a_coeffs(1) = a_mid;
a_coeffs(2) = 1/2*a_rad;
a_grid = coeffs2grid(a_coeffs);
q_grid = coeffs2grid(q_coeffs);

dlambda_q_grid = real( eval_Cheb( dlambda_coeffs, (q_grid-p_mid)/p_rad ) );
F_grid = dlambda_q_grid - a_grid;

err = max(abs(F_grid));
if show
    disp(err)
end
it = 0;
while (err > tol) && (it < it_max) && (err < 10^10)
    DF_grid = real( eval_Cheb( d2lambda_coeffs, (q_grid-p_mid)/p_rad ) );
    q_grid = q_grid - F_grid./DF_grid;
    dlambda_q_grid = real( eval_Cheb( dlambda_coeffs, (q_grid-p_mid)/p_rad ) );
    F_grid = dlambda_q_grid - a_grid;

    err = max(abs(F_grid));
    if show
        disp(err)
    end
    it = it+1;
end

q_coeffs = grid2coeffs(q_grid);

if show
    fprintf('\n')
end

if err > tol || isnan(err)
    warning('\nNewton method may not have converged, the residual error is %e',err)
end