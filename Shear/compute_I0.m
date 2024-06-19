function [I0,iI0] = compute_I0(lambda_coeffs,problem,irmin)

% Computing I(0) = - inf Lambda. A rigorous enclosure of I(0) is provided
% if a second output is asked for.

p_mid = i2f(problem.continuation.p_mid);
p_rad = i2f(problem.continuation.p_rad);
dlambda_coeffs = 1/p_rad * derCheb(lambda_coeffs);
d2lambda_coeffs = 1/p_rad * derCheb(dlambda_coeffs);


%%% Finding the approximate zero of Lambda'
y = p_mid;
it_max = 1e2;
tol = 1e-12;
it = 0;
dlambda_y = eval_Cheb(dlambda_coeffs,(y-p_mid)/p_rad);
while abs(dlambda_y)>tol && it<it_max
    y = y - dlambda_y/eval_Cheb(d2lambda_coeffs,(y-p_mid)/p_rad);
    dlambda_y = eval_Cheb(dlambda_coeffs,(y-p_mid)/p_rad);
    it = it+1;
end
if abs(dlambda_y)>tol
    fprintf("\nNewton's method may not have converged, final residual error: %g\n",abs(dlambda_y))
end
if abs(dlambda_y)>1e-5
    fprintf("No zero of Lambda' found, unable to compute I(0)\n")
    I0 = NaN;
    if nargout>1
        iI0 = NaN;
    end
else
    I0 = - eval_Cheb(lambda_coeffs,(y-p_mid)/p_rad);
    fprintf("\nApproximate value of I(0): %g\n",I0)  

    if nargout > 1
        %%% Getting a guaranteed enclosure of the zero y of Lambda', and then of
        %%% Lambda'(y)
        ilambda_coeffs = intval(lambda_coeffs);
        idlambda_coeffs = derCheb(ilambda_coeffs) / p_rad;
        ieta = problem.proof.eta;
        y_resc = (y-p_mid)/p_rad;
        if abs(y_resc) >= 1 
             fprintf("The obtained approximate zero of Lambda' is not within the values of k on which Lambda was validated.\n")
             fprintf("We are therefore unable to rigorously enclose I(0) :(\n")
             iI0 = NaN;
        else
            %%% Finding a suitable lower bound for y (by finding ymin such
            %%% that Lambda'(ymin)is guaranteed to be <= 0)
            radius = 1e-16;
            y_resc_min = max(y_resc-radius,-1);
            while sup(eval_Cheb(idlambda_coeffs,y_resc_min)) > -constant_derivative(ieta,1,y_resc_min,y_resc_min)*irmin && radius <1
                radius = 2*radius;
                y_resc_min = max(y_resc-radius,-1);
            end
            if sup(eval_Cheb(idlambda_coeffs,y_resc_min)) > -constant_derivative(ieta,1,y_resc_min,y_resc_min)*irmin
                fprintf("Unable to enclose the 0 of Lambda', and therefore to rigorously enclose I(0) :(\n")
                iI0 = NaN;
                return
            end
    
            %%% Finding a suitable upper bound for y (by finding ymax such
            %%% that Lambda'(ymax)is guaranteed to be >= 0)
            radius = 1e-16;
            y_resc_max = min(y_resc+radius,1);
            while inf(eval_Cheb(idlambda_coeffs,y_resc_max)) < constant_derivative(ieta,1,y_resc_max,y_resc_max)*irmin && radius <1
                radius = 2*radius;
                y_resc_max = min(y_resc+radius,1);
            end
            if inf(eval_Cheb(idlambda_coeffs,y_resc_max)) < constant_derivative(ieta,1,y_resc_max,y_resc_max)*irmin
                fprintf("Unable to enclose the zero of Lambda', and therefore to rigorously enclose I(0) :(\n")
                iI0 = NaN;
            else
                % [y_resc_min,y_resc_max] contains the (rescaled) zero of Lambda'
                fprintf("I(0) rigorously enclosed:\n")
                iI0 = - eval_Cheb(ilambda_coeffs,infsup(y_resc_min,y_resc_max)) + midrad(0,irmin);
                infsup(iI0)
            end
        end
    end
end