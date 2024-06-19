function X = Newton_eigenproblem(X,problem,tol,it_max,show)

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
FX = F_eigenproblem(X,problem);   
err = norm(FX,1);
if show
    disp(err)
end
it = 0;
while (err > tol) && (it < it_max) && (err < 10^10)
    X = X - DF_eigenproblem(X,problem)\FX;
    X = symmetrize(X);
    FX = F_eigenproblem(X,problem);
    err = norm(FX,1);
    if show
        disp(err)
    end
    it = it+1;
end

if show
    fprintf('\n')
end

if err > tol || isnan(err)
    warning('\nNewton method may not have converged, the residual error is %e',err)
end