function [alpha,b,sigma,p,p_mid,p_rad,p_grid] = obtain_parameters(problem)

alpha = problem.alpha;
b = problem.b;
sigma = problem.sigma;

if isfield(problem,'p')
    p = problem.p;
else
    p = [];
end

if isfield(problem.continuation,'p_mid')
    p_mid = problem.continuation.p_mid;
else
    p_mid = [];
end

if isfield(problem.continuation,'p_rad')
    p_rad = problem.continuation.p_rad;
else
    p_rad = [];
end

if isfield(problem.continuation,'p_grid')
    p_grid = problem.continuation.p_grid;
else
    p_grid = [];
end