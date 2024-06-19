function plot_sols(X_grid,problem,nb_Cheb,nb_Four,nb_Cheb2)

if nargin<5
    nb_Cheb2 = 1e3;
end
[~,~,~,~,p_mid,p_rad] = obtain_parameters(problem);

u_grid = X_grid(1:end-1,:);
u_coeffs = grid2coeffs(u_grid);
lambda_grid = X_grid(end,:);
lambda_coeffs = grid2coeffs(lambda_grid);

figure
pts_Cheb = linspace(-1,1,nb_Cheb);
pts_Four = linspace(0,2*pi,nb_Four);
u_eval = eval_Four( eval_Cheb(u_coeffs,pts_Cheb), pts_Four, true);
[X_Four,X_Cheb] = meshgrid(p_mid+p_rad*pts_Cheb,pts_Four);
surf(X_Four,X_Cheb,u_eval)
xlabel('$p$', 'Interpreter', 'latex')
ylabel('$\phi$', 'Interpreter', 'latex')
title('Eigenfunctions')
set(gca,'FontSize',15)
axis tight

figure
pts_Cheb = linspace(-1,1,nb_Cheb2);
lambda_eval = eval_Cheb(lambda_coeffs,pts_Cheb);
plot(p_mid+pts_Cheb*p_rad,lambda_eval, 'k', 'linewidth', 2)
xlabel('$p$', 'Interpreter', 'latex')
ylabel('$\Lambda(p)$', 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight
drawnow