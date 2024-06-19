function out = prove_positivityFFT(u,rmin,N_FFT,K_FFT,N_FFT_max,K_FFT_max)

% Rigorously enclose the range of u, using interval arithmetic, Taylor
% expansions and the "interval FFT" trick, so as to check that u is indeed
% positive.

sz = size(u);
N = (sz(1)-1)/2;
K = sz(2)-1;

% Making sure that N_FFT and K_FFT are at least large enough to contain u
N_FFT_input = N_FFT;
K_FFT_input = K_FFT;
if N_FFT<2*N+1
    N_FFT = 2^ceil(log2(2*N+1));
end
if K_FFT<K
    K_FFT = 2^ceil(log2(K));
end
N_FFT_mid = N_FFT/2+1;

ipi = intval('pi');
delta_N = midrad(0,sup(ipi/(N_FFT)));
delta_K = midrad(0,sup(ipi/(2*K_FFT)));

% Widening corresponding to taking intervals in angle and parameter
delta_u = exp(1i*( (-N:N)'*delta_N + (0:K)*delta_K )); 

%% order 0
% Widened u
u_int = intval(zeros(N_FFT,K_FFT+1));
u_int(N_FFT_mid+(-N:N),1:K+1) = u.*delta_u;

% Validated lower bound for the approximate eigenfunctions
u_eval_direct = real(coeffs2grid_2d(u_int)); 
u_min = min(inf(u_eval_direct),[],'all');

if u_min > rmin
    out = true;
    fprintf("Rigorous lower bound for the eigenfunctions: %g\n",u_min-rmin)
    fprintf("Positivity of the eigenfunctions proven :)\n")
else
    %% order 2
    dudFour = derFour(u);
    dudCheb = derCheb(u);
    d2udFour2 = derFour(dudFour);
    d2udCheb2 = derCheb(dudCheb);
    d2udFourCheb = derCheb(dudFour);
        
    % Exact u
    u_mid = intval(zeros(N_FFT,K_FFT+1));
    u_mid(N_FFT_mid+(-N:N),1:K+1) = u;
    
    % Exact du
    dudFour_mid = intval(zeros(N_FFT,K_FFT+1));
    dudFour_mid(N_FFT_mid+(-N:N),1:K+1) = dudFour;
    dudCheb_mid = intval(zeros(N_FFT,K_FFT+1));
    dudCheb_mid(N_FFT_mid+(-N:N),1:K+1) = dudCheb;
    
    % Widened d2u
    d2udFour2_int = intval(zeros(N_FFT,K_FFT+1));
    d2udFour2_int(N_FFT_mid+(-N:N),1:K+1) = d2udFour2.*delta_u;
    d2udFourCheb_int = intval(zeros(N_FFT,K_FFT+1));
    d2udFourCheb_int(N_FFT_mid+(-N:N),1:K+1) = d2udFourCheb.*delta_u;
    d2udCheb2_int = intval(zeros(N_FFT,K_FFT+1));
    d2udCheb2_int(N_FFT_mid+(-N:N),1:K+1) = d2udCheb2.*delta_u;
    
    % Intervals
    dint_Four = delta_N*speye(N_FFT);
    x_grid = cos((K_FFT:-1:0)*ipi/K_FFT);
    radx_grid = [0, sup((x_grid(2:end)-x_grid(1:end-1))/2), 0]; %The two zeros account for the fact that we only need half intervals at the end points
    dint_Cheb = spdiags(infsup(-radx_grid(1:end-1),radx_grid(2:end))',0,K_FFT+1,K_FFT+1);
    
    % Validated lower bound for the approximate eigenfunctions
    u0_eval = real(coeffs2grid_2d(u_mid));
    dudFour_eval = dint_Four*real(coeffs2grid_2d(dudFour_mid));
    dudCheb_eval = real(coeffs2grid_2d(dudCheb_mid))*dint_Cheb;
    d2udFour2_eval = dint_Four.^2*real(coeffs2grid_2d(d2udFour2_int));
    d2udFourCheb_eval = dint_Four*real(coeffs2grid_2d(d2udFourCheb_int))*dint_Cheb;
    d2udCheb2_eval = real(coeffs2grid_2d(d2udCheb2_int))*dint_Cheb.^2;
    u_eval_taylor = u0_eval + dudFour_eval + dudCheb_eval + 1/2 * (d2udFour2_eval + 2*d2udFourCheb_eval + d2udCheb2_eval);

    u_min = max( u_min, min(inf(u_eval_taylor),[],'all') );

    if u_min > rmin
        out = true;
        fprintf("Rigorous lower bound for the eigenfunctions: %g\n",u_min-rmin)
        fprintf("Positivity of the eigenfunctions proven :)\n")
    else
        fprintf("Can only prove that the approximate eigenfunctions are above: %g\n",u_min-rmin)
        % u_eval = intersect(u_eval_direct,u_eval_taylor);
        % problematic_pts = ( inf(u_eval) <= rmin );
        % figure
        % spy(problematic_pts)
        % xlabel('$k$', 'Interpreter', 'latex')
        % ylabel('$\phi$', 'Interpreter', 'latex')
        % title('Locations where positivity could not be proven')
        % set(gca,'FontSize',15)
        % drawnow
        if N_FFT_input < N_FFT_max || K_FFT_input < K_FFT_max
            N_FFT = min(2*N_FFT_input,N_FFT_max);
            K_FFT = min(2*K_FFT_input,K_FFT_max);
            fprintf("Trying again to prove positivity with a finer discretization...\n")
            out = prove_positivityFFT(u,rmin,N_FFT,K_FFT,N_FFT_max,K_FFT_max);
        else
            out = false;
            fprintf("Unable to prove the positivity of the eigenfunctions, and therefore to guarantee that we have the correct eigenvalues.\n")
        end
    end
end

        