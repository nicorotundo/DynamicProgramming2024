%{
--------------------------------------------------------------------------------------------------------------------------
Purpose: Coding part of problem set 5
Created: Nico Rotundo 2024-12
--------------------------------------------------------------------------------------------------------------------------
%}

%% Define parameters and initialize grids

%{ 
--------------------------------------------------------------------------------------------------------------------------
Define parameters
--------------------------------------------------------------------------------------------------------------------------
%} 

clear all;

% Relative risk aversion coefficient 
sigma = 2;

% Productivity
z_u = 1;
z_e = 2;
z = [z_u, z_e];

% Transition rates 
lambda_u = 1/3;
lambda_e = 1/3;
lambda = [lambda_u, lambda_e];

% Discount rate 
rho = 0.05;

% Capital depreciation rate
delta = 0.05;

% Capital share in production
alpha = 1/3; 

% TFP
A = .1; 

%{ 
--------------------------------------------------------------------------------------------------------------------------
Define capital grid
--------------------------------------------------------------------------------------------------------------------------
%} 

% Borrowing constraint minimum
k_min = 0; 

% Borrowing constraint maximum
k_max = 20; 

% Number of grid points
num_points = 1000; 

k_grid = linspace(k_min, k_max, num_points)';
dk = (k_max-k_min)/(num_points-1);

kk = [k_grid, k_grid]; % num_points*2 matrix

%{ 
--------------------------------------------------------------------------------------------------------------------------
Define utility function and its derivative 
--------------------------------------------------------------------------------------------------------------------------
%}

% Utility function (CRRA), handling vector inputs
U = @(c) (c.^(1 - sigma)) / (1 - sigma); 

% Derivative of utility, handling vector inputs
U_prime = @(c) c.^(-sigma);

% Inverse of derivative of utility 
U_prime_inv = @(Vp) (Vp).^(-1 / sigma); 

%{ 
--------------------------------------------------------------------------------------------------------------------------
Firm production technology and FOCs
--------------------------------------------------------------------------------------------------------------------------
%}

% Production technology Y = A K^alpha L^{1-alpha}
Y = @(K,L) K.^alpha*L.^(1-alpha);

% F.O.C.:
w_foc = @(K,L) (1-alpha)*A*(K/L).^alpha;
r_foc = @(K,L) alpha*A*(L/K).^(1-alpha); 
K_partial = @(L,r) L*((alpha*A/(r+delta)).^(1/(1-alpha))) ;

%{ 
--------------------------------------------------------------------------------------------------------------------------
Tuning parameters (general)
--------------------------------------------------------------------------------------------------------------------------
%} 

% Step size: can be arbitrarily large in implicit method
Delta = 1000;

% The maximum number of value function iterations
max_iterations_vf = 100;

% Tolerance for value function iterations
tolerance = 10^(-6);

 %{ 
--------------------------------------------------------------------------------------------------------------------------
Tuning parameters (interest rate iteration)
--------------------------------------------------------------------------------------------------------------------------
%} 

% The maximum number of interest rate iterations
num_iterations_r = 1000;

% Tolerance for interest rate iterations
tolerance_S = 10^(-5);

% Parameters for MIT shock
T = 100;
It = 1000;
nu = 0.2; 

%% Initialize matrices before iterating

%{ 
--------------------------------------------------------------------------------------------------------------------------
Define differential operator
--------------------------------------------------------------------------------------------------------------------------
%} 


% Forward; satisfies Df*V=dVf
Df = zeros(num_points, num_points);
for i = 1:num_points-1
    Df(i,i) = -1/dk; Df(i,i+1) = 1/dk;
end
Df = sparse(Df);

% Backward; satisfies Db*V=dV
Db = zeros(num_points, num_points);
for i = 2:num_points
    Db(i,i-1) = -1/dk; Db(i,i) = 1/dk;
end
Db = sparse(Db);

%{ 
--------------------------------------------------------------------------------------------------------------------------
Define A-switch matrix 
--------------------------------------------------------------------------------------------------------------------------
%} 

A_switch = [speye(num_points).*(-lambda(1)), speye(num_points).*lambda(1);
            speye(num_points).*lambda(2), speye(num_points).*(-lambda(2))];


%{ 
--------------------------------------------------------------------------------------------------------------------------
Calculate labor, capital, and wage levels 
--------------------------------------------------------------------------------------------------------------------------
%} 
    
% Labor
L = (z_e*lambda_u + z_u*lambda_e)/(lambda_e+lambda_u);

%{ 
--------------------------------------------------------------------------------------------------------------------------
Define initial guesses
--------------------------------------------------------------------------------------------------------------------------
%} 

% Set bounds on interest rate
r_min = 0.01;
r_max = 0.04;

% Define the number of points * 2 matrix
z = ones(num_points, 1).*z; 

% Initialize guess for, capital demand, and wage
r = (r_min+r_max)/2;
K = K_partial(L,r);
w = w_foc(K,L);

% Initial guess for value function
V_0 = U(w.*z + r.*kk)./rho; 
V = V_0;

%% Value function iteration
for nr=1:num_iterations_r
        
    r_r(nr) = r;
    
    % Use the value function solution from the previous interest rate iteration 
    % as the initial guess for the next iteration
    if nr>1
        v0 = V_r(:,:,nr-1);
        V = v0;
    end

    % Compute capital demand and wage for the equlibrium
    K = K_partial(L,r);
    w = w_foc(K,L);
    
    % Loop over number of interations for the value function
    for n=1:max_iterations_vf

        % Derivative of the forward value function 
        dVf = Df*V;

        % Derivative of the backward value function 
        dVb = Db*V;

        % Boundary condition on backwards value function i.e., the borrowing constraint; a>=a_min
        dVb(1,:) = U_prime(w.*z(1,:) + r.*kk(1,:)); 

        % Boundary condition on forward value function; a<=a_max 
        dVf(end,:) = U_prime(w.*z(end,:) + r.*kk(end,:)); 

        % Indicator whether value function is concave; for stability purposes
        I_concave = dVb > dVf; 

        % Compute optimal consumption using forward derivative
        cf = U_prime_inv(dVf);

        % Compute optimal consumption using backward derivative
        cb = U_prime_inv(dVb);
        
        % Compute optimal savings using forward derivative
        sf = w.*z + r.*kk - cf;

        % Compute optimal savings using backward derivative
        sb = w.*z + r.*kk - cb;

        % Upwind scheme
        If = sf>0;
        Ib = sb<0;
        I0 = 1-If-Ib;
        dV0 = U_prime(w.*z + r.*kk); % If sf<=0<=sb, set s=0

        dV_upwind = If.*dVf + Ib.*dVb + I0.*dV0;

        c = U_prime_inv(dV_upwind);

        % Update value function
        V_stacked = V(:);

        % Update consumption function
        c_stacked = c(:); 

        % A = SD
        SD_u = spdiags(If(:,1).*sf(:,1), 0, num_points, num_points)*Df + spdiags(Ib(:,1).*sb(:,1), 0, num_points, num_points)*Db; 
        SD_e = spdiags(If(:,2).*sf(:,2), 0, num_points, num_points)*Df + spdiags(Ib(:,2).*sb(:,2), 0, num_points, num_points)*Db; 
        SD = [SD_u, sparse(num_points, num_points);
            sparse(num_points, num_points), SD_e]; 

        % P = A + A_switch
        P = SD + A_switch;

        % B = [(rho + 1/Delta)*num_points - P]
        B = (rho + 1/Delta)*speye(2*num_points) - P; 

        % b = U(c) + 1/Delta*V
        b = U(c_stacked) + (1/Delta)*V_stacked;

        % V = B\b;
        V_update = B\b; 
        V_change = V_update - V_stacked;
        V = reshape(V_update, num_points, 2); 

        % Convergence criterion
        dist(n) = max(abs(V_change));
        
        if dist(n)<tolerance
            disp('Value function converged. Iteration = ')
            disp(n)
        break

        end
    end

    %% KF Equation

    % Solve for 0=gdot=P'*g
    PT = P';
    gdot_stacked = zeros(2*num_points,1);

    % Fix one value to obtain a non-singular matrix, otherwise matrix is singular
    i_fix = 1;
    gdot_stacked(i_fix)=.1;

    row_fix = [zeros(1,i_fix-1),1,zeros(1,2*num_points-i_fix)];
    AT(i_fix,:) = row_fix;

    g_stacked = PT\gdot_stacked; 

    % Normalization
    g_sum = g_stacked'*ones(2*num_points,1)*dk;
    g_stacked = g_stacked./g_sum;

    % Reshape
    gg = reshape(g_stacked, num_points, 2);

    % Compute interest-rate dependent variables for this iteration
    g_r(:,:,nr) = gg;
    kdot(:,:,nr) = w.*z + r.*kk - c;
    V_r(:,:,nr) = V;
    dV_r(:,:,nr) = dV_upwind;
    c_r(:,:,nr) = c;
    
    S(nr) = gg(:,1)'*k_grid*dk + gg(:,2)'*k_grid*dk;

    %% Update interest rate using Newton's method

    % Difference between savings and capital
    S_excess(nr) = K - S(nr);

    % Update interest rate only if iteration is past one 
    if nr > 1

        % Change in savings
        dS_excess = S_excess(nr) - S_excess(nr - 1);
        
        % Change in interest rate
        dr = r - r_old;

        % Derivative of savings with respect to interest rate
        S_excess_prime = dS_excess / dr;

        % Save current iteration interest rate
        r_old = r;
        
        % update interest rate
        r = r - S_excess(nr) / S_excess_prime;

    else
        r_old = r;

        r = r + 1e-4; 
    end
        
    % Check convergence 
    if sqrt(sum((r-r_old).^2))<tolerance
        disp(r)
        break 
    end

end

%% Plots 

%{ 
--------------------------------------------------------------------------------------------------------------------------
Question 6: Optimal consumption
--------------------------------------------------------------------------------------------------------------------------
%} 

% Initialize plot
set(gca, 'FontSize', 18)
    
    % Employed
    plot(k_grid, c_r(:,2,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', [41/255, 182/255, 164/255])
    hold on

    % Unmployed
    plot(k_grid, c_r(:,1,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', [250/255, 165/255, 35/255])
    
    hold off

    grid
        xlabel('Capital, k','FontSize', 14)
        
        set(gca, 'TickDir', 'out'); 
            box off; 
        
        ylabel('Consumption, c_j(k)','FontSize', 14)
        
        xlim([k_min k_max])

        legend(sprintf('Employed'), ...
            sprintf('Unemployed'), ...
            sprintf('r=%.4f', r), 'Location', 'best', 'FontSize', 14)
        
% Export graph
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_5/output/question_6_consumption.pdf', 'ContentType', 'vector');

%{ 
--------------------------------------------------------------------------------------------------------------------------
Question 6: Optimal savings
--------------------------------------------------------------------------------------------------------------------------
%} 

% Initialize plot
set(gca, 'FontSize', 18)
    
    % Employed
    plot(k_grid, kdot(:,2,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', [41/255, 182/255, 164/255])
    hold on

    % Unemployed
    plot(k_grid, kdot(:,1,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', [250/255, 165/255, 35/255])

    hold off

    grid
        xlabel('Capital, k', 'FontSize', 14)
        
        ylabel('Saving, s_j(k)', 'FontSize', 14)
        
        set(gca, 'TickDir', 'out'); 
            box off; 
        
        xlim([k_min k_max])
        
        legend(sprintf('Employed'), ...
            sprintf('Unemployed'), ...
            sprintf('r=%.4f', r), 'Location', 'best', 'FontSize', 14)

% Export graph
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_5/output/question_6_savings.pdf', 'ContentType', 'vector');

%{ 
--------------------------------------------------------------------------------------------------------------------------
Question 6: Stationary distribution
--------------------------------------------------------------------------------------------------------------------------
%} 

% Initialize plot
set(gca, 'FontSize', 18)
    
    % Employed
    plot(k_grid, g_r(:,2,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', [41/255, 182/255, 164/255])
    hold on

    % Unemployed
    plot(k_grid, g_r(:,1,nr), 'LineWidth', 2, 'LineStyle', '-', 'Color', [250/255, 165/255, 35/255])

    hold off

    grid
        xlabel('Capital, k', 'FontSize', 14)
        
        ylabel('Densities, g_j(k)', 'FontSize', 14)
        
        set(gca, 'TickDir', 'out'); 
            box off; 
        
        xlim([-.01 1])
        
        legend(sprintf('Employed'), ...
            sprintf('Unemployed'), ...
            sprintf('r=%.4f', r), 'Location', 'best', 'FontSize', 14)

% Export graph
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_5/output/question_6.pdf', 'ContentType', 'vector');

%{ 
--------------------------------------------------------------------------------------------------------------------------
Transition dynamics
--------------------------------------------------------------------------------------------------------------------------
%} 

% Initialize time grid
tt = linspace(0,T,It);

% Time step
dt = (tt(end)-tt(1))/(It-1);

% Initialize TFP grid
A_t = ones(It,1);
A_t(1) = A*0.97;

% Initialize arrays for general equilibrium
for it = 2:It
    dA_t = nu*(A-A_t(it-1))*dt;
    A_t(it) = A_t(it-1) + dA_t;
end

% Store values of general equilibrium
V_T = V;
r_T = r;
w_T = w;
K_T = K;
gg_T = gg;

% Set terminal condition for value function
V_t = zeros(num_points, 2, It);
V_t(:,:,It) = V_T;

% Initialize stationary distribution
g_stacked = zeros(2*num_points, It);
g_0 = gg(:);
g_stacked(:,1) = g_0;
gg_t = zeros(num_points, 2, It);
gg_t(:,:,1) = gg;

% Initialize array for capital supply
K_supply = zeros(It,1);

% Initial guesses for parameters 
L_t = L.*ones(It,1);
K_t = K_partial(L,r_T).*ones(It,1);
r_t = alpha.*A_t.*(L_t./K_t).^(1-alpha)-delta.*ones(It,1);
w_t = (1-alpha).*A_t.*(K_t./L_t).^alpha;

% Relaxation parameter for update in K
relax = 0.15; 

% Loop over iterations for interest rate 
for n=1:num_iterations_r

    % For given r(t), solve the HJB equation backwards with terminal condition
    for t=It:-1:1
        
        % Set value function at time t
        V = V_t(:,:,t);

        % Derivative of the forward value function 
        dVf = Df*V;

        % Derivative of the backward value function 
        dVb = Db*V;

        % Boundary condition on backwards value function i.e., the borrowing constraint; a>=a_min
        dVb(1,:) = U_prime(w_t(t).*z(1,:) + r_t(t).*kk(1,:)); 

        % Boundary condition on forward value function; a<=a_max
        dVf(end,:) = U_prime(w_t(t).*z(end,:) + r_t(t).*kk(end,:)); % a<=a_max is enforced which helps stability of the algorithm

        I_concave = dVb > dVf; % indicator whether value function is concave (problems arise if this is not the case)

        % Compute optimal consumption using forward derivative
        cf = U_prime_inv(dVf);

        % Compute optimal consumption using backward derivative
        cb = U_prime_inv(dVb);
        
        % Compute optimal savings using forward derivative
        sf = w_t(t).*z + r_t(t).*kk - cf;

        % Compute optimal savings using backward derivative
        sb = w_t(t).*z + r_t(t).*kk - cb;

        % Upwind scheme
        If = sf>0;
        Ib = sb<0;
        I0 = 1-If-Ib;
        dV0 = U_prime(w_t(t).*z + r_t(t).*kk); % If sf<=0<=sb, set s=0

        dV_upwind = If.*dVf + Ib.*dVb + I0.*dV0;

        c = U_prime_inv(dV_upwind);

        % Update value function
        V_stacked = V(:);

        % Update consumption function
        c_stacked = c(:); 
        
        % A = SD
        SD_u = spdiags(If(:,1).*sf(:,1), 0, num_points, num_points)*Df + spdiags(Ib(:,1).*sb(:,1), 0, num_points, num_points)*Db; 
        SD_e = spdiags(If(:,2).*sf(:,2), 0, num_points, num_points)*Df + spdiags(Ib(:,2).*sb(:,2), 0, num_points, num_points)*Db; 
        SD = [SD_u, sparse(num_points, num_points);
            sparse(num_points, num_points), SD_e]; 

        % P = A + A_switch
        P = SD + A_switch;
        P_t{t} = P;

        % B = [(rho + 1/Delta)*num_points - P]
        B = (rho + 1/dt)*speye(2*num_points) - P; 

        % b = U(c) + 1/Delta*V
        b = U(c_stacked) + (1/dt)*V_stacked;

        % V = B\b;
        V_update = B\b; 

        if t>1
            V_t(:,:,t-1) = reshape(V_update, num_points, 2); % num_points*2 matrix
        end

    end

    for t=1:It
        PT = P_t{t}';
        g_stacked(:,t+1) = (speye(2*num_points) - dt*PT)\g_stacked(:,t);
        gg_t(:,:,t+1) = reshape(g_stacked(:, t+1), num_points, 2);
        K_supply(t) = g_stacked(:,t)'*kk(:)*dk;
    end
        
    % Update capital 
    K_t = relax .* K_supply +(1 - relax).* K_t;
    
    % Update the interest rate 
    r_t = alpha.*A_t.*(L_t./K_t).^(1-alpha)-delta.*ones(It,1);

    % Update the wage
    w_t = (1-alpha).*A_t.*(K_t./L_t).^alpha;
    
    Kdist(n) = max(abs(K_supply - K_t));
        disp(['ITERATION = ', num2str(n)])
        disp(['Convergence criterion = ', num2str(Kdist(n))])
    if Kdist(n)<tolerance_S
        break 
    end
end

%{ 
--------------------------------------------------------------------------------------------------------------------------
Question 7: TFP sequence
--------------------------------------------------------------------------------------------------------------------------
%} 

% Initialize plot
set(gca, 'FontSize', 18)
    
    plot(tt, A_t, 'LineWidth', 2, 'LineStyle', '-', 'Color', [41/255, 182/255, 164/255])
        
    grid
        
        xlabel('Time, t', 'FontSize', 14)
        
        ylabel('TFP, A_t', 'FontSize', 14)

        set(gca, 'TickDir', 'out'); 
            box off; 

        xlim([min(tt) max(tt)]);
        
        legend('TFP', 'Location', 'best', 'FontSize', 14)

% Export graph
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_5/output/question_7.pdf', 'ContentType', 'vector');

%{ 
--------------------------------------------------------------------------------------------------------------------------
Question 8: interest rate
--------------------------------------------------------------------------------------------------------------------------
%} 

% Initialize plot
set(gca, 'FontSize', 18)

    plot(tt, r_t, 'LineWidth', 2, 'LineStyle', '-', 'Color', [41/255, 182/255, 164/255]); 
      
    grid
        
        xlabel('Time, t', 'FontSize', 18);
            
        ylabel('Interest rate, r_t', 'FontSize', 18);

        set(gca, 'TickDir', 'out'); 
            box off; 

        xlim([min(tt) max(tt)]);
        
        legend('Interest Rate', 'Location', 'best', 'FontSize', 14)
    

% Export graph
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_5/output/question_8_interest_rate.pdf', 'ContentType', 'vector');

%{ 
--------------------------------------------------------------------------------------------------------------------------
Question 8: wages
--------------------------------------------------------------------------------------------------------------------------
%} 

% Initialize plot
set(gca, 'FontSize', 18)

    plot(tt, w_t, 'LineWidth', 2, 'LineStyle', '-', 'Color', [41/255, 182/255, 164/255]); 
      
    grid
        
        xlabel('Time, t', 'FontSize', 18);
            
        ylabel('Wage, w_t', 'FontSize', 18);

        set(gca, 'TickDir', 'out'); 
            box off; 

        xlim([min(tt) max(tt)]);
        
        legend('Wages', 'Location', 'best', 'FontSize', 14)

% Export graph
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_5/output/question_8_wages.pdf', 'ContentType', 'vector');


%{ 
--------------------------------------------------------------------------------------------------------------------------
Question 8: Capital
--------------------------------------------------------------------------------------------------------------------------
%} 

% Initialize plot
set(gca, 'FontSize', 18)

    plot(tt, K_t, 'LineWidth', 2, 'LineStyle', '-', 'Color', [41/255, 182/255, 164/255]); 
      
    grid
        
        xlabel('Time, t', 'FontSize', 18);
            
        ylabel('Capital, K_t', 'FontSize', 18);

        set(gca, 'TickDir', 'out'); 
            box off; 

        xlim([min(tt) max(tt)]);
        
        legend('Capital', 'Location', 'best', 'FontSize', 14)

% Export graph
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_5/output/question_8_capital.pdf', 'ContentType', 'vector');
