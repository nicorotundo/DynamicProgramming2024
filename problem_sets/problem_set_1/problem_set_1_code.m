%{
--------------------------------------------------------------------------------------------------------------------------
Purpose: Coding part of problem set 1
Created: Nico Rotundo 2024-10-28
--------------------------------------------------------------------------------------------------------------------------
%}

%{ 
--------------------------------------------------------------------------------------------------------------------------
Define parameters
--------------------------------------------------------------------------------------------------------------------------
%}

% Terminal consumption at time T
C_T = 2.0;

% Initial interest rate 
r_0 = 0.05;

% Growth rate of interest rate 
alpha = 0.01;

%{
% Interest rate function as a handle function of time t
function r = interest_rate_function(t, r_0, alpha)
    r = r_0 + alpha * t;
end
%}

% Inverse of the intertemporal elasticity of substitution
theta = 2.0;

% Time discount factor
rho = 0.03;

% Terminal time 
T = 10;

% Display the result of interest rate function at specific time to check 
%fprintf('The interest rate at t = %d is: %.3f\n', 20, interest_rate_function(20, r_0, alpha));

% Checking how assert statements work in this matlab
%assert(interest_rate_function(20, r_0, alpha) == .25);


%{ 
--------------------------------------------------------------------------------------------------------------------------
1. Solve the equation analytically using the integrating factor method. Then, plot the analytical solution as a function 
    of time t and consumption C(t)

    - solved equation: C(t) = C_T*exp[(1/theta)*[(r_0-rho)(t-T)+(alpha/2)(t^2-T^2)]]
--------------------------------------------------------------------------------------------------------------------------
%}

% Time vector from 0 to T with 100 points
t = linspace(0, T, 100);

% Define value for consumption at time t 
C_t_analytical = C_T * exp((1 / theta) * ((r_0 - rho) * (t - T) + (alpha / 2) * (t.^2 - T^2)));

% Generate figure for analytical solution
figure;
plot(t, C_t_analytical, 'LineWidth', 2, 'Color', [41/255, 182/255, 164/255]);
    xlabel('Time (t)');
    set(gca, 'TickDir', 'out'); 
    box off; 
    ylabel('Consumption C(t)');
    title('Analytical Solution of C(t) over Time');
    grid on;

% Export
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_1/output/question_1.pdf', 'ContentType', 'vector');

%{ 
--------------------------------------------------------------------------------------------------------------------------
2. Solve the equation numerically using the finite difference method with 100 time steps. Plot the numerical and analytical 
    solutions together in a single figure for comparison.

    - equation: C(t-dt) = (theta-dt(r_t-rho))/theta * C(t)
--------------------------------------------------------------------------------------------------------------------------
%}

% Time step (\triangle t)
dt = 0.1;           

% Number of time steps
N = T / dt;          

% Vector to store consumption values
C_t_numeric = zeros(1, N);   

% Set terminal consumption
C_t_numeric(end) = C_T;

% Time vector from 0 to T with 100 points
t = linspace(1, T, N);

% Time stepping using the finite difference method
for n = N:-1:2
    
    % Compute consumption at the next time step using backward difference
    C_t_numeric(n-1) = (theta - dt*(r_0 + alpha*n-rho))/theta * C_t_numeric(n);

end

% Generate figure for numeric solution
figure;
    plot(t, C_t_numeric, 'LineWidth', 2, 'Color', [250/255, 165/255, 35/255]);
    xlabel('Time (t)');
    set(gca, 'TickDir', 'out'); 
    box off; 
    ylabel('Consumption C(t)');
    title('Numerical Solution of C(t) over Time');
    grid on;

% Generate both analytical and numeric sultions on the same plot 
figure;
   
    % Plot analytical solution
    plot(t, C_t_analytical, 'LineWidth', 2, 'Color', [41/255, 182/255, 164/255]);

    hold on;

    % Plot numerical solution
    plot(t, C_t_numeric, 'LineWidth', 2, 'Color', [250/255, 165/255, 35/255]);

    % Add labels and title
    xlabel('Time (t)');
    ylabel('Consumption C(t)');
    title('Analytical vs Numerical Solutions of C(t) over Time');

    % Add a legend
    legend({'Analytical Solution', 'Numerical Solution'}, 'Location', 'southeast');

    % Adjust plot settings
    set(gca, 'TickDir', 'out');
    box off;
    grid on;

    hold off;

% Export
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_1/output/question_2.pdf', 'ContentType', 'vector');

%{ 
--------------------------------------------------------------------------------------------------------------------------
3. Repeat the numerical solution using 10 time steps. Plot the numerical and analytical solutions together in a single 
    figure.
--------------------------------------------------------------------------------------------------------------------------
%}

% Time step (\triangle t)
dt = 1;           

% Number of time steps
N = T / dt;          

% Vector to store consumption values
C_t_numeric = zeros(1, N);   

% Set terminal consumption
C_t_numeric(end) = C_T;

% Time vector from 0 to T with 100 points
t_10 = linspace(1, T, N);

% Time stepping using the finite difference method
for n = N:-1:2
    
    % Compute consumption at the next time step using backward difference
    C_t_numeric(n-1) = (theta - dt*(r_0 + alpha*n-rho))/theta * C_t_numeric(n);

end

% Generate both analytical and numeric sultions on the same plot 
figure;
   
    % Plot analytical solution
    plot(t, C_t_analytical, 'LineWidth', 2, 'Color', [41/255, 182/255, 164/255]);

    hold on;

    % Plot numerical solution
    plot(t_10, C_t_numeric, 'LineWidth', 2, 'Color', [250/255, 165/255, 35/255]);

    % Add labels and title
    xlabel('Time (t)');
    ylabel('Consumption C(t)');
    title('Analytical vs Numerical Solutions of C(t) over Time');

    % Add a legend
    legend({'Analytical Solution', 'Numerical Solution'}, 'Location', 'southeast');

    % Adjust plot settings
    set(gca, 'TickDir', 'out');
    box off;
    grid on;

    hold off;

% Export
exportgraphics(gcf, '/Users/nicorotundo/Documents/GitHub/DynamicProgramming2024/problem_sets/problem_set_1/output/question_3.pdf', 'ContentType', 'vector');
