%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: HJB_ramsey_implicit
% 
% Author: Kiyea Jin
% Date: Oct 27, 2024
%
% Description:
% This MATLAB script implements implicit method to solve the HJB equation
% of the deterministic Ramsey Growth Model using mixed method. Boundary 
% conditions are not yet considered in this script.
%
% Reference:
% HJB_NGM_implicit.m by Benjamin Moll
% ramsey_implicit.m by Pontus Rendahl
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-gamma))/(1-gamma)
% - Production function: f(k) = A*k^alpha
% - Relative risk aversion coefficient (gamma): 2
% - Discount rate (rho): 0.03
% - Depreciation rate (delta): 0.025
% - Elasticity of output with respect to capital (alpha): 1/3
% - Total fator productivity (A): 1
% - Delta = 1000 (Can be arbitrarily large in implicit method)
%
% Code Structure:
% 1. DEFINE PARAMETERS
% 2. INITIALIZE GRID POINTS
% 3. PRE-ITERATION INITIALIZATION
% 4. VALUE FUNCTION ITERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters_demo();

%% 2. INITIALIZE GRID POINTS

% Steady-state level of capital: f'(kss)=rho+delta
kss = ((p.rho+p.delta)/(p.A*p.alpha))^(1/(p.alpha-1));

% log(k_min) = log(kss)-p.klim
k_min = kss*exp(-p.klim); 
k_max = kss*exp(p.klim);

k = linspace(k_min, k_max, p.I)';
dk = (k_max-k_min)/(p.I-1);

%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the differential operator D such that DV=dV

    D = zeros(p.I, p.I);
    
    % Forward differencing for i=1
    D(1,1) = -1/dk; D(1,2) = 1/dk;
    
    % Backward differencing for i=I
    D(end,end-1) = -1/dk; D(end,end) = 1/dk;
    
    % Central differencing for i=2,...,I-1
    for i = 2:p.I-1
        D(i,i-1) = -0.5/dk; D(i,i+1) = 0.5/dk;
    end

% 3-2. Guess an initial value of the value function

    v0 = p.u(p.f(k) - p.delta*k)/p.rho;
    V = v0;   

% 3-3. Guess an initial value of the consumption
% Notes: This script does not consider the boundary conditions, so we'll
% impose an initial value of consumption 
    
    c0 = p.f(k);
    c = c0;

%% 4. VALUE FUNCTION ITERATION

tic;

for n = 1:p.maxit

    % 4-3. Compute the optimal savings
        s = p.f(k) - p.delta*k - c;

    % 4-4. Update the value function: V^(n+1) = [(rho+1/Delta)*I - SD]^(-1)[u(c) + 1/Delta*V^n]
    
        % B = [(rho+1/Delta)*I - SD]
        S = diag(s);
        B = (p.rho + 1/p.Delta)*eye(p.I) - S*D;

        % b = [u(c) + 1/Delta*V^n]
        b = p.u(c) + 1/p.Delta*V;

        % V^(n+1) = B^(-1)*b
        V_update = B\b;
        
        % Update the value function
        V_change = V_update - V;
        V = V_update;

    % 4-1. Compute the derivative of the value function
        dV = D*V;

    % 4-2. Compute the optimal consumption
        c = p.inv_mu(dV);

    % 4-5. Check convergence
          
        dist(n) = max(abs(V_change));

        if dist(n)<p.tol
        disp('Value function converged. Iteration = ')
        disp(n)
        break
        end
end;

toc;

%% 5. Saddle path

% Saddle path
figure;
p1 = plot(k, c, 'linewidth', 2);
set(gca, 'FontSize', 18)
xlabel('Capital, k','FontSize', 18)
ylabel('Consumption, c','FontSize',18)

% dk/dt = 0
k_null_c = p.f(k) - p.delta*k;
hold on
p2 = plot(k, k_null_c, 'linewidth', 2);
yy = get(gca, 'yLim');

% dc/dt = 0
p3 = plot([kss kss], yy, 'linewidth', 2);

legend1 = legend([p1,p2,p3], 'Saddle path', '\Delta k=0', '\Delta c=0');
set(legend1, 'Location', 'best', 'FontSize', 18)
hold off
