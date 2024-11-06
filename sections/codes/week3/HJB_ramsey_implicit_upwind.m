%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: HJB_ramsey_implicit_upwind
% 
% Author: Kiyea Jin
% Date: Nov 1, 2024
%
% Description:
% This MATLAB script implements implicit method to solve the HJB equation
% of the deterministic Ramsey Growth Model using upwind scheme.
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
% - Try with rho = delta = 0.05
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

p = define_parameters();

%% 2. INITIALIZE GRID POINTS

% Steady-state level of capital: f'(kss)=rho+delta
kss = ((p.rho+p.delta)/(p.A*p.alpha))^(1/(p.alpha-1));

% log(k_min) = log(kss)-p.klim
k_min = kss*exp(-p.klim); 
k_max = kss*exp(p.klim);

k = linspace(k_min, k_max, p.I)';
dk = (k_max-k_min)/(p.I-1);

%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the forward and backward differential operator 
% Df such that Df*V=dVf and Db such that Db*V=dVb

    Df = zeros(p.I, p.I);
    for i = 1:p.I-1
        Df(i,i) = -1/dk; Df(i,i+1) = 1/dk;
    end
    %Df = sparse(Df);

    Db = zeros(p.I, p.I);
    for i = 2:p.I
        Db(i,i-1) = -1/dk; Db(i,i) = 1/dk;
    end
    %Db = sparse(Db);

% 3-2. Guess an initial value of the value function

    v0 = p.u(p.f(k))/p.rho;
    V = v0;   

% 3-3. Pre-allocate arrays for solutions

    dVf = zeros(p.I,1);
    dVb = zeros(p.I,1);
    c = zeros(p.I,1);

%% 4. VALUE FUNCTION ITERATION

tic;

for n = 1:p.maxit

    % 4-1. Compute the derivative of the value function
        dVf = Df*V;
        dVb = Db*V;

    % BOUNDARY CONDITIONS
        dVf(end) = p.mu(p.f(k_max) - p.delta*k_max); % k<=k_max is enforced which helps stability of the algorithm
        dVb(1)   = p.mu(p.f(k_min) - p.delta*k_min); % k>=k_min is enforced which helps stability of the algorithm

    % 4-2. Compute the optimal consumption
        cf = p.inv_mu(dVf);
        cb = p.inv_mu(dVb);

    % 4-3. Compute the optimal savings
        sf = p.f(k) - p.delta*k - cf;
        sb = p.f(k) - p.delta*k - cb;
   
    % UPWIND SCHEME
        If = sf>0;
        Ib = sb<0;
        I0 = 1-If-Ib;
        dV0 = p.mu(p.f(k) - p.delta*k);

        dV_upwind = dVf.*If + dVb.*Ib + dV0.*I0;

        c = p.inv_mu(dV_upwind);

    % 4-4. Update the value function: V^(n+1) = [(rho+1/Delta)*I - SD]^(-1)[u(c) + 1/Delta*V^n]
    
        % B = [(rho+1/Delta)*I - SD]
        Sf = diag(sf.*If);
        Sb = diag(sb.*Ib);
        SD = Sf*Df + Sb*Db;

        B = (p.rho + 1/p.Delta)*eye(p.I) - SD;

        % b = [u(c) + 1/Delta*V^n]
        b = p.u(c) + 1/p.Delta*V;

        % V^(n+1) = B^(-1)*b
        V_update = B\b;
        
        % Update the value function
        V_change = V_update - V;
        V = V_update;

    % 4-5. Check convergence
          
        dist(n) = max(abs(V_change));

        if dist(n)<p.tol
        disp('Value function converged. Iteration = ')
        disp(n)
        break
        end

end

toc;

%% 5. SADDLE PATH

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
