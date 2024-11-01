%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script: ramsey.m
%
% Author: Kiyea Jin
% Date: October 15, 2024
%
% Description:
% This script solves the Ramsey-Cass-Koopmans model
% using FD scheme with shooting algorithm.
%
% The model consists of two equations:
%   (dc/dt)/c = (f'(k) - rho - theta*g)/theta
%   dk/dt = f(k) - c - (n+g)k
% where an initial condition for capital K0 is provided,
% and intertemporal budget constraint with equality is imposed as a terminal condition.
%
% Parameters:
% - Discount rate (rho): 0.03
% - Inverse of IES (theta): 1
% - Technology growth rate (g): 0.02
% - Population growth rate (n): 0.02
% - Capital share (alpha): 1/3
% - TFP (A): 1
% - Initial boundary condition: K0 = 10
%
% Code Structure:
% 1. DEFINE PARAMETERS
% 2. INITIALIZE GRID POINTS
% 3. STEADY STATES
% 4. SHOOTING ALGORITHM
% 5. PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters();

%% 2. INITIALIZE GRID POINTS

t = linspace(p.tmin, p.tmax, p.I)';
dt = (p.tmax-p.tmin)/(p.I-1);

%% 3. STEADY STATES
% kss = ((rho + theta*g)/(alpha*A))^(1/(alpha-1))
% css = Ak^alpha - (n+g)*kss

kss = ((p.rho + p.theta * p.g) / (p.alpha * p.A))^(1 / (p.alpha - 1));
css = p.f(kss) - (p.n+p.g) * kss;

%% 4. SHOOTING ALGORITHM

        % 4-1. Find the initial value of consumption that satisfies terminal boundary condition
    
            tic;
            % Objective function that calculates the difference between terminal capital k(T) and steady-state kss
            % Note: k(T)-kss is a function of c0
            diff = @(c0) terminal_condition(c0, p.k0, kss, p.f, p.f_prime, p.rho, p.theta, p.g, p.n, dt, p.I);
            
            % Guess an initial value of consumption
            c0_guess = 1;

            % Use fsolve to find the initial consumption c0 that makes k(T) = k_ss
            % Note: X = fsolve(FUN, X0, options) starts at the matrix X0 and tries to solve the equations in FUN.
            % Set OPTIONS = optimoptions('fsolve','Algorithm','trust-region'), and then pass OPTIONS to fsolve.
            options = optimoptions('fsolve', 'TolFun', p.tol, 'Display', 'iter');
            c0 = fsolve(diff, c0_guess, options);

    % 4-2. Forward simulate with the updated initial consumption

    [k, c] = forward_simulate(c0, p.k0, p.f, p.f_prime, p.rho, p.theta, p.g, p.n, dt, p.I);
    toc;

%% 5. PLOT

% 5-1. Evolution of capital and consumption

figure;
subplot(2,1,1);
plot(t, k, 'r-', 'LineWidth', 2);
xlabel('Time'); ylabel('Capital k(t)');
title('Capital Accumulation over Time');

subplot(2,1,2);
plot(t, c, 'b-', 'LineWidth', 2);
xlabel('Time'); ylabel('Consumption c(t)');
title('Consumption Growth over Time');

%% 5-2. Phase diagram

% Saddle path
figure;
p1 = plot(k, c, 'linewidth', 2);
set(gca, 'FontSize', 18)
xlabel('Capital, k','FontSize', 18)
ylabel('Consumption, c','FontSize',18)

% dk/dt = 0
k_null_c = p.f(k) - (p.n+p.g)*k;
hold on
p2 = plot(k, k_null_c, 'linewidth', 2);
yy = get(gca, 'yLim');

% dc/dt = 0
p3 = plot([kss kss], yy, 'linewidth', 2);

legend1 = legend([p1,p2,p3], 'Saddle path', '\Delta k=0', '\Delta c=0');
set(legend1, 'Location', 'best', 'FontSize', 18)
hold off