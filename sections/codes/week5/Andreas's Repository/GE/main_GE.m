%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: main_GE
% 
% Author: Kiyea Jin
% Date: Nov 16, 2024
%
% Description:
% This MATLAB script solves the general equilibrium of the Huggett model.
% The implementation leverages the SparseEcon repository developed by 
% Andreas Schaab and Allen T. Zhang (available at https://github.com/schaab-lab/SparseEcon).
%
% Reference:
% - Huggett_equilibrium_iterate.m by Benjamin Moll
% - Codes developed during the 2024 Coding Workshop by Andreas Schaab
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-sigma))/(1-sigma)
% - Elasticity of intertemporal substitution (sigma): 2
% - Discount rate (rho): 0.05
% - Income: z = [z_u, z_e] = [0.1, 0.2];
% - Lambda: la = [la_u, la_e] = [1.5, 1];
% - Discrete grid of asset levels (a): -0.15 to 5
% - Borrowing constraint: a>=-0.15
% - Delta = 1000; (Can be arbitrarily large in implicit method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc;

%% 0. Add the path of SparseEcon folder

addpath(genpath('../SparseEcon/lib'));
% addpath(genpath('SparseEcon/lib'));

%% 1. DEFINE PARAMETERS

p = define_parameters();

%% 2. INITIALIZE GRID POINTS

% a = linspace(p.amin, p.amax, p.I)';
% da = (p.amax-p.amin)/(p.I-1);
% 
% aa = [a, a]; % I*2 matrix

%{
setup_grid.m: Initializes grid structure

function G = setup_grid(n, surplus, min, max, varargin)

INPUTS:
- n:       Level of sparse grid
- surplus: Dimension-specific level surplus e.g., [0, 1, -1]
- min:     Dimension-specific minimums in economic units
- max:     Dimension-specific maximums in economic units

VARIABLE INPUTS:
- NamedDims: Cell of vectors of named dimensions
- Names:     Names for named dimensions
- DxxDims:   Vector of dimensions to compute dxx operators
- DxyDims:   (n x 2) matrix of dimensions to compute dxy operators

%}

G = setup_grid(p.l, 0, p.amin, p.amax, 'NamedDims', {1}, 'Names', {'a'});

% Notes:  G.a is the grid points. G.J is the number of grid points.

%% NEWTON'S METHOD

% Guess an initial value of the interest rate
r0 = p.r;

% Define asset supply S(r) as a function of r
S = @(r) stationary(r, G, p);

% Solve for equilibrium interest rate using fsolve
options = optimoptions("fsolve", "FunctionTolerance", p.tol_S);
r_equilibrium = fsolve(S, r0, options);

% Rerun stationary using the equilibrium interest rate
[~, ss] = stationary(r_equilibrium, G, p);

%% 6. GRAPHS 

% 6-1. Optimal consumption 

set(gca, 'FontSize', 18)
plot(G.a, ss.c, 'LineWidth', 2)
grid
xlabel('Wealth, a','FontSize', 14)
ylabel('Consumption, c_j(a)','FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

% 6-2. Optimal savings 

% adot = G.income - c;

set(gca, 'FontSize', 18)
plot(G.a, ss.s, G.a, zeros(1,G.J), '--k', 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Saving, s_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

% 6-3. Value function

set(gca, 'FontSize', 18)
plot(G.a, ss.V, 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Value function, V_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

% 6-4. Wealth distribution

set(gca, 'FontSize', 14)
plot(G.a, ss.gg, 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Densities, g_j(a)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([-0.15,-0.15], yy, '--k', 'LineWidth', 2)
text(-0.15, yy(1)-0.02*(yy(2) - yy(1)), '$\underline{a}$', 'HorizontalAlignment', 'center', 'FontSize', 15, 'Interpreter', 'latex')
xlim([-0.2 1])
legend('Unemployed', 'Employed', 'Borrowing Constraint', 'Location', 'best', 'FontSize', 14)