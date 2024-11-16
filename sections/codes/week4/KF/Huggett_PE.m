%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: HJB_Huggett_PE
% 
% Author: Kiyea Jin
% Date: Nov 7, 2024
%
% Description:
% This MATLAB script solves the HJB equation and the KF equation 
% of the Huggett model.
%
% Reference: Huggett_partialeq.m by Benjamin Moll
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-sigma))/(1-sigma)
% - Elasticity of intertemporal substitution (sigma): 1.2
% - Interest rate (r) : 0.035
% - Discount rate (rho): 0.05
% - Income: z = [z_u, z_e] = [0.1, 0.2];
% - Lambda: la = [la_u, la_e] = [1.5, 1];
% - Discrete grid of asset levels (a): -0.02 to 2
% - Borrowing constraint: a>=-0.02
% - Delta = 1000; (Can be arbitrarily large in implicit method)
%
% Code Structure:
% 1. DEFINE PARAMETERS
% 2. INITIALIZE GRID POINTS
% 3. PRE-ITERATION INITIALIZATION
% 4. VALUE FUNCTION ITERATION
% 5. KF EQUATION
% 6. GRAPHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters();

%% 2. INITIALIZE GRID POINTS

a = linspace(p.amin, p.amax, p.I)';
da = (p.amax-p.amin)/(p.I-1);

aa = [a, a]; % I*2 matrix

%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the forward and backward differential operator 
% Df such that Df*V=dVf and Db such that Db*V=dVb

    Df = zeros(p.I, p.I);
    for i = 1:p.I-1
        Df(i,i) = -1/da; Df(i,i+1) = 1/da;
    end
    Df = sparse(Df);

    Db = zeros(p.I, p.I);
    for i = 2:p.I
        Db(i,i-1) = -1/da; Db(i,i) = 1/da;
    end
    Db = sparse(Db);

% 3-2. Construct A_switch matrix

    A_switch = [speye(p.I).*(-p.lambda(1)), speye(p.I).*p.lambda(1);
                speye(p.I).*p.lambda(2), speye(p.I).*(-p.lambda(2))];

%A_switch = zeros(2*I, 2*I);
%for i=1:I
%    A_switch(i,i) = -lambda(1);
%    A_switch(i,i+I) = lambda(1);
%    A_switch(i+I,i) = lambda(2);
%    A_switch(i+I,i+I) = -lambda(2);
%end

% 3-3. Guess an initial value of the value function
    zz = ones(p.I, 1).*p.zz; % I*2 matrix
    
    % The value function of "staying put" 
    v0 = p.u(zz + p.r.*aa)./p.rho; 
    V = v0;

%% 4. VALUE FUNCTION ITERATION

for n=1:p.maxit

    % 4-1. Compute the derivative of the value function 
    dVf = Df*V;
    dVb = Db*V;

    % 4-2. Boundary conditions
    dVb(1,:) = p.mu(zz(1,:) + p.r.*aa(1,:)); % a>=a_min is enforced (borrowing constraint)
    dVf(end,:) = p.mu(zz(end,:) + p.r.*aa(end,:)); % a<=a_max is enforced which helps stability of the algorithm

    I_concave = dVb > dVf; % indicator whether value function is concave (problems arise if this is not the case)

    % 4-3. Compute the optimal consumption
    cf = p.inv_mu(dVf);
    cb = p.inv_mu(dVb);
    
    % 4-4. Compute the optimal savings
    sf = zz + p.r.*aa - cf;
    sb = zz + p.r.*aa - cb;

    % 4-5. Upwind scheme
    If = sf>0;
    Ib = sb<0;
    I0 = 1-If-Ib;
    dV0 = p.mu(zz + p.r.*aa); % If sf<=0<=sb, set s=0

    dV_upwind = If.*dVf + Ib.*dVb + I0.*dV0;

    c = p.inv_mu(dV_upwind);

    % 4-6. Update value function: 
    % Vj^(n+1) = [(rho + 1/Delta)*I - (Sj^n*Dj^n+A_switch)]^(-1)*[u(cj^n) + 1/Delta*Vj^n]
    
    V_stacked = V(:); % 2I*1 matrix
    c_stacked = c(:); % 2I*1 matrix

    % A = SD
    SD_u = spdiags(If(:,1).*sf(:,1), 0, p.I, p.I)*Df + spdiags(Ib(:,1).*sb(:,1), 0, p.I, p.I)*Db; % I*I matrix
    SD_e = spdiags(If(:,2).*sf(:,2), 0, p.I, p.I)*Df + spdiags(Ib(:,2).*sb(:,2), 0, p.I, p.I)*Db; % I*I matrix
    SD = [SD_u, sparse(p.I, p.I);
         sparse(p.I, p.I), SD_e]; % 2I*2I matrix
   
    % P = A + A_switch
    P = SD + A_switch;

    % B = [(rho + 1/Delta)*I - P]
    B = (p.rho + 1/p.Delta)*speye(2*p.I) - P; 

    % b = u(c) + 1/Delta*V
    b = p.u(c_stacked) + (1/p.Delta)*V_stacked;

    % V = B\b;
    V_update = B\b; % 2I*1 matrix
    V_change = V_update - V_stacked;
    V = reshape(V_update, p.I, 2); % I*2 matrix

    % 3-6. Convergence criterion
    dist(n) = max(abs(V_change));
    if dist(n)<p.tol
       disp('Value function converged. Iteration = ')
       disp(n)
       break
    end
end

toc;

%% 5. KF EQUATION

%% 5-1. SOLVE FOR 0=gdot=P'*g

PT = P';
PT_eigs = PT;
gdot_stacked = zeros(2*p.I,1);

% Since the eigenvector is only defined up to a scalar, we need to fix one
% value; otherwise matrix is singular.
i_fix = 1;

gdot_stacked(i_fix)=.1;

row_fix = [zeros(1,i_fix-1),1,zeros(1,2*p.I-i_fix)];
PT(i_fix,:) = row_fix;

g_stacked = PT\gdot_stacked; 

% Normalization
g_sum = g_stacked'*ones(2*p.I,1)*da;
g_stacked = g_stacked./g_sum;

% Reshape
gg = reshape(g_stacked, p.I, 2);

%% 5-2. KF EQUATION SOLVED WITH EIGS
% Notes: [V, D] = eigs(A, k, sigma) returns the k largest (in magnitude) eigenvalues D and 
% corresponding eigenvectors V of matrix A closest to the target sigma. 
% - A: The input matrix.
% - k: The number of eigenvalues (and corresponding eigenvectors) to compute.
% - If SIGMA is a real or complex scalar including 0, eigs finds the eigenvalues closest to SIGMA.

[g_stacked_eigs, eigval] = eigs(PT_eigs, 1, 0);
g_sum_eigs = g_stacked_eigs'*ones(2*p.I,1)*da;
g_stacked_eigs = g_stacked_eigs./g_sum_eigs;


%% 6. GRAPHS 

set(gca,'FontSize',14)
plot(dist,'LineWidth',2)
grid
xlabel('Iteration')
ylabel('||V^{n+1} - V^n||')

% Verr = c.^(1-s)/(1-s) + dV_Upwind.*(zz + r.*aa - c) + ones(I,1)*la.*(V_switch - V) - rho.*V;
% 
% set(gca,'FontSize',14)
% plot(a,Verr,'LineWidth',2)
% grid
% xlabel('k')
% ylabel('Error in HJB Equation')
% xlim([amin amax])

% 6-1. Optimal consumption 

set(gca, 'FontSize', 18)
plot(a, c, 'LineWidth', 2)
grid
xlabel('Wealth, a','FontSize', 14)
ylabel('Consumption, c_j(a)','FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

% 6-2. Optimal savings 

adot = zz + p.r.*aa - c;

set(gca, 'FontSize', 18)
plot(a, adot, a, zeros(1,p.I), '--k', 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Saving, s_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

% 6-3. Value function

set(gca, 'FontSize', 18)
plot(a, V, 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Value function, V_j(a)', 'FontSize', 14)
xlim([p.amin p.amax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)

% 6-4. Wealth distribution

set(gca, 'FontSize', 14)
plot(a, gg, 'LineWidth', 2)
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Densities, g_j(a)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([-0.02,-0.02], yy, '--k', 'LineWidth', 2)
text(-0.02, yy(1)-0.02*(yy(2) - yy(1)), '$\underline{a}$', 'HorizontalAlignment', 'center', 'FontSize', 15, 'Interpreter', 'latex')
xlim([-0.03 1])
legend('Unemployed', 'Employed', 'Borrowing Constraint', 'Location', 'best', 'FontSize', 14)