function p = define_parameters_initial()

% This function sets the parameters required to compute the initial equilibrium
% for use in the Huggett_transition.m script.

%% Economic Parameters
    
    % Relative risk aversion
    p.gamma = 2;

    % Discount rate
    p.rho = 0.05;
    
    % Income process
    p.z_u = 0.1;
    p.z_e = 0.2;
    p.zz = [p.z_u, p.z_e];

    % Probability density   
    p.lambda_u = 0.6;
    p.lambda_e = 0.5;
    p.lambda = [p.lambda_u, p.lambda_e];
    
%% Economic Functions
    
    % Utility funtion
    p.u = @(c) c.^(1-p.gamma)/(1-p.gamma);

    % Marginal utility function
    p.mu = @(c) c.^(-p.gamma);

    % FOC: mu(c)=dV -> c=inv_mu(dV)
    p.inv_mu = @(dV) dV.^(-1/p.gamma);

%% Grid Parmaters

    p.amin = -0.15;
    p.amax = 4;

    % The number of grid points
    p.I = 1000;

%% Tuning parameters

    % Step size: can be arbitrarily large in implicit method
    p.Delta = 1000;

    % The maximum number of value function iterations
    p.maxit = 100;

    % Tolerance for value function iterations
    p.tol = 10^(-6);

    %% TUNING PARAMETERS FOR INTEREST RATE ITERATION
    
    % The maximum number of interest rate iterations
    p.Nr = 40;

    % Tolerance for interest rate iterations
    p.tol_S = 10^(-8);

end