function p = define_parameters()

% This function defines the parameters needed for the HJB_Huggett.m script

%% Economic Parameters

    % Discount rate
    p.rho = 0.05;

    % Relative risk aversion coefficient
    p.sigma = 2;

    % Interest rate
    p.r = 0.03;

    % Income
    p.z_u = 0.1;
    p.z_e = 0.2;
    p.zz = [p.z_u, p.z_e];

    % Transition rates
    p.lambda_u = 0.02;
    p.lambda_e = 0.03;
    p.lambda = [p.lambda_u, p.lambda_e];
    
%% Economic Functions
    
    % Utility function
        % if sigma == 1
        % p.u = @(c) log(c);
        % else
        % p.u = @(c) (c.^(1-sigma))./(1-sigma);
        % end
    p.u = @(c) (c.^(1-p.sigma))./(1-p.sigma);

    % Marginal utility function
    p.mu = @(c) c.^(-p.sigma);

    % Inverse of marginal utility
    % Note: FOC: mu(c) = dv(a) -> c = inv_mu(dv)
    p.inv_mu = @(dv) dv.^(-1/p.sigma);

%% Grid Paramters

    % The lower bound of the state space (borrowing constraint)
    p.amin = -0.02;

    % The upper bound of the state space
    p.amax = 2;

    % The number of grid points
    p.I = 500;

%% Tuning Parameters
    
    % The maximum number of iteration
    p.maxit = 100;

    % Convergence criterion
    p.tol = 1e-8;

    % Step size (Delta)
    p.Delta = 1000;

end