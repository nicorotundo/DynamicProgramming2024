function p = define_parameters()

% This function defines the parameters needed for the HJB_ramsey_implicit.m script

%% Economic Parameters

    % Discount rate
    p.rho = 0.03;

    % Relative risk aversion coefficient
    p.gamma = 2;

    % Depreciation rate
    p.delta = 0.025;

    % Capital share
    p.alpha = 1/3;

    % TFP
    p.A = 1;

%% Economic Functions
    
    % Utility function
        % if gamma == 1
        % p.u = @(c) log(c);
        % else
        % p.u = @(c) (c.^(1-gamma))./(1-gamma);
        % end
    p.u = @(c) (c.^(1-p.gamma))./(1-p.gamma);

    % Marginal utility function
    p.mu = @(c) c.^(-p.gamma);

    % Inverse of marginal utility
    % Note: FOC: mu(c) = dv(a) -> c = inv_mu(dv)
    p.inv_mu = @(dv) dv.^(-1/p.gamma);

    % Production function
    p.f = @(k) p.A * k.^p.alpha;

%% Grid Paramters

    % The number of time steps
    p.I = 2000;

    % Higher klim implies broader coverage
    p.klim = 1.5;

%% Tuning Parameters
    
    % The maximum number of iteration
    p.maxit = 1000;

    % Convergence criterion, and
    p.tol = 1e-8;

    % Step size (Delta)
    p.Delta = 1000;

end