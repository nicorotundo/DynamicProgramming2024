function difference = terminal_condition(c0, k0, kss, f, f_prime, rho, theta, g, n, dt, I)

% This function calculates the difference between terminal capital k(T) and
% steady-state kss.

    [k, ~] = forward_simulate(c0, k0, f, f_prime, rho, theta, g, n, dt, I);
    k_T = k(end);  % Terminal capital k(T)
    
    % The difference between terminal k(T) and steady-state k_ss
    difference = k_T - kss;  
end