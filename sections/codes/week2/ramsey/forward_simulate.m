function [k, c] = forward_simulate(c0, k0, f, f_prime, rho, theta, g, n, dt, I)

% This function solves the two differential equations using forward simulation.

    % Pre-allocate arrays for solution

    k = zeros(I, 1);  
    c = zeros(I, 1);  
    k(1) = k0;
    c(1) = c0;
    
    for i = 1:I-1
        
        % Euler equation for consumption growth: (c(i+1)-c(i))/dt = c(i)*(f'(k(i)-rho-theta*g)/theta
        c(i+1) = c(i) + dt * (f_prime(k(i)) - rho - theta * g) / theta * c(i);
        
        % Capital accumulation equation:(k(i+1)-k(i))/dt = f(k(i))-c(i)-(n+g)k(i)
        k(i+1) = k(i) + dt * (f(k(i)) - c(i) - (n + g) * k(i));
    end

end