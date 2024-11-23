function [V, c, s, P] = HJB(r, G, p)

%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the forward and backward differential operator 
% Df such that Df*V=dVf and Db such that Db*V=dVb
    % 
    % Df = zeros(p.I, p.I);
    % for i = 1:p.I-1
    %     Df(i,i) = -1/da; Df(i,i+1) = 1/da;
    % end
    % Df = sparse(Df);
    % 
    % Db = zeros(p.I, p.I);
    % for i = 2:p.I
    %     Db(i,i-1) = -1/da; Db(i,i) = 1/da;
    % end
    % Db = sparse(Db);

    %% DIFFERENTIAL OPERATOR + BOUNDARY CONDITIONS

    %{
    gen_FD.m: Constructs FD operators for given boundary conditions

    function G = gen_FD(G, BC, name)

    INPUTS:
    - G: Grid struct
    - BC: Boundary condition inputs. This is a (d x 1) cell, where each
          element is a struct with the following properties
        - left.type: type of left boundary condition
        - right.type: type of right boundary condition
        - left.f: values associated with left boundary condition
        - right.f: values associated with right boundary condition
    - name: (Optional) name associated with boundary condition, useful when
            using multiple boundary conditions for the same grid
    %}

    G.income = p.zz + r.*G.a; 

    left_bound  = p.mu(G.income(G.a == p.amin, :));
    right_bound = p.mu(G.income(G.a == p.amax, :));
    for j = 1:2
        BC{1}.left.type = 'VNB'; BC{1}.right.type = 'VNF'; % Von Nuemann Boundary Condition
        BC{1}.left.f  = @(points) left_bound(j) * ones(size(points, 1), 1);
        BC{1}.right.f = @(points) right_bound(j) * ones(size(points, 1), 1);
        G = gen_FD(G, BC, num2str(j));
    end

% 3-2. Construct A_switch matrix

    A_switch = [speye(G.J).*(-p.lambda(1)), speye(G.J).*p.lambda(1);
                speye(G.J).*p.lambda(2), speye(G.J).*(-p.lambda(2))];

% 3-3. Guess an initial value of the value function
    
    % The value function of "staying put"
    
    v0 = p.u(G.income)./p.rho; 
    V = v0;

%% 4. VALUE FUNCTION ITERATION

for n=1:p.maxit

    %% 4-1. Compute the derivative of the value function 
    % dVf = Df*V;
    % dVb = Db*V;

    % function deriv = deriv_sparse(G, f, k, operator, name)

    dVf =zeros(G.J, 2);
    dVb =zeros(G.J, 2);
    for j = 1:2
        dVf(:, j) = deriv_sparse(G, V(:,j), 1, 'D1F', num2str(j));
        dVb(:, j) = deriv_sparse(G, V(:,j), 1, 'D1B', num2str(j));
    end
 
    %% 4-2. Boundary conditions
    % dVb(1,:) = p.mu(zz(1,:) + r.*aa(1,:)); % a>=a_min is enforced (borrowing constraint)
    % dVf(end,:) = p.mu(zz(end,:) + r.*aa(end,:)); % a<=a_max is enforced which helps stability of the algorithm
    
    % dVf(G.a == p.amax, :) = p.mu(G.income(G.a == p.amax, :));
    % dVb(G.a == p.amin, :) = p.mu(G.income(G.a == p.amin, :));

    I_concave = dVb > dVf; % indicator whether value function is concave (problems arise if this is not the case)

    % 4-3. Compute the optimal consumption
    cf = p.inv_mu(dVf);
    cb = p.inv_mu(dVb);
    
    % 4-4. Compute the optimal savings
    sf = G.income - cf;
    sb = G.income - cb;

    % 4-5. Upwind scheme
    If = sf>0;
    Ib = sb<0;
    I0 = 1-If-Ib;
    dV0 = p.mu(G.income); % If sf<=0<=sb, set s=0

    dV_upwind = If.*dVf + Ib.*dVb + I0.*dV0;

    c = p.inv_mu(dV_upwind);
    s = G.income - c;

    % 4-6. Update value function: 
    % Vj^(n+1) = [(rho + 1/Delta)*I - (Sj^n*Dj^n+A_switch)]^(-1)*[u(cj^n) + 1/Delta*Vj^n]
    
    V_stacked = V(:); % 2I*1 matrix
    c_stacked = c(:); % 2I*1 matrix

    % A = SD
    % SD_u = spdiags(If(:,1).*sf(:,1), 0, p.I, p.I)*Df + spdiags(Ib(:,1).*sb(:,1), 0, p.I, p.I)*Db; % I*I matrix
    % SD_e = spdiags(If(:,2).*sf(:,2), 0, p.I, p.I)*Df + spdiags(Ib(:,2).*sb(:,2), 0, p.I, p.I)*Db; % I*I matrix

    %% Generator: A
    % Aa{1} = FD_operator(G, drift, volatility, dimension, type);
    s = If.*sf + Ib.*sb;
    
    Aa{1} = FD_operator(G, s(:,1), zeros(G.J, 1), 1, '1'); % earning type 1
    Aa{2} = FD_operator(G, s(:,2), zeros(G.J, 1), 1, '2'); % earning type 2

    SD = [Aa{1}, sparse(G.J, G.J);
         sparse(G.J, G.J), Aa{2}]; % 2I*2I matrix 

    % SD = blkdiag(Aa{1}, Aa{2});

    % P = A + A_switch
    P = SD + A_switch;

    % B = [(rho + 1/Delta)*I - P]
    B = (p.rho + 1/p.Delta)*speye(2*G.J) - P; 

    % b = u(c) + 1/Delta*V
    b = p.u(c_stacked) + (1/p.Delta)*V_stacked;

    % V = B\b;
    V_update = B\b; % 2I*1 matrix
    V_change = V_update - V_stacked;
    V = reshape(V_update, G.J, 2); % I*2 matrix

    % 3-6. Convergence criterion
    dist(n) = max(abs(V_change));
    if dist(n)<p.tol
       disp('Value function converged. Iteration = ')
       disp(n)
       break
    end
end

end