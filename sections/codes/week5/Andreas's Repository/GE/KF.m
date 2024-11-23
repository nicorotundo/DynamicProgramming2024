function [gg] = KF(P, G)

%% 5. KF EQUATION

% 5-1. Solve for 0=gdot=P'*g

PT = P';
gdot_stacked = zeros(2*G.J,1);

% need to fix one value, otherwise matrix is singular
i_fix = 1;
gdot_stacked(i_fix)=.1;

row_fix = [zeros(1,i_fix-1),1,zeros(1,2*G.J-i_fix)];
AT(i_fix,:) = row_fix;

g_stacked = PT\gdot_stacked; 

% 5-2. Normalization

g_sum = g_stacked'*ones(2*G.J,1)*G.da;
g_stacked = g_stacked./g_sum;

% 5-3. Reshape

gg = reshape(g_stacked, G.J, 2);

end