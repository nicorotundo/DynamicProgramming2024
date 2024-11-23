function [S, ss] = stationary(r, G, p)

%% HJB

[V, c, s, P] = HJB(r, G, p);

%% KF

[gg] = KF(P, G);

%% OUTPUT

S = gg(:,1)'*G.a*G.da + gg(:,2)'*G.a*G.da;

ss.V = V;
ss.gg = gg;
ss.c = c;
ss.s = s;

end