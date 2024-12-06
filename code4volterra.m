% Input
Estart = -10; E0 = -1;
Tend = 40;
dtau = 0.01;

% LHS equations
E = @(t) Estart + t;
f = @(t) sqrt(pi)./(1+exp(-(E(t) - E0)));

% Integral equations
Kfun = @(t, tau) 2*sqrt(t - tau);

% Grid
tvec = 0:dtau:Tend; tvecs = tvec; tvecs(1) = [];

% Set up system
LHSvec = f(tvecs)';

N = length(tvecs); K = zeros(N,N); checkmat = K;
for ii = 1:N
    for jj = 1:ii
        K(ii, jj) = min(jj, 2) * Kfun(tvecs(ii), tvec(jj));
    end
end
% Solve for I' (did not implement forward substitution yet
Idiff = K\LHSvec;

% Solve for I
I = zeros(N,1); I(2) = dtau / 2 * (Idiff(1) + Idiff(2));
for ii = 3:N
    I(ii) = I(ii - 1) + dtau / 2 * (Idiff(ii - 1) + Idiff(ii));
end

% Plot results
plot(tvec(1:(end-1)), I*2/dtau)
    
        
