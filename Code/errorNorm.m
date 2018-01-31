function res = errorNorm(x)
% x = [thetas0; thetadots0; T];
z0 = x(1:(end-1));
T = z0(end);
n = length(z0)/2;
% parameters
p.l = ones(1, n); 
p.m = ones(1, n);
p.I = p.m .* (p.l).^2 ./ 12;
p.g = 9.8;

opts.RelTol = 1e-12; opts.AbsTol = 1e-12;
[t, zarray] = ode45(str2func(['pendulum_lagrange_', num2str(n)]), [0, 2*T], z0, opts, p);
f1 = interp1(t, zarray, T)';
f2 = interp1(t, zarray, 2*T)';

Q = diag([ones(1, n).*100, ones(1, n).*100]);

% res = norm(z0 - f1) + norm(z0 - f2) + norm(z0 - f3);
res = (z0 - f1)'*Q*(z0 - f1) + (z0 - f2)'*Q*(z0 - f2);
end