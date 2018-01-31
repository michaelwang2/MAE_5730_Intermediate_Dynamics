%% Lagrange Pendulum
clear; clc; close;
%%%%%%%%%
n = 3; %%% number of links
%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
animation = 't'; %% % to play animation or display energy conservation
periodic = 'f';  %% % to simulate periodic motion
energy = 't';    %% % to display energy conservation next to animation
%%%%%%%%%%%%%%%%%%%
method = 'Lagrange';

% parameters
tspan = 0:0.01:50;
opts.RelTol = 1e-13; opts.AbsTol = 1e-13;
p.l = ones(1, n); 
p.m = ones(1, n);
p.I = p.m .* (p.l).^2 ./ 12;
p.g = 9.8;

% integration
% z = [thetas, thetadots]
if periodic == 't'
    load(['PendPer', num2str(n), '_4.mat']);
    z0 = x(1:(end-1));
else
    z0 = zeros(n*2, 1); 
    z0(1:n) = pi/2;
end
[t, Lzarray] = ode45(str2func(['pendulum_lagrange_', num2str(n)]), tspan, z0, opts, p);
[Ex, Ey] = thetas2pos(n, t, p, Lzarray);

if periodic == 't'
    save(['Lagrange_Pendulum_Periodic_', num2str(n),'_4'], 't', 'Lzarray');
else
    save(['Lagrange_Pendulum_', num2str(n)], 't', 'Lzarray');
end

%% plot / animate for Lagrange Pendulum Only (to prevent animation slow down)
nl = sum(p.l);
if animation == 't'
    animatePendulum(t, Ex, Ey, Lzarray, p, method, energy);
else
    % conservation of energy
    [T, V, H, M] = energyPendulum(t, Lzarray(:, 1:n), Lzarray(:, (n+1):end), p);
    TT = sum(T, 2);
    VV = sum(V, 2);
    figure; 
    plot(t, TT+VV, 'k-');
    grid on; box on; 
    xlabel('Time [s]');
    yy = ylabel('Energy [J]', 'Rotation', 0);
    set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    title('Total Energy Lagrange');
end

%% N-E Pendulum 3-linked
clear; clc; close;
%%%%%%%%%%%%%%%%%%%
animation = 't'; %% % to play animation or display energy conservation
periodic = 'f';  %% % to simulate periodic motion
energy = 't';    %% % to display energy conservation next to animation
%%%%%%%%%%%%%%%%%%%
method = 'Newton-Euler';

% parameters
tspan = 0:0.01:50;
opts.RelTol = 1e-13; opts.AbsTol = 1e-13;
p.l = ones(1, 3);
p.m = ones(1, 3);
p.I = p.m .* (p.l).^2 ./ 12;
p.g = 9.8;

% integration
% z = [thetas, thetadots]
if periodic == 't'
    load(['PendPer', num2str(3), '.mat']);
    z0 = x(1:(end-1));
else
    z0 = zeros(6, 1); 
    z0(1:3) = pi/2;
end
[t, NEzarray] = ode45(@pendulum_NE_3, tspan, z0, opts, p);
[Ex, Ey] = thetas2pos(3, t, p, NEzarray);

save('Newton_Euler_Pendulum_3', 't', 'NEzarray');

% animation
if animation == 't'
    animatePendulum(t, Ex, Ey, NEzarray, p, method, energy);
else
    % conservation of energy
    [T, V, H, M] = energyPendulum(t, NEzarray(:, 1:3), NEzarray(:, 4:6), p);
    TT = sum(T, 2);
    VV = sum(V, 2);
    figure; 
    plot(t, TT+VV, 'k-');
    grid on; box on; 
    xlabel('Time [s]');
    yy = ylabel('Energy [J]', 'Rotation', 0);
    set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    title('Total Energy Newton-Euler');
end

%% DAE Pendulum 3-linked
clear; clc; close;
%%%%%%%%%%%%%%%%%%%
animation = 't'; %% % to play animation or display energy conservation
%%%%%%%%%%%%%%%%%%%
Title = 'DAE: 3-Linked Pendulum';

% parameters
tspan = 0:0.01:50;
opts.RelTol = 1e-13; opts.AbsTol = 1e-13;
p.l = ones(1, 3); 
p.m = ones(1, 3);
p.I = p.m .* (p.l).^2 ./ 12;
p.g = 9.8;

% integration
% z = [x1 y1 x1d y1d x2 y2 x2d y2d x3 y3 x3d y3d th1 th2 th3 thd1 thd2 thd3]
z0 = zeros(18, 1);
thetas = ones(3, 1).*(pi/2);
z0(13:15) = thetas;
z0(1:2) = p.l(1)/2 * [cos(thetas(1)); sin(thetas(1))];
z0(5:6) = p.l(1) * [cos(thetas(1)); sin(thetas(1))] + p.l(2)/2 * [cos(thetas(2)); sin(thetas(2))];
z0(9:10) = p.l(1) * [cos(thetas(1)); sin(thetas(1))] + p.l(2) * [cos(thetas(2)); sin(thetas(2))] + p.l(3)/2 * [cos(thetas(3)); sin(thetas(3))];
[t, DAEzarray] = ode45(@pendulum_DAE_3, tspan, z0, opts, p);
pos1 = DAEzarray(:, 1:2);
pos2 = DAEzarray(:, 5:6);
pos3 = DAEzarray(:, 9:10);
ths = DAEzarray(:, 13:15);
thds = DAEzarray(:, 16:18);
[Ex, Ey] = thetas2pos(3, t, p, [ths, thds]);

save('DAE_Pendulum_3', 't', 'DAEzarray');

% animation
if animation == 't'
    [KE, PE] = animatePendulumDAE(t, DAEzarray, p, Title, z0, animation);
else
    [KE, PE] = animatePendulumDAE(t, DAEzarray, p, Title, z0, animation);
    TT = sum(KE, 2);
    VV = sum(PE, 2);
    figure; 
    plot(t, TT+VV, 'k-');
    grid on; box on; 
    xlabel('Time [s]');
    yy = ylabel('Energy [J]', 'Rotation', 0);
    set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    title('Total Energy DAE');
end

%% DAE 4-bar Linkage 
clear; clc; close;
%%%%%%%%%%%%%%%%%%%
animation = 't'; %% % to play animation or display energy conservation
%%%%%%%%%%%%%%%%%%%
Title = 'DAE: 4-Bar Linkage';

% parameters
tspan = 0:0.01:50;
opts.RelTol = 1e-13; opts.AbsTol = 1e-13;
p.l = [1, 2, 1];
p.m = ones(1, 3);
p.I = p.m .* (p.l).^2 ./ 12;
p.g = 9.8;

% integration
% z = [x1 y1 x1d y1d x2 y2 x2d y2d x3 y3 x3d y3d th1 th2 th3 thd1 thd2 thd3]
z0 = zeros(18, 1);
del = deg2rad(30.01);
thetas = [del; pi - del; del];
z0(13:15) = thetas;
z0(1:2) = p.l(1)/2 * [cos(thetas(1)); sin(thetas(1))];
z0(5:6) = p.l(1) * [cos(thetas(1)); sin(thetas(1))] + p.l(2)/2 * [cos(thetas(2)); sin(thetas(2))];
z0(9:10) = p.l(1) * [cos(thetas(1)); sin(thetas(1))] + p.l(2) * [cos(thetas(2)); sin(thetas(2))] + p.l(3)/2 * [cos(thetas(3)); sin(thetas(3))];
[t, DAEzarray] = ode45(@FourBarLinkage_DAE, tspan, z0, opts, p);
pos1 = DAEzarray(:, 1:2);
pos2 = DAEzarray(:, 5:6);
pos3 = DAEzarray(:, 9:10);
ths = DAEzarray(:, 13:15);
thds = DAEzarray(:, 16:18);
[Ex, Ey] = thetas2pos(3, t, p, [ths, thds]);

save('DAE_4BarLinkage', 't', 'DAEzarray');

% animation
if animation == 't'
    [KE, PE] = animatePendulumDAE(t, DAEzarray, p, Title, z0, animation);
else
    [KE, PE] = animatePendulumDAE(t, DAEzarray, p, Title, z0, animation);
    TT = sum(KE, 2);
    VV = sum(PE, 2);
    figure; 
    plot(t, TT+VV, 'k-');
    grid on; box on; 
    xlabel('Time [s]');
    yy = ylabel('Energy [J]', 'Rotation', 0);
    set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    title('Total Energy 4-Bar Linkage');
end

%% Periodic Solution for n-linked Pendulum
clear; clc; close;
n = 4;

% parameters
p.l = ones(1, n); 
p.m = ones(1, n);
p.I = p.m .* (p.l).^2 ./ 12;
p.g = 9.8;

% constraints
% x = [thetas0; thetadots0; T];
A = []; b = []; Aeq = []; beq = []; nonlcon = [];
th = -pi/3; thd = -2*pi/2; 
lb = [ones(n, 1).*(th); ones(n, 1).*(thd); 3];
ub = -lb; ub(end) = 10;

% solver options
options = optimoptions('fmincon','Display','iter');
options.Algorithm = 'active-set';
options.MaxIter = 5e5; % Max number of iterations
options.MaxFunEvals  = 5e5; % Max number of function evaluations
options.TolX = 1e-12;
options.TolCon = 1e-12;
options.TolFun = 1e-12;

% initial guess (randomize within bounds)
x0 = zeros(n*2 + 1, 1);
% x0(1:n) = rand(n, 1).*(ub(1)*2) - ub(1).*ones(n, 1);
% x0((n+1):(end-1)) = rand(n, 1).*(ub(n+1)*2) - ub(n+1).*ones(n, 1);
x0(1:n) = normrnd(zeros(n, 1), ub(1)/3);
x0((n+1):(end-1)) = normrnd(zeros(n, 1), ub(n+1)/3);
x0(end) = rand()*(ub(end) - lb(end)) + lb(end); 

% fmincon
[x,fval] = fmincon(@errorNorm,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
save(['PendPer', num2str(n), '_2.mat'], 'x');
