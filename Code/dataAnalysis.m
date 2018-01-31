%% Data Analysis for Triple Pendulum
clear; close; clc;
Lag = load('Lagrange_Pendulum_3.mat');
NE = load('Newton_Euler_Pendulum_3.mat');
DAE_pend = load('DAE_Pendulum_3.mat');

t = Lag.t;
lag_ths = Lag.Lzarray(:, 1:3);
NE_ths = NE.NEzarray(:, 1:3);
DAE_ths = DAE_pend.DAEzarray(:, 13:15);
lag_thds = Lag.Lzarray(:, 4:6);
NE_thds = NE.NEzarray(:, 4:6);
DAE_thds = DAE_pend.DAEzarray(:, 16:18);

% cutoff point
tend = 19;
ind = find(t > tend, 1, 'first');
lag_ths = lag_ths(1:ind, :);
NE_ths = NE_ths(1:ind, :);
DAE_ths = DAE_ths(1:ind, :);
lag_thds = lag_thds(1:ind, :);
NE_thds = NE_thds(1:ind, :);
DAE_thds = DAE_thds(1:ind, :);
t = t(1:ind);

figure; 
for i = 1:3
    subplot(3, 2, 2*i-1); hold on;
    if i == 1
        title('Error Between Lagrange and Newton-Euler');
    end
    plot(t, rad2deg(abs(wrapToPi(lag_ths(:, i)) - wrapToPi(NE_ths(:, i)))), 'k-');
    yy = ylabel(['$$\delta\theta_', num2str(i), '[\circ]$$'], 'Rotation', 0, 'Interpreter', 'latex');
    set(yy, 'Units', 'Normalized', 'Position', [-0.18, 0.35, 0]);
    grid on; box on;
    if i == 3
        xlabel('Time [s]');
    end
    hold off;
    subplot(3, 2, 2*i); hold on;
    plot(t, rad2deg(abs(wrapToPi(lag_thds(:, i)) - wrapToPi(NE_thds(:, i)))), 'k-');
    yy = ylabel(['$$\delta\dot{\theta}_', num2str(i), '[\circ/s]$$'], 'Rotation', 0, 'Interpreter','latex');
    set(yy, 'Units', 'Normalized', 'Position', [-0.18, 0.35, 0]);
    grid on; box on;
    if i == 3
        xlabel('Time [s]');
    end
    hold off;
end


figure; 
for i = 1:3
    subplot(3, 2, 2*i-1); hold on;
    if i == 1
        title('Error Between Lagrange and DAE');
    end
    plot(t, rad2deg(abs(wrapToPi(lag_ths(:, i)) - wrapToPi(DAE_ths(:, i)))), 'k-');
    yy = ylabel(['$$\delta\theta_', num2str(i), '[\circ]$$'], 'Rotation', 0, 'Interpreter', 'latex');
    set(yy, 'Units', 'Normalized', 'Position', [-0.18, 0.35, 0]);
    grid on; box on;
    if i == 3
        xlabel('Time [s]');
    end
    hold off;
    subplot(3, 2, 2*i); hold on;
    plot(t, rad2deg(abs(wrapToPi(lag_thds(:, i)) - wrapToPi(DAE_thds(:, i)))), 'k-');
    yy = ylabel(['$$\delta\dot{\theta}_', num2str(i), '[\circ/s]$$'], 'Rotation', 0, 'Interpreter','latex');
    set(yy, 'Units', 'Normalized', 'Position', [-0.18, 0.35, 0]);
    grid on; box on;
    if i == 3
        xlabel('Time [s]');
    end
    hold off;
end

figure; 
for i = 1:3
    subplot(3, 2, 2*i-1); hold on;
    if i == 1
        title('Error Between DAE and Newton-Euler');
    end
    plot(t, rad2deg(abs(wrapToPi(DAE_ths(:, i)) - wrapToPi(NE_ths(:, i)))), 'k-');
    yy = ylabel(['$$\delta\theta_', num2str(i), '[\circ]$$'], 'Rotation', 0, 'Interpreter', 'latex');
    set(yy, 'Units', 'Normalized', 'Position', [-0.18, 0.35, 0]);
    grid on; box on;
    if i == 3
        xlabel('Time [s]');
    end
    hold off;
    subplot(3, 2, 2*i); hold on;
    plot(t, rad2deg(abs(wrapToPi(DAE_thds(:, i)) - wrapToPi(NE_thds(:, i)))), 'k-');
    yy = ylabel(['$$\delta\dot{\theta}_', num2str(i), '[\circ/s]$$'], 'Rotation', 0, 'Interpreter','latex');
    set(yy, 'Units', 'Normalized', 'Position', [-0.18, 0.35, 0]);
    grid on; box on;
    if i == 3
        xlabel('Time [s]');
    end
    hold off;
end

%% plot periodic data
clear; close; clc;
ll = load('Lagrange_Pendulum_Periodic_3_4.mat');
%ll2 = load('Lagrange_Pendulum_Periodic_4_2.mat');
len = length(ll.t)/2;
n = length(ll.Lzarray(1, :))/2;
figure;
for i = 1:n
    subplot(n, 2, 2*i-1); hold on;
    if i == 1
        title('Angle');
    end
    plot(ll.t(1:round(len/2)), rad2deg(ll.Lzarray(1:round(len/2), i)), 'k-');
    %plot(ll2.t(1:round(len/2)), rad2deg(ll2.Lzarray(1:round(len/2), i)), 'b-');
    yy = ylabel(['$$\theta_', num2str(i), '[\circ]$$'], 'Rotation', 0, 'Interpreter', 'latex');
    set(yy, 'Units', 'Normalized', 'Position', [-0.18, 0.35, 0]);
    grid on; box on;
    if i == n
        xlabel('Time [s]');
    end
    hold off;
    subplot(n, 2, 2*i); hold on;
    if i == 1
        title('Angular Velocity');
    end
    plot(ll.t(1:round(len/2)), rad2deg(ll.Lzarray(1:round(len/2), i+n)), 'k-');
    %plot(ll2.t(1:round(len/2)), rad2deg(ll2.Lzarray(1:round(len/2), i)), 'b-');
    yy = ylabel(['$$\dot{\theta}_', num2str(i), '[\circ/s]$$'], 'Rotation', 0, 'Interpreter','latex');
    set(yy, 'Units', 'Normalized', 'Position', [-0.18, 0.35, 0]);
    grid on; box on;
    if i == n
        xlabel('Time [s]');
    end
    hold off;
end


