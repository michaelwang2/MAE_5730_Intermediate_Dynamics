function [KE, PE] = animatePendulumDAE(t, zarray, p, Title, z0, animation)
% z = [x1 y1 x1d y1d x2 y2 x2d y2d x3 y3 x3d y3d th1 th2 th3 thd1 thd2 thd3]
n = 3; % number of links
l = sum(p.l);

% positions and velocities
pos1 = zarray(:, 1:2);
vel1 = zarray(:, 3:4);
pos2 = zarray(:, 5:6);
vel2 = zarray(:, 7:8);
pos3 = zarray(:, 9:10);
vel3 = zarray(:, 11:12);
ths = zarray(:, 13:15);
thds = zarray(:, 16:18);
vels = zeros(length(zarray(:, 1)), 2, 3);
vels(:, :, 1) = vel1;
vels(:, :, 2) = vel2;
vels(:, :, 3) = vel3;
poss = zeros(length(zarray(:, 1)), 2, 3);
poss(:, :, 1) = pos1;
poss(:, :, 2) = pos2;
poss(:, :, 3) = pos3;

% calculate KE and PE
KE = zeros(length(zarray(:, 1)), n);
PE = zeros(length(zarray(:, 1)), n);
d = [p.l(1)/2, p.l(1) + p.l(2)/2, p.l(1) + p.l(2) + p.l(3)/2];
for i = 1:n
    h = d(i).*ones(length(zarray(:, 1)), 1) - poss(:, 1, i);
    KE(:, i) = 0.5*p.m(i).*(vels(:, 1, i).^2 + vels(:, 2, i).^2) + 0.5*p.I(i).*(thds(:, i).^2);
    PE(:, i) = p.m(i)*p.g.*h;
end
TT = sum(KE, 2);
VV = sum(PE, 2);

if animation == 'f'
    return;
end

% animate
ff = figure; hold on;
subplot(1,2,2);
hold on;
tt = plot(t, TT, 'g-');
v = plot(t, VV, 'b-');
tv = plot(t, TT+VV, 'k-');
grid on; box on; 
title('System Energy');
xlabel('Time [s]');
yy = ylabel('[J]','Rotation',0);
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
legend([tt,v,tv], 'KE', 'PE', 'Total', 'Location', 'Best');
hold off;

subplot(1,2,1); 
tic; hold on;
cur_time = toc;
title([Title, ' (g = ', num2str(p.g), ' m/s^2)']);
xlabel('X [m]');
yy = ylabel('Y [m]', 'Rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
plot(0, 0, 'b.', 'MarkerSize', 30);
if Title(6:10) == '4-Bar'
    ttt = z0(15);
    xxx = z0(9);
    yyy = z0(10);
    origin = plot(yyy + p.l(3)*sin(ttt)/2, -(xxx + p.l(3)*cos(ttt)/2), 'b.', 'MarkerSize', 30);
    fourth = plot([0, yyy + p.l(3)*sin(ttt)/2], [0, -(xxx + p.l(3)*cos(ttt)/2)], 'k-', 'LineWidth', 6);
end
while cur_time <= t(end)
    cur_time = toc;
    x1 = interp1(t, pos1(:, 1), cur_time);
    y1 = interp1(t, pos1(:, 2), cur_time);
    x2 = interp1(t, pos2(:, 1), cur_time);
    y2 = interp1(t, pos2(:, 2), cur_time);
    x3 = interp1(t, pos3(:, 1), cur_time);
    y3 = interp1(t, pos3(:, 2), cur_time);
    th1 = interp1(t, ths(:, 1), cur_time);
    th2 = interp1(t, ths(:, 2), cur_time);
    th3 = interp1(t, ths(:, 3), cur_time);
    ex = [x1, x2, x3] + [cos(th1)*p.l(1)/2, cos(th2)*p.l(2)/2, cos(th3)*p.l(3)/2];
    ey = [y1, y2, y3] + [sin(th1)*p.l(1)/2, sin(th2)*p.l(2)/2, sin(th3)*p.l(3)/2];
    
    pp = plot([0, ey], [0, -ex], 'k-', 'LineWidth', 6);
    axis equal; box on; grid on;
    if Title(6:10) == '4-Bar'
        cx = (yyy + p.l(3)*sin(ttt)/2)/2;
        cy = (-(xxx + p.l(3)*cos(ttt)/2))/2;
        dd = 1.2;
        axis([cx - l*dd/2, cx + l*dd/2, cy - l*dd/2, cy + l*dd/2]);
    else
        axis([-l l -l l]);
    end
    drawnow;
    
    if ~ishandle(ff)
        close all;
        break;
    end
    
    delete(pp);
end
if Title(6:10) == '4-Bar'
    delete(origin);
    delete(fourth);
end
end