function animatePendulum(t, Ex, Ey, zarray, p, method, energy)
N = length(Ex(:,1));

% conservation of energy
if energy == 't'
    [T, V, H, M] = energyPendulum(t, zarray(:, 1:N), zarray(:, (N+1):end), p);
    TT = sum(T, 2);
    VV = sum(V, 2);
end
l = sum(p.l);
g = p.g; 

ff = figure; hold on;
if energy == 't'
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
end

tic; hold on;
cur_time = toc;
title([method, ': ', num2str(N), '-Linked Pendulum (g = ', num2str(g), ' m/s^2)']);
xlabel('X [m]');
yy = ylabel('Y [m]', 'Rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
plot(0, 0, 'b.', 'MarkerSize', 30);
while cur_time <= t(end)
    cur_time = toc;
    ex = interp1(t, Ex', cur_time);
    ey = interp1(t, Ey', cur_time);
    
    pp = plot([0, ey], [0, -ex], 'k-', 'LineWidth', 6);
    axis equal; box on; grid on;
    axis([-l l -l l]);
    drawnow;
    
    if ~ishandle(ff)
        close all;
        break;
    end
    
    delete(pp);
end
end