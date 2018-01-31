function [Ex, Ey] = thetas2pos(n, t, p, zarray)
Ex = zeros(n, length(t)); Ex(1, :) = p.l(1).*cos(zarray(:, 1)');
Ey = zeros(n, length(t)); Ey(1, :) = p.l(1).*sin(zarray(:, 1)');
for i = 2:n
    Ex(i, :) = Ex(i-1, :) + p.l(i).*cos(zarray(:, i)');
    Ey(i, :) = Ey(i-1, :) + p.l(i).*sin(zarray(:, i)');
end
end