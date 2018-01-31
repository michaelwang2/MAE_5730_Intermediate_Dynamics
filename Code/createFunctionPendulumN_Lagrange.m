function createFunctionPendulumN_Lagrange(N)
% use Lagrange Equations on N-Linked Pendulum to create RHS 
% p.m = 1xN of masses of links
% p.I = 1xN of moment of inertias about Gi of links
% p.l = 1xN of lengths of links
% p.g = 1x1 gravitational constant

% use equationsToMatrix to extract thetadd's 
% check conservation of total energy and angular momentum for entire system
% add in unconservative forces to lagrange formulation
% make lengths small and sinusoidally excite end of n-linked pendulum.
% examine how frequency of driven function corresponds to modal shapes
% maybe insert natural spring constant (add to potential energy term) and
% see how it corresponds to 

% create symbols
syms g     real
thetas = sym('th', [1 N], 'real');
thetadots = sym('thd', [1 N], 'real');
thetadotdots = sym('thdd', [1 N], 'real');
masses = sym('m', [1 N], 'real');
lengths = sym('l', [1 N], 'real');
inertias = sym('I', [1 N], 'real');

% set up Lagrange equations one-by-one 
T = 0; V = 0;
for a = 1:N
    % Lagrange Equation for Link a
    
    % initialize position an velocity vectors
    r_G = zeros(3, 1);
    v_G_I = zeros(3, 1);
    d = 0;
    
    % find location of center of mass of current link in inertial
    % coordinates
    for b = 1:a
        % temporary vectors
        tempR_B = [lengths(b); 0; 0];
        tempV_B = [0; thetadots(b)*lengths(b); 0];
        
        % Rotation matrix from Link Frame to Inertial Frame
        R_B2I = [cos(thetas(b)), -sin(thetas(b)), 0; sin(thetas(b)), cos(thetas(b)), 0; 0, 0, 1];
        if b == a
            r_G = r_G + 0.5.*R_B2I*tempR_B;
            v_G_I = v_G_I + 0.5.*R_B2I*tempV_B;
            d = d + 0.5*lengths(b);
        else
            r_G = r_G + R_B2I*tempR_B;
            v_G_I = v_G_I + R_B2I*tempV_B;
            d = d + lengths(b);
        end
    end
    % height of the center of mass of the current link
    h = d - r_G(1); 
    
    % L = T - V
    T = T + 0.5*inertias(a)*(thetadots(a)^2) + 0.5*masses(a)*(v_G_I.' * v_G_I);
    V = V + masses(a)*g*h;
end
% Lagrange
L = T - V;

% set up lagrange equations
dL_dth = jacobian(L, thetas).';
dL_dthd = jacobian(L, thetadots).';
dL_dthd_dt = jacobian(dL_dthd, [thetas, thetadots]) * [thetadots.'; thetadotdots.'];
eqns = dL_dthd_dt - dL_dth;

% equationsToMatrix
A = simplify(jacobian(eqns, thetadotdots)); [r,c] = size(A);
b = simplify(-(eqns - A*(thetadotdots.')));

% create file 
fid = fopen(['pendulum_lagrange_', num2str(N), '.m'], 'w');
fprintf(fid, ['function zdot = pendulum_lagrange_', num2str(N), '(t, z, p)\n']);
% z = [thetas, thetadots]
for i = 1:N
    fprintf(fid, [char(thetas(i)), ' = z(', num2str(i), '); ']);
end
fprintf(fid, '\n');
for i = (N+1):2*N
    fprintf(fid, [char(thetadots(i-N)), ' = z(', num2str(i), '); ']);
end
fprintf(fid, '\n');
for i = 1:N
    fprintf(fid, [char(masses(i)), ' = p.m(', num2str(i), '); ']);
end
fprintf(fid, '\n');
for i = 1:N
    fprintf(fid, [char(lengths(i)), ' = p.l(', num2str(i), '); ']);
end
fprintf(fid, '\n');
for i = 1:N
    fprintf(fid, [char(inertias(i)), ' = p.I(', num2str(i), '); ']);
end
fprintf(fid, '\n');
fprintf(fid, 'g = p.g;\n\n');
fprintf(fid, ['A = zeros(', num2str(r), ', ', num2str(c), ');\n']);
fprintf(fid, ['b = zeros(', num2str(r), ', 1);\n']);
for i = 1:r
    for j = 1:c
        fprintf(fid, ['A(', num2str(i), ', ', num2str(j), ') = ', char(A(i, j)), ';\n']);
    end
end
for i = 1:r
    fprintf(fid, ['b(', num2str(i), ') = ', char(b(i)), ';\n']);
end
fprintf(fid, 'vec = A\\b;\n\n');
for i = 1:N
    fprintf(fid, [char(thetadotdots(i)), ' = vec(', num2str(i), ');\n']);
end
fprintf(fid, '\n');
fprintf(fid, 'zdot = [');
for i = 1:N
    fprintf(fid, [char(thetadots(i)), '; ']);
end
for i = 1:N
   fprintf(fid, [char(thetadotdots(i)), '; ']);
end
fprintf(fid, '];\n');
fprintf(fid, 'end\n');
fclose(fid);
end