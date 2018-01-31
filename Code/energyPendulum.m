function [T, V, H, M] = energyPendulum(t, thetas, thetadots, p)
N = length(thetas(1, :));
lengths = p.l;
inertias = p.I;
g = p.g;
masses = p.m;
T = zeros(length(t), N);
V = zeros(length(t), N);
H = zeros(length(t), N);
M = zeros(length(t), N);
for i = 1:length(t)
    for a = 1:N
        % initialize position an velocity vectors
        r_G = zeros(3, 1);
        v_G_I = zeros(3, 1);
        d = 0;

        % find location of center of mass of current link in inertial
        % coordinates
        for b = 1:a
            % temporary vectors
            tempR_B = [lengths(b); 0; 0];
            tempV_B = [0; thetadots(i, b)*lengths(b); 0];

            % Rotation matrix from Link Frame to Inertial Frame
            R_B2I = [cos(thetas(i, b)), -sin(thetas(i, b)), 0; sin(thetas(i, b)), cos(thetas(i, b)), 0; 0, 0, 1];
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

        % KE and PE
        T(i, a) = 0.5*inertias(a)*(thetadots(i, a)^2) + 0.5*masses(a)*(v_G_I.' * v_G_I);
        V(i, a) = masses(a)*g*h;
        
        % Angular momentum about O
        tempH = [0; 0; inertias(a)*thetadots(i, a)] + masses(a).*cross(r_G, v_G_I);
        H(i, a) = tempH(3);
        
        % External torque about O
        tempM = cross(r_G, [masses(a)*g; 0; 0]);
        M(i, a) = tempM(3);
    end
end
end