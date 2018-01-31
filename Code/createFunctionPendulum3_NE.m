function createFunctionPendulum3_NE()
% creates RHS function for 3-linked Pendulum using Newton Euler Method
% create symbols
syms th1 th2 th3 thd1 thd2 thd3 thdd1 thdd2 thdd3     real
syms g m1 m2 m3 I1 I2 I3 l1 l2 l3                     real

% unit vectors
i = [1; 0; 0]; j = [0; 1; 0]; k = [0; 0; 1];

% rotation vectors from link frame to inertial frame
R_B12I = [cos(th1), -sin(th1), 0;...
          sin(th1), cos(th1), 0;...
          0       , 0       , 1];
R_B22I = [cos(th2), -sin(th2), 0;...
          sin(th2), cos(th2), 0;...
          0       , 0       , 1];
R_B32I = [cos(th3), -sin(th3), 0;...
          sin(th3), cos(th3), 0;...
          0       , 0       , 1];

% center of mass positions in inertial coordinates (about O, E1, E2)
rG1_O = R_B12I*[l1/2; 0; 0];
rE1_O = R_B12I*[l1; 0; 0];
rG2_E1 = R_B22I*[l2/2; 0; 0];
rE2_E1 = R_B22I*[l2; 0; 0];
rG2_O = R_B12I*[l1; 0; 0] + R_B22I*[l2/2; 0; 0];
rG3_E2 = R_B32I*[l3/2; 0; 0];
rG3_E1 = R_B22I*[l2; 0; 0] + R_B32I*[l3/2; 0; 0];
rG3_O = R_B12I*[l1; 0; 0] + R_B22I*[l2; 0; 0] + R_B32I*[l3/2; 0; 0];

% angular velocities and angular accelerations
w_B1_I = [0; 0; thd1];
w_B2_I = [0; 0; thd2];
w_B3_I = [0; 0; thd3];
alpha_B1_I = [0; 0; thdd1];
alpha_B2_I = [0; 0; thdd2];
alpha_B3_I = [0; 0; thdd3];

% center of mass accelerations in inertial frame
a_G1_I = cross(alpha_B1_I, rG1_O) + cross(w_B1_I, cross(w_B1_I, rG1_O));
a_E1_I = cross(alpha_B1_I, rE1_O) + cross(w_B1_I, cross(w_B1_I, rE1_O));
a_G2_I = a_E1_I + cross(alpha_B2_I, rG2_E1) + cross(w_B2_I, cross(w_B2_I, rG2_E1));
a_E2_I = a_E1_I + cross(alpha_B2_I, rE2_E1) + cross(w_B2_I, cross(w_B2_I, rE2_E1));
a_G3_I = a_E2_I + cross(alpha_B3_I, rG3_E2) + cross(w_B3_I, cross(w_B3_I, rG3_E2));

% torque about O
M1_O = cross(rG1_O, [m1*g; 0; 0]);
M2_O = cross(rG2_O, [m2*g; 0; 0]);
M3_O = cross(rG3_O, [m3*g; 0; 0]);

% torque about E1
M2_E1 = cross(rG2_E1, [m2*g; 0; 0]);
M3_E1 = cross(rG3_E1, [m3*g; 0; 0]);

% torque about E2
M3_E2 = cross(rG3_E2, [m3*g; 0; 0]);

% angular momentum balance about O
Hd_O = m1.*cross(rG1_O, a_G1_I) + m2.*cross(rG2_O, a_G2_I) + m3.*cross(rG3_O, a_G3_I)...
    + I1.*alpha_B1_I + I2.*alpha_B2_I + I3.*alpha_B3_I;
eq1 = Hd_O - (M1_O + M2_O + M3_O); eq1 = eq1(3);

% angular momentum balance about E1
Hd_E1 = m2.*cross(rG2_E1, a_G2_I) + m3.*cross(rG3_E1, a_G3_I)...
    + I2.*alpha_B2_I + I3.*alpha_B3_I;
eq2 = Hd_E1 - (M2_E1 + M3_E1); eq2 = eq2(3);

% angular momentum balance about E2
Hd_E2 = m3.*cross(rG3_E2, a_G3_I)...
    + I3.*alpha_B3_I;
eq3 = Hd_E2 - M3_E2; eq3 = eq3(3);

% solve for thetadotdots
[r1, r2, r3] = solve(eq1, eq2, eq3, thdd1, thdd2, thdd3);

% write the file 
fid = fopen('pendulum_NE_3.m', 'w');
fprintf(fid, 'function zdot = pendulum_NE_3(t,z,p)\n');
fprintf(fid, 'm1 = p.m(1); m2 = p.m(2); m3 = p.m(3);\n');
fprintf(fid, 'I1 = p.I(1); I2 = p.I(2); I3 = p.I(3);\n');
fprintf(fid, 'l1 = p.l(1); l2 = p.l(2); l3 = p.l(3);\n');
fprintf(fid, 'g = p.g;\n');
fprintf(fid, 'th1 = z(1); th2 = z(2); th3 = z(3);\n');
fprintf(fid, 'thd1 = z(4); thd2 = z(5); thd3 = z(6);\n\n');
fprintf(fid, ['thdd1 = ', char(simplify(r1)), ';\n']);
fprintf(fid, ['thdd2 = ', char(simplify(r2)), ';\n']);
fprintf(fid, ['thdd3 = ', char(simplify(r3)), ';\n\n']);
fprintf(fid, 'zdot = [thd1; thd2; thd3; thdd1; thdd2; thdd3];\n');
fprintf(fid, 'end');
fclose(fid);
end