function createFunction4BarLinkage()
% create the RHS function file of a 4-bar linkage using the DAE method
% create symbols
syms x1 y1 x1d y1d x1dd y1dd                          real
syms x2 y2 x2d y2d x2dd y2dd                          real
syms x3 y3 x3d y3d x3dd y3dd                          real
syms xend yend                                        real
syms th1 th2 th3 thd1 thd2 thd3 thdd1 thdd2 thdd3     real
syms g m1 m2 m3 I1 I2 I3 l1 l2 l3                     real
syms R0x R0y R1x R1y R2x R2y R3x R3y                  real

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

% position vectors
r_G1_O = [x1; y1; 0];
r_O1_O = r_G1_O + R_B12I*[l1/2; 0; 0];
r_G2_O = [x2; y2; 0];
r_G2_O1 = r_G2_O - r_O1_O;
r_O2_O1 = r_G2_O + R_B22I*[l2/2; 0; 0];
r_G3_O = [x3; y3; 0];
r_G3_O1 = r_G3_O - r_O1_O;
r_G3_O2 = r_G3_O - r_O2_O1 - r_O1_O;

% Link 1
F1x = R0x + m1*g - R1x;
F1y = R0y - R1y;
M1 = cross(R_B12I*[l1/2; 0; 0], [-R1x; - R1y; 0]) + cross(R_B12I*[-l1/2; 0; 0], [R0x;  R0y; 0]);
M1 = M1(3);

% Link 2
F2x = R1x + m2*g - R2x;
F2y = R1y - R2y;
M2 = cross(R_B22I*[-l2/2; 0; 0], [R1x; R1y; 0]) + cross(R_B22I*[l2/2; 0; 0], [-R2x; -R2y; 0]);
M2 = M2(3);

% Link 3
F3x = R2x + m3*g - R3x;
F3y = R2y - R3y;
M3 = cross(R_B32I*[-l3/2; 0; 0], [R2x; R2y; 0]) + cross(R_B32I*[l3/2; 0; 0], [-R3x; -R3y; 0]);
M3 = M3(3);

% build equations
eqs(1) =  F1x./m1 - x1dd;
eqs(2) = F1y./m1 - y1dd;
eqs(3) =  F2x./m2 - x2dd;
eqs(4) = F2y./m2 - y2dd;
eqs(5) =  F3x./m3 - x3dd;
eqs(6) = F3y./m3 - y3dd;
eqs(7) = M1./I1 - thdd1;
eqs(8) = M2./I2 - thdd2;
eqs(9) = M3./I3 - thdd3;

% Link 1 Constraint
eqs(10) = x1dd - l1/2 * (-thdd1*sin(th1) - (thd1^2)*cos(th1));
eqs(11) = y1dd - l1/2 * (thdd1*cos(th1) - (thd1^2)*sin(th1));

% Link2 Constraint
eqs(12) = x1dd + l1/2 * (-thdd1*sin(th1) - (thd1^2)*cos(th1)) - (x2dd - l2/2 * (-thdd2*sin(th2) - (thd2^2)*cos(th2)));
eqs(13) = y1dd + l1/2 * (thdd1*cos(th1) - (thd1^2)*sin(th1)) - (y2dd - l2/2 * (thdd2*cos(th2) - (thd2^2)*sin(th2)));

% Link 3 Constraint
eqs(14) = x2dd + l2/2 * (-thdd2*sin(th2) - (thd2^2)*cos(th2)) - (x3dd - l3/2 * (-thdd3*sin(th3) - (thd3^2)*cos(th3)));
eqs(15) = y2dd + l2/2 * (thdd2*cos(th2) - (thd2^2)*sin(th2)) - (y3dd - l3/2 * (thdd3*cos(th3) - (thd3^2)*sin(th3)));
eqs(16) = x3dd + l3/2 * (-thdd3*sin(th3) - (thd3^2)*cos(th3));
eqs(17) = y3dd + l3/2 * (thdd3*cos(th3) - (thd3^2)*sin(th3));

% form equations
vars = [x1dd; y1dd; x2dd; y2dd; x3dd; y3dd; thdd1; thdd2; thdd3; R0x; R0y; R1x; R1y; R2x; R2y; R3x; R3y];
A = jacobian(eqs.', vars);
b = -(eqs.' - A*vars);
vec = simplify(A\b);

% create MATLAB RHS file
% z = [x1 y1 x1d y1d x2 y2 x2d y2d x3 y3 x3d y3d th1 th2 th3 thd1 thd2 thd3]
fid = fopen('FourBarLinkage_DAE.m', 'w');
fprintf(fid, 'function zdot = FourBarLinkage_DAE(t,z,p)\n');
fprintf(fid, 'x1 = z(1); y1 = z(2); x1d = z(3); y1d = z(4);\n');
fprintf(fid, 'x2 = z(5); y2 = z(6); x2d = z(7); y2d = z(8);\n');
fprintf(fid, 'x3 = z(9); y3 = z(10); x3d = z(11); y3d = z(12);\n');
fprintf(fid, 'th1 = z(13); th2 = z(14); th3 = z(15);\n');
fprintf(fid, 'thd1 = z(16); thd2 = z(17); thd3 = z(18);\n');
fprintf(fid, 'l1 = p.l(1); l2 = p.l(2); l3 = p.l(3);\n');
fprintf(fid, 'I1 = p.I(1); I2 = p.I(2); I3 = p.I(3);\n');
fprintf(fid, 'm1 = p.m(1); m2 = p.m(2); m3 = p.m(3);\n');
fprintf(fid, 'g = p.g;\n\n');
for i = 1:9
    fprintf(fid, [char(vars(i)), ' = ', char(vec(i)), ';\n']);
end
fprintf(fid, '\nzdot = [x1d; y1d; x1dd; y1dd; x2d; y2d; x2dd; y2dd; x3d; y3d; x3dd; y3dd; thd1; thd2; thd3; thdd1; thdd2; thdd3];\n');
fprintf(fid, 'end\n');
fclose(fid);
end