function zdot = pendulum_lagrange_4(t, z, p)
th1 = z(1); th2 = z(2); th3 = z(3); th4 = z(4); 
thd1 = z(5); thd2 = z(6); thd3 = z(7); thd4 = z(8); 
m1 = p.m(1); m2 = p.m(2); m3 = p.m(3); m4 = p.m(4); 
l1 = p.l(1); l2 = p.l(2); l3 = p.l(3); l4 = p.l(4); 
I1 = p.I(1); I2 = p.I(2); I3 = p.I(3); I4 = p.I(4); 
g = p.g;

A = zeros(4, 4);
b = zeros(4, 1);
A(1, 1) = I1 + (l1^2*m1)/4 + l1^2*m2 + l1^2*m3 + l1^2*m4;
A(1, 2) = (l1*l2*cos(th1 - th2)*(m2 + 2*m3 + 2*m4))/2;
A(1, 3) = (l1*l3*cos(th1 - th3)*(m3 + 2*m4))/2;
A(1, 4) = (l1*l4*m4*cos(th1 - th4))/2;
A(2, 1) = (l1*l2*cos(th1 - th2)*(m2 + 2*m3 + 2*m4))/2;
A(2, 2) = I2 + (l2^2*m2)/4 + l2^2*m3 + l2^2*m4;
A(2, 3) = (l2*l3*cos(th2 - th3)*(m3 + 2*m4))/2;
A(2, 4) = (l2*l4*m4*cos(th2 - th4))/2;
A(3, 1) = (l1*l3*cos(th1 - th3)*(m3 + 2*m4))/2;
A(3, 2) = (l2*l3*cos(th2 - th3)*(m3 + 2*m4))/2;
A(3, 3) = I3 + (l3^2*m3)/4 + l3^2*m4;
A(3, 4) = (l3*l4*m4*cos(th3 - th4))/2;
A(4, 1) = (l1*l4*m4*cos(th1 - th4))/2;
A(4, 2) = (l2*l4*m4*cos(th2 - th4))/2;
A(4, 3) = (l3*l4*m4*cos(th3 - th4))/2;
A(4, 4) = I4 + (l4^2*m4)/4;
b(1) = -(l1*(g*m1*sin(th1) + 2*g*m2*sin(th1) + 2*g*m3*sin(th1) + 2*g*m4*sin(th1) + l2*m2*thd2^2*sin(th1 - th2) + 2*l2*m3*thd2^2*sin(th1 - th2) + 2*l2*m4*thd2^2*sin(th1 - th2) + l3*m3*thd3^2*sin(th1 - th3) + 2*l3*m4*thd3^2*sin(th1 - th3) + l4*m4*thd4^2*sin(th1 - th4)))/2;
b(2) = -(l2*(g*m2*sin(th2) + 2*g*m3*sin(th2) + 2*g*m4*sin(th2) - l1*m2*thd1^2*sin(th1 - th2) - 2*l1*m3*thd1^2*sin(th1 - th2) - 2*l1*m4*thd1^2*sin(th1 - th2) + l3*m3*thd3^2*sin(th2 - th3) + 2*l3*m4*thd3^2*sin(th2 - th3) + l4*m4*thd4^2*sin(th2 - th4)))/2;
b(3) = (l3*(l1*m3*thd1^2*sin(th1 - th3) - 2*g*m4*sin(th3) - g*m3*sin(th3) + 2*l1*m4*thd1^2*sin(th1 - th3) + l2*m3*thd2^2*sin(th2 - th3) + 2*l2*m4*thd2^2*sin(th2 - th3) - l4*m4*thd4^2*sin(th3 - th4)))/2;
b(4) = (l4*m4*(l1*thd1^2*sin(th1 - th4) - g*sin(th4) + l2*thd2^2*sin(th2 - th4) + l3*thd3^2*sin(th3 - th4)))/2;
vec = A\b;

thdd1 = vec(1);
thdd2 = vec(2);
thdd3 = vec(3);
thdd4 = vec(4);

zdot = [thd1; thd2; thd3; thd4; thdd1; thdd2; thdd3; thdd4; ];
end
