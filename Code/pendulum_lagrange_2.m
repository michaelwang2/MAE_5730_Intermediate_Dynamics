function zdot = pendulum_lagrange_2(t, z, p)
th1 = z(1); th2 = z(2); 
thd1 = z(3); thd2 = z(4); 
m1 = p.m(1); m2 = p.m(2); 
l1 = p.l(1); l2 = p.l(2); 
I1 = p.I(1); I2 = p.I(2); 
g = p.g;

thdd1 = -(2*l1*(g*l2^2*m2^2*sin(th1) + 4*I2*g*m1*sin(th1) + 8*I2*g*m2*sin(th1) + g*l2^2*m2^2*sin(th1 - 2*th2) + l2^3*m2^2*thd2^2*sin(th1 - th2) + g*l2^2*m1*m2*sin(th1) + 4*I2*l2*m2*thd2^2*sin(th1 - th2) + l1*l2^2*m2^2*thd1^2*sin(2*th1 - 2*th2)))/(16*I1*I2 + 2*l1^2*l2^2*m2^2 + 4*I2*l1^2*m1 + 4*I1*l2^2*m2 + 16*I2*l1^2*m2 - 2*l1^2*l2^2*m2^2*cos(2*th1 - 2*th2) + l1^2*l2^2*m1*m2);
thdd2 = (2*l2*m2*(4*I1*l1*thd1^2*sin(th1 - th2) - 4*I1*g*sin(th2) + g*l1^2*m1*sin(2*th1 - th2) + 2*g*l1^2*m2*sin(2*th1 - th2) + l1^3*m1*thd1^2*sin(th1 - th2) + 4*l1^3*m2*thd1^2*sin(th1 - th2) - 2*g*l1^2*m2*sin(th2) + l1^2*l2*m2*thd2^2*sin(2*th1 - 2*th2)))/(16*I1*I2 + 2*l1^2*l2^2*m2^2 + 4*I2*l1^2*m1 + 4*I1*l2^2*m2 + 16*I2*l1^2*m2 - 2*l1^2*l2^2*m2^2*cos(2*th1 - 2*th2) + l1^2*l2^2*m1*m2);

zdot = [thd1; thd2; thdd1; thdd2; ];
end
