% Robot Controls - State Feedback Control for RRbot Manipulator
% Author: Krutarth Trivedi | ktrivedi@wpi.edu

clear; 
clc; 
close all;

global A; global B; global K;

% define symbols
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot t1 t2 I1 I2 'real'
syms m1 r1 m2 l1 r2 g 'real' 'positive'

% physical parameters of the robot
m1_ = 1; r1_ = 0.45; l1_ = 1; I1_ = 0.084;
m2_ = 1; r2_ = 0.45; l2_ = 1; I2_ = 0.084;
g_ = 9.81;

% -------- Part A : Equilibrium points of the robot -------------------%

% Symbolic EOM derived in programming assignment 1
eom_1 = I1*theta1_ddot - t1 + (I2*(2*theta1_ddot + 2*theta2_ddot))/2 - (m2*(2*(l1*cos(theta1) + r2*cos(theta1 + theta2))*(r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)^2 + l1*sin(theta1)*theta1_dot^2 - l1*theta1_ddot*cos(theta1) - r2*cos(theta1 + theta2)*(theta1_ddot + theta2_ddot)) - 2*(l1*sin(theta1) + r2*sin(theta1 + theta2))*(l1*cos(theta1)*theta1_dot^2 + l1*theta1_ddot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_ddot + theta2_ddot) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)^2)))/2 + (m1*(2*r1^2*theta1_ddot*cos(theta1)^2 + 2*r1^2*theta1_ddot*sin(theta1)^2))/2 - g*m2*(l1*sin(theta1) + r2*sin(theta1 + theta2)) - g*m1*r1*sin(theta1);
eom_2 = (m2*(2*r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)*(r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot) + l1*cos(theta1)*theta1_dot) - 2*r2*cos(theta1 + theta2)*(l1*sin(theta1)*theta1_dot + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot))*(theta1_dot + theta2_dot)))/2 - (m2*(2*r2*cos(theta1 + theta2)*(r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)^2 + l1*sin(theta1)*theta1_dot^2 - l1*theta1_ddot*cos(theta1) - r2*cos(theta1 + theta2)*(theta1_ddot + theta2_ddot)) - 2*r2*sin(theta1 + theta2)*(l1*cos(theta1)*theta1_dot^2 + l1*theta1_ddot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_ddot + theta2_ddot) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)^2) + 2*r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)*(r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot) + l1*cos(theta1)*theta1_dot) - 2*r2*cos(theta1 + theta2)*(l1*sin(theta1)*theta1_dot + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot))*(theta1_dot + theta2_dot)))/2 - t2 + (I2*(2*theta1_ddot + 2*theta2_ddot))/2 - g*m2*r2*sin(theta1 + theta2);

% Substituting the values
eom_1 = subs(eom_1,[theta1_dot, theta2_dot, theta1_ddot, theta2_ddot, t1,t2],[0,0,0,0,0 ,0]);
eom_2 = subs(eom_2,[theta1_dot, theta2_dot, theta1_ddot, theta2_ddot, t1, t2],[0,0,0,0,0,0]);

sol = solve([eom_1 == 0, eom_2 == 0], [theta1, theta2]);

disp("Equilibirum points of the robot: ");
disp(" ");
disp("theta 1 = ");
disp(sol.theta1);
disp(" ");
disp("theta 2 = ");
disp(sol.theta2);
disp(" ");
disp("theta1_dot = 0");
disp(" ");
disp("theta2_dot = 0");
disp(" ");

% ----------- Part B : Jacobian Linearization -------------%

x = [theta1, theta2, theta1_dot, theta2_dot];

dX1 = theta1_dot;
dX2 = theta2_dot;
dX3 = (I2*t1 - I2*t2 + m2*r2^2*t1*cos(theta1 + theta2)^2 - m2*r2^2*t2*cos(theta1 + theta2)^2 + m2*r2^2*t1*sin(theta1 + theta2)^2 - m2*r2^2*t2*sin(theta1 + theta2)^2 + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot^2 + l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta2_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta2_dot^2 - l1*m2*r2*t2*sin(theta1)*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1)*cos(theta1 + theta2)^2 - l1*m2*r2*t2*cos(theta1)*cos(theta1 + theta2) - l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 + I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot^2 - I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 + I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta2_dot^2 - I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta2_dot^2 + g*m1*m2*r1*r2^2*sin(theta1)*cos(theta1 + theta2)^2 + g*m1*m2*r1*r2^2*sin(theta1)*sin(theta1 + theta2)^2 - g*l1*m2^2*r2^2*cos(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2) + 2*l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot*theta2_dot - 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot*theta2_dot + l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot^2 + l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta2_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta2_dot^2 - l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta1_dot^2 + l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta1_dot^2 + l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 + 2*l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot*theta2_dot - 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot*theta2_dot + 2*I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot*theta2_dot - 2*I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot*theta2_dot)/(I1*I2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*sin(theta1 + theta2)^2 - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2));
dX4 = (I1*t2 - I2*t1 + I2*t2 - m2*r2^2*t1*cos(theta1 + theta2)^2 + m2*r2^2*t2*cos(theta1 + theta2)^2 - m2*r2^2*t1*sin(theta1 + theta2)^2 + m2*r2^2*t2*sin(theta1 + theta2)^2 + l1^2*m2*t2*cos(theta1)^2 + m1*r1^2*t2*cos(theta1)^2 + l1^2*m2*t2*sin(theta1)^2 + m1*r1^2*t2*sin(theta1)^2 - I2*g*l1*m2*sin(theta1) - I2*g*m1*r1*sin(theta1) + I1*g*m2*r2*sin(theta1 + theta2) - l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot^2 - l1^3*m2^2*r2*cos(theta1)^3*sin(theta1 + theta2)*theta1_dot^2 + l1^3*m2^2*r2*sin(theta1)^3*cos(theta1 + theta2)*theta1_dot^2 - l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta2_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta2_dot^2 - l1*m2*r2*t1*sin(theta1)*sin(theta1 + theta2) + 2*l1*m2*r2*t2*sin(theta1)*sin(theta1 + theta2) - g*l1*m2^2*r2^2*sin(theta1)*cos(theta1 + theta2)^2 + g*l1^2*m2^2*r2*cos(theta1)^2*sin(theta1 + theta2) - l1*m2*r2*t1*cos(theta1)*cos(theta1 + theta2) + 2*l1*m2*r2*t2*cos(theta1)*cos(theta1 + theta2) + 2*l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 + l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta2_dot^2 - g*l1^2*m2^2*r2*cos(theta1)*sin(theta1)*cos(theta1 + theta2) - I1*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot^2 + I1*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 - I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot^2 + I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 - I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta2_dot^2 + I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta2_dot^2 - g*m1*m2*r1*r2^2*sin(theta1)*cos(theta1 + theta2)^2 + g*m1*m2*r1^2*r2*cos(theta1)^2*sin(theta1 + theta2) - g*m1*m2*r1*r2^2*sin(theta1)*sin(theta1 + theta2)^2 + g*m1*m2*r1^2*r2*sin(theta1)^2*sin(theta1 + theta2) + l1^3*m2^2*r2*cos(theta1)^2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 + g*l1*m2^2*r2^2*cos(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2) - l1^3*m2^2*r2*cos(theta1)*sin(theta1)^2*sin(theta1 + theta2)*theta1_dot^2 - 2*l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot*theta2_dot + 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot*theta2_dot - l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot^2 - l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta2_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta2_dot^2 + 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta1_dot^2 + l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta2_dot^2 - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta1_dot^2 - l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta2_dot^2 - 2*l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 - l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta2_dot^2 - 2*l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot*theta2_dot + 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot*theta2_dot + 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta1_dot*theta2_dot - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta1_dot*theta2_dot - g*l1*m1*m2*r1*r2*sin(theta1)^2*sin(theta1 + theta2) - 2*l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot*theta2_dot - l1*m1*m2*r1^2*r2*cos(theta1)^3*sin(theta1 + theta2)*theta1_dot^2 + l1*m1*m2*r1^2*r2*sin(theta1)^3*cos(theta1 + theta2)*theta1_dot^2 + 2*l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot*theta2_dot - 2*I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot*theta2_dot + 2*I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot*theta2_dot - g*l1*m1*m2*r1*r2*cos(theta1)*sin(theta1)*cos(theta1 + theta2) + l1*m1*m2*r1^2*r2*cos(theta1)^2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 - l1*m1*m2*r1^2*r2*cos(theta1)*sin(theta1)^2*sin(theta1 + theta2)*theta1_dot^2)/(I1*I2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*sin(theta1 + theta2)^2 - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2));

x_dot = [dX1, dX2, dX3, dX4];
u = [t1, t2];

% jacobian linearization
A = jacobian(x_dot, x);     %4x4
B = jacobian(x_dot, u);     %4x2

% substituting the physical parameters of the robot
A = subs(A,[m1, m2 ,l1, r1, r2, I1, I2, g, theta1_dot, theta2_dot], [m1_, m2_ ,l1_, r1_, r2_, I1_, I2_, g_, 0, 0]);
B = subs(B,[m1, m2 ,l1, r1, r2, I1, I2, g, theta1_dot, theta2_dot], [m1_, m2_ ,l1_, r1_, r2_, I1_, I2_, g_, 0, 0]);

disp("A: ");
disp(A);
disp("B: ");
disp(B);

% ------------ Part C : Linearization around equilibrium points --------%
fprintf("Stability around: (0,0,0,0) \n");
A1 = subs(A, [theta1, theta2], [sol.theta1(1), sol.theta2(1)]);
A1 = double(A1);
eig_A1 = eig(A1);
fprintf("Eigenvalues are: %.3f, %.3f, %.3f, %.3f \n",eig_A1(1), eig_A1(2),eig_A1(3),eig_A1(4));
numberOfPositiveRealParts = size(find(eig_A1 > 0), 1);
if (numberOfPositiveRealParts > 0)
    disp("The linearized system is unstable around the equilibirum point (0,0,0,0)");
elseif (find(all(eig_A1)))
    disp("The linearized system is marginally stable around the equilibirum point (0,0,0,0)");
else
    disp("The linearized system is asymptotically stable around the equilibirum point (0,0,0,0)");
end

disp(" ");
fprintf('Stability around: (pi,0,0,0) \n');
A2 = subs(A, [theta1, theta2], [sol.theta1(2), sol.theta2(2)]);
A2 = double(A2);
eig_A2 = eig(A2);
fprintf('Eigenvalues are: %.3f%+.3fj, %.3f%+.3fj, %.3f%+.3fj, %.3f%+.3fj \n',real(eig_A2(1)),imag(eig_A2(1)), real(eig_A2(2)),imag(eig_A2(2)),real(eig_A2(3)),imag(eig_A2(3)),real(eig_A2(4)),imag(eig_A2(4)));
numberOfPositiveRealParts = size(find(eig_A2 > 0), 1);
if (numberOfPositiveRealParts > 0)
    disp("The linearized system is unstable around the equilibirum point (pi,0,0,0)");
elseif (find(all(eig_A2)))
    disp("The linearized system is marginally stable around the equilibirum point (pi,0,0,0)");
else
    disp("The linearized system is asymptotically stable around the equilibirum point (pi,0,0,0)");
end

disp(" ");
fprintf('Stability around: (0,pi,0,0) \n');
A3 = subs(A, [theta1, theta2], [sol.theta1(3), sol.theta2(3)]);
A3 = double(A3);
eig_A3 = eig(A3);
fprintf('Eigenvalues are: %.3f%+.3fj, %.3f%+.3fj, %.3f%+.3fj, %.3f%+.3fj \n',real(eig_A3(1)),imag(eig_A3(1)), real(eig_A3(2)),imag(eig_A3(2)),real(eig_A3(3)),imag(eig_A3(3)),real(eig_A3(4)),imag(eig_A3(4)));
numberOfPositiveRealParts = size(find(eig_A3 > 0), 1);
if (numberOfPositiveRealParts > 0)
    disp("The linearized system is unstable around the equilibirum point (0,pi,0,0)");
elseif (find(all(eig_A3)))
    disp("The linearized system is marginally stable around the equilibirum point (0,pi,0,0)");
else
    disp("The linearized system is asymptotically stable around the equilibirum point (0,pi,0,0)");
end

disp(" ");

%-------- Investigating the controllability for "Upward" Configuration ---%

A = subs(A, [theta1, theta2], [sol.theta1(1), sol.theta2(1)]);
B = subs(B, [theta1, theta2], [sol.theta1(1), sol.theta2(1)]);

% converting symbols to double
A = double(A);
B = double(B);

% computing the controllability Matrix and check for controllability
C = ctrb(A,B);
rank_C = rank(C);
fprintf("Rank of Controllability Matrix: %0.0f \n",rank_C);

if rank_C == size(A,1)
    disp("As Rank of C equals to n, the linearized system is controllable for Upward Configuration");
else
    disp("As Rank of C is not equal to n, the linearized system is not controllable for Upward Configuration");
end

disp(" ");
% ---------- State-Feedback Controller ------------ %

% choosing the desired eigenvalues
desiredEigen = [-1.1, -2.2, -3.3, -4.4];
disp("Selected Eigenvalues are: [-1.1, -2.2, -3.3, -4.4]");
disp(" ");
% placing the desired eigenvalues to obtain the gain matrix
K = place(A,B,desiredEigen);
disp("Calculated Gain Matrix: ");
display(K);

% -------- Simulation and Plotting ---------%

% simulate the system for 10 sec for given control inputs using ODE45
T = 10;
y0 = [deg2rad(30), deg2rad(45), 0 ,0];
[t,y] = ode45(@ode_rrbot, [0,T], y0);

% reconstruct the applied control input
u = -K*y';
t1 = u(1,:);
t2 = u(2,:);

% plot the trajectories
figure('Name','Trajectories', 'NumberTitle','off');
subplot(3,2,1)
plot(t,(y(:,1)),'b');
title('theta1')
xlabel('T');
ylabel('rad');

subplot(3,2,2)
plot(t,(y(:,2)),'b')
title('theta2')
xlabel('T');
ylabel('rad');

subplot(3,2,3)
plot(t,(y(:,3)),'b')
title('theta1-dot')
xlabel('T');
ylabel('rad/s');

subplot(3,2,4);
plot(t,(y(:,4)),'b');
title('theta2-dot')
xlabel('T');
ylabel('rad/s');

subplot(3,2,5);
plot(t,t1,'b');
title('t1')
xlabel('s');
ylabel('Nm');

subplot(3,2,6);
plot(t,t2,'b');
title('t2')
xlabel('s');
ylabel('Nm');

fprintf('Eigenvalues are selected such that: \n\n');
fprintf('Torque1: %.3f < u1 < %.3f \n\n',min(t1),max(t1));
fprintf('Torque2: %.3f < u2 < %.3f \n\n', min(t2),max(t2));