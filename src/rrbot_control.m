% Robot Controls - State Feedback Control for RRbot Manipulator
% Author: Krutarth Trivedi | ktrivedi@wpi.edu

clear; close; clc;

% ROS Setup
rosinit;

% global variable to store Gain Matrix
global K;

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(30), deg2rad(45)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;

sample_theta1 = [];
sample_theta2 = [];
sample_theta1_dot = [];
sample_theta2_dot = [];
sample_time = [];
sample_tau1 = [];
sample_tau2 = [];

while(t < 10)

    t = toc;
    % read the joint states
    jointData = receive(JointStates);

    % implementing the state feedback controller  
    theta1 = jointData.Position(1,1);
    theta2 = jointData.Position(2,1);
    theta1_dot = jointData.Velocity(1,1);
    theta2_dot = jointData.Velocity(2,1);
    X = [theta1, theta2, theta1_dot, theta2_dot];
    u = -K*X';
    tau1.Data = u(1);
    tau2.Data = u(2);
    
    send(j1_effort,tau1);
    send(j2_effort,tau2);

    % sample the time, joint state values, and calculated torques here to be plotted at the end
    sample_time = [sample_time, t];
    sample_theta1 = [sample_theta1, theta1];
    sample_theta2 = [sample_theta2, theta2];
    sample_theta1_dot = [sample_theta1_dot, theta1_dot];
    sample_theta2_dot = [sample_theta2_dot, theta2_dot];
    sample_tau1 = [sample_tau1, u(1)];
    sample_tau2 = [sample_tau2, u(2)];
    
end

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

% disconnect from roscore
rosshutdown;

% plot the trajectories
figure('Name','State Trajectories', 'NumberTitle','off');
subplot(3,2,1)
plot(sample_time,sample_theta1,'b');
title('theta1')
xlabel('T');
ylabel('rad');

subplot(3,2,2)
plot(sample_time,sample_theta2,'b')
title('theta2')
xlabel('T');
ylabel('rad');

subplot(3,2,3)
plot(sample_time,sample_theta1_dot,'b')
title('theta1-dot')
xlabel('T');
ylabel('rad/s');

subplot(3,2,4);
plot(sample_time,sample_theta2_dot,'b');
title('theta2-dot')
xlabel('T');
ylabel('rad/s');

subplot(3,2,5);
plot(sample_time,sample_tau1,'b');
title('tau1')
xlabel('s');
ylabel('Nm');

subplot(3,2,6);
plot(sample_time,sample_tau2,'b');
title('tau2')
xlabel('s');
ylabel('Nm');