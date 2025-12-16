clear; clc; close all;
% distance from body center to coxa pivot
a1 = [135; 80; 135; 135; 80; 135];
% angle from body x-axis to coxa pivot
t1 = [pi/6; pi/2; 5*pi/6; 7*pi/6; 3*pi/2; 11*pi/6];
% coxa, femur, tibia lengths
a2 = 40; a3 = 80; a4 = 125;

% set joint limits
jointLimits = [-pi/4     pi/4;    % t2
               -pi/3   2*pi/3;    % t3
               -5*pi/6   pi/10;]; % t4

nLegs = length(a1);
for j=1:nLegs
    legs(j).linkParams(1) = createLink(a1(j), 0,    0, t1(j), 0,     [0; 0; 0],     0, zeros(3,3), 0);
    legs(j).linkParams(2) = createLink(a2,    0, pi/2,    [], 0, [0.035; 0; 0], 0.072, zeros(3,3), 1);
    legs(j).linkParams(3) = createLink(a3,    0,    0,    [], 0, [0.040; 0; 0], 0.025, zeros(3,3), 1);
    legs(j).linkParams(4) = createLink(a3,    0,    0,    [], 0, [0.025; 0; 0], 0.080, zeros(3,3), 1);
end 
nActiveJoints = 3;

totalTime = 20;
wayPointTimes = linspace(0,totalTime,5);
transPercent  = 0.50;
nStepCycles   = 3;
nTimeSteps    = 300;
ts = linspace(0, nStepCycles*totalTime, nTimeSteps);

nPlotPoints = int16(nTimeSteps/nStepCycles);
ps = zeros(nActiveJoints,nPlotPoints);
vs = zeros(nActiveJoints,nPlotPoints);
as = zeros(nActiveJoints,nPlotPoints);
Q  = zeros(nActiveJoints, nTimeSteps, nLegs);
Qd  = zeros(nActiveJoints, nTimeSteps, nLegs);
Qdd  = zeros(nActiveJoints, nTimeSteps, nLegs);
torques  = zeros(nActiveJoints, nTimeSteps, nLegs);

% --- define tripod groups and phase offset ---
tripodA     = [1 3 5];
tripodB     = [2 4 6];
phaseOffset = totalTime/2;

% define gait parameters
bodyHeight       = 50;
footLiftHeight   = 40;
stepLength       = 100;
footSweepDist    = 40;
footDistFromCoxa = 100;

for i=1:nLegs
    trajectory = buildGaitTrajectory(wayPointTimes, bodyHeight, ...
                        footLiftHeight, stepLength, footSweepDist, ...
                        footDistFromCoxa, legs(i).linkParams);
    for j=1:nTimeSteps
        t = mod(ts(j), totalTime);

        if ismember(i, tripodB)
            t = t + phaseOffset;
            if t > totalTime
                t = t - totalTime;
            end
        end

        [p,v,a] = constAccelInterp(t, trajectory, transPercent);
        if j<= nPlotPoints && i == 5
            ps(:, j) = p;
            vs(:, j) = v;
            as(:, j) = a;
        end
    
        [t2, t3, t4] = Solve_IK(legs(i).linkParams, p);
        Q(:,j,i) = [t2;t3;t4];
    end
end

% Find the joint accelerations and velocities using time differences
for j = 1:length(ts)
    if j ~= length(ts)
        dt = ts(j+1) - ts(j);
        
        for i = 1:nLegs
            for k = 1:nActiveJoints
                Qd(k,:,i)  = gradient(Q(k,:,i), dt);
                Qdd(k,:,i) = gradient(Qd(k,:,i), dt);
            end
        end
    end
end

%% Dynamics of the walking

mass = 2; %kg, robot total mass
servoTorque = 10; %kgcm

% Create boundary conditions
bc.base_angular_velocity = [0;0;0];
bc.base_angular_acceleration = [0;0;0];
bc.base_linear_acceleration = [0;0; 9.81]; %Gravity compenstation
bc.distal_force = [0;0; mass/nLegs]; % force per leg is mass of robot divided by nLeg
bc.distal_torque = [0;0;0];


% Find joint torques for each leg
for i=1:nLegs
    for j = 1:nTimeSteps
        torques(:,j,i) = newtonEuler(legs(i).linkParams(2:4), Q(:,j,i), Qd(:,j,i), Qdd(:,j,i), bc);
    end
end

torques = torques / (100/9.81); %Convert to kg*cm



[speed, weight] = robotSpeedWeight(torques, totalTime, stepLength, mass, servoTorque); %speed in m/s weight in kg

disp('Maximum Speed (m/s) = ');
disp(speed)
disp('Maximum Extra Payload (kg) = ');
disp(weight)


%% Video and Plotting

% generateVideo(Q, legs);
Q = Q(:,1:nPlotPoints,5);
Qd = Qd(:,1:nPlotPoints,5);
Qdd = Qdd(:,1:nPlotPoints,5);
torques = torques(:,1:nPlotPoints,5);
ts = ts(1:nPlotPoints);


figure;
plot3(ps(1,:), ps(2,:), ps(3,:), 'LineWidth', 2);
grid on; axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('End-Effector Trajectory');

figure;
plot(ts, vs(1,:), 'LineWidth', 2); hold on;
plot(ts, vs(2,:), 'LineWidth', 2);
plot(ts, vs(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Joint Angle Velocity (rad/s)');
ylim([-110, 110]);
legend('Coxa','Femur','Tibia');
title('End-Effector Velocity for a Single Leg');
grid on;

figure;
plot(ts, as(1,:), 'LineWidth', 2); hold on;
plot(ts, as(2,:), 'LineWidth', 2);
plot(ts, as(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Joint Angle Acceleration (rad/s^2)');
ylim([-1100, 1100]);
legend('Coxa','Femur','Tibia');
title('End-Effector Acceleration for a Single Leg');
grid on;

figure;
plot(ts, Q(1,:), 'LineWidth', 2); hold on;
plot(ts, Q(2,:), 'LineWidth', 2);
plot(ts, Q(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Joint Angle (rad)');
legend('Coxa','Femur','Tibia');
title('Joint Actuator Angles for a Single Leg');
grid on;

figure;
plot(ts, Qd(1,:), 'LineWidth', 2); hold on;
plot(ts, Qd(2,:), 'LineWidth', 2);
plot(ts, Qd(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Joint Angle (rad/s)');
legend('Coxa','Femur','Tibia');
title('Joint Actuator Velocities for a Single Leg');
grid on;

figure;
plot(ts, Qdd(1,:), 'LineWidth', 2); hold on;
plot(ts, Qdd(2,:), 'LineWidth', 2);
plot(ts, Qdd(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Joint Angle (rad/s^2)');
legend('Coxa','Femur','Tibia');
title('Joint Actuator Accelerations for a Single Leg');
grid on;


figure;
plot(ts, torques(1,:), 'LineWidth', 2); hold on;
plot(ts, torques(2,:), 'LineWidth', 2);
plot(ts, torques(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Joint Angle (rad)');
legend('Coxa','Femur','Tibia');
title('Joint Actuator Torques for a Single Leg');
grid on


