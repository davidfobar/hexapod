% Orientation of body with fixed feet

clear; clc; close all;

%% Setup Hexapod

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
    legs(j).linkParams(1) = createLink(a1(j), 0, 0, t1(j),0,0,0,0,0);
    legs(j).linkParams(2) = createLink(a2, 0, pi/2, [],0,0,0,0,0);
    legs(j).linkParams(3) = createLink(a3, 0, 0, [],0,0,0,0,0);
    legs(j).linkParams(4) = createLink(a3, 0, 0, [],0,0,0,0,0);
end
nActiveJoints = 3;

wayPointTimes = [0, 0.5, 1.0, 1.5, 2.0];
transPercent  = 0.20;
nStepCycles   = 3;
nTimeSteps    = 300;
totalTime     = wayPointTimes(end);
ts = linspace(0, nStepCycles*totalTime, nTimeSteps);

nPlotPoints = int16(nTimeSteps/nStepCycles);
ps = zeros(nActiveJoints,nPlotPoints);
vs = zeros(nActiveJoints,nPlotPoints);
as = zeros(nActiveJoints,nPlotPoints);
Q  = zeros(nActiveJoints, nTimeSteps, nLegs);

%% Define Orientation Parameters
bodyHeight       = 50; % mm
roll             = 20; % Degrees
pitch            = 0; % Degrees
yaw              = 0; % Degrees
bodyOrientation  = rpy2Rot(deg2rad(roll), deg2rad(pitch), deg2rad(yaw)); %RPY Euler Angle
bodyHT0          = [eye(3), [0;0;bodyHeight]; 0,0,0,1];
bodyHT1          = [bodyOrientation, [0; 0; bodyHeight+40]; 0,0,0,1];
footDistFromCoxa = 100;

%% Calculate Trajectories and Joint angles
for i=1:nLegs
    trajectory = buildOrientationTrajectory(wayPointTimes, bodyHT0, bodyHT1, ... 
                                    footDistFromCoxa, legs(i).linkParams);
    for j=1:nTimeSteps
        t = ts(j);
        if t > totalTime
            t = totalTime;
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

bodyHTs = zeros(4,4,nTimeSteps);

totalBodyTime = nStepCycles * totalTime;

for i = 1:nTimeSteps
    s = (ts(i) / totalBodyTime);  
    s = min(max(s,0),1);
    bodyHTs(:,:,i) = interpolateHT(bodyHT0, bodyHT1, s);
end

%% Plotting
generateOrientationVideo(Q, legs, bodyHTs)
Q = Q(:,1:nPlotPoints,5);
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
