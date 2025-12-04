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

bodyPos = zeros(3,nTimeSteps);
bodyRPY = zeros(3,nTimeSteps);

%raise height from 20 to 50 and gradually roll during motion
z_mean = 50;
z_amp  = 5;
for i = 1:nTimeSteps
    t_local = mod(ts(i), totalTime);
    bodyPos(3,i) = z_mean + z_amp * sin(2*pi*(t_local/totalTime));
    bodyRPY(:,i) = [5*pi/180 * sin(2*pi*ts(i)/totalTime); 0;0;];
end


nPlotPoints = int16(nTimeSteps/nStepCycles);
ps = zeros(nActiveJoints,nPlotPoints);
vs = zeros(nActiveJoints,nPlotPoints);
as = zeros(nActiveJoints,nPlotPoints);
Q  = zeros(nActiveJoints, nTimeSteps, nLegs);

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
    trajectory = buildGaitTrajectory(wayPointTimes, ...
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

        [p_world,v,a] = constAccelInterp(t, trajectory, transPercent);

        % Calculate the rotation matrix based on the body roll, pitch, and yaw
        roll = bodyRPY(1,j);
        pitch = bodyRPY(2,j);
        yaw = bodyRPY(3,j);
        R_WB = rotZ(yaw)*rotY(pitch)*rotX(roll);
        p_WB = bodyPos(:,j);

        % move the trajectory to the world frame
        p_body = R_WB' * (p_world' - p_WB);

        if j<= nPlotPoints && i == 5
            ps(:, j) = p_world;
            vs(:, j) = v;
            as(:, j) = a;
        end
    
        [t2, t3, t4] = Solve_IK(legs(i).linkParams, p_body);
        Q(:,j,i) = [t2;t3;t4];
    end
end

generateVideo(Q, legs, bodyRPY, bodyPos);
Q = Q(:,1:nPlotPoints,5);
ts = ts(1:nPlotPoints);


%%
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

