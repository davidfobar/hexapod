clear; clc;
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
for i=1:nLegs
    legs(i).linkParams(1) = createLink(a1(i), 0, 0, t1(i),0,0,0,0,0);
    legs(i).linkParams(2) = createLink(a2, 0, pi/2, [],0,0,0,0,0);
    legs(i).linkParams(3) = createLink(a3, 0, 0, [],0,0,0,0,0);
    legs(i).linkParams(4) = createLink(a3, 0, 0, [],0,0,0,0,0);
end

T = eye(4);
T(1:3,4) = [30; -150; 0];
t_IK = Solve_IK(legs(5).linkParams, T);

transPercent = 0.20;
trajectory = [0.0   0 -180 -20;
              0.5 -50 -180 -20;
              1.0 -20 -220 20;
              1.5  20 -220 20;
              2.0  50 -180 -20;
              2.5   0 -180 -20];
ts = linspace(0, 2.5, 100);
ps = zeros(3,100);
vs = zeros(3,100);
as = zeros(3,100);
Q = zeros(3, 100, nLegs); % 3 joints per leg, 100 timesteps
for i=1:100
    [p,v,a] = constAccelInterp(ts(i), trajectory, transPercent);
    ps(:, i) = p;
    vs(:, i) = v;
    as(:, i) = a;

    [t2, t3, t4] = Solve_IK(legs(5).linkParams, p);
    Q(:,i,5) = [t2;t3;t4];
end

figure;
plot3(ps(1,:), ps(2,:), ps(3,:), 'LineWidth', 2);
grid on; axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('End-Effector Trajectory');

figure;
plot(ts, Q(1,:,5), 'LineWidth', 2); hold on;
plot(ts, Q(2,:,5), 'LineWidth', 2);
plot(ts, Q(3,:,5), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Joint Angle (rad)');
legend('Coxa','Femur','Tibia');
grid on;

generateVideo(Q, legs, 'title', 'test');