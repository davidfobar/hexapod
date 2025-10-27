clear; clc;
addpath("/home/dave/PhD-Classes/ROBO554/functions")

a1 = 135; a2 = 40; a3 = 80; a4 = 125;
t1 = pi/6;
linkParams = {t1, a1, a2, a3, a4};

jointLimits = [-pi/4     pi/4;    % t2
               -pi/3   2*pi/3;    % t3
               -5*pi/6   pi/10;]; % t4

%theta_i = [t1, jointLimits(1,1), jointLimits(2,1), jointLimits(3,1)];
%theta_f = [t1, jointLimits(1,2), jointLimits(2,2), jointLimits(3,2)];

theta_i = [t1, jointLimits(1,1), deg2rad(60), deg2rad(-120)];
theta_f = [t1, jointLimits(1,2), deg2rad(60), deg2rad(-120)];

ti = 0.0;
tf = 2.0;
Ti = Solve_FK(linkParams, theta_i);
Tf = Solve_FK(linkParams, theta_f);

nTimeSteps = 60;
tSteps = linspace(ti, tf, nTimeSteps);

% cubic interpolation of joint parameter space with zero velocity assumption
a = compute_Cubic_Interpolation_Coef_for_Zero_Velocity(theta_i, theta_f);
t = cubic_Interpolation_Eval(a, ti, tSteps, tf);
%generatePlot(t, tSteps, 'Parameter-space (θ,d) cubic interpolation', 'parameterSpace');

generateVideo(t, linkParams, 'Motion from limit to limit - Leg 1', 'lim2lim');
