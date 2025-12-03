clear; clc;
addpath("/home/dave/PhD-Classes/ROBO554/functions")
syms a1 a2 a3 a4 t1 t2 t3 t4
linkParams = {t1, a1, a2, a3, a4};


leg.linkParams(1) = createLink(a1, 0, 0, t1,0,0,0,0,0);
leg.linkParams(2) = createLink(a2, 0, pi/2, t2,0,0,0,0,0);
leg.linkParams(3) = createLink(a3, 0, 0, t3,0,0,0,0,0);
leg.linkParams(4) = createLink(a3, 0, 0, t4,0,0,0,0,0);

H = dhFwdKine(leg.linkParams, [0,0,0,0]);
H = simplify(H, 'IgnoreAnalyticConstraints', true);
H = vpa(H, 3);
disp(H);