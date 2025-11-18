% List of each of the a's (link lengths)
a_list = [1 1 1 1];

% Assumed 6 positions of each of the legs (all theta 1)
th1_list = [pi/4 pi/2 3*pi/2 -3*pi/2 -pi/2 -pi/4];

% Thetas 2-4 (Generalized to be the same for each leg)
th2_4list = zeros(1,3);

% List of each of the leg's thetas.
th_list1 = [th1_list(1) th2_4list];
th_list2 = [th1_list(2) th2_4list];
th_list3 = [th1_list(3) th2_4list];
th_list4 = [th1_list(4) th2_4list];
th_list5 = [th1_list(5) th2_4list];
th_list6 = [th1_list(6) th2_4list];

% Jacobians for each of the 6 legs at the same thetas and a's
J1 = hexJacobian(th_list1, a_list);
J2 = hexJacobian(th_list2, a_list);
J3 = hexJacobian(th_list3, a_list);
J4 = hexJacobian(th_list4, a_list);
J5 = hexJacobian(th_list5, a_list);
J6 = hexJacobian(th_list6, a_list);

% Rate of each of the joints (th1 is always zero since its rigid
q_rate = [0;1;1;1];

% Vecloity of the end effector = J*qdot (6x1)
V_ee = J1*q_rate;