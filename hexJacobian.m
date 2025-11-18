% hexJacobian: Computes the jacobian of one leg of a Hexpod.
%
% [Jhex] = hexJacobian(theta_list, a_list)
% Computes the 6×4 geometric Jacobian for a single 4-DOF leg, using
% standard DH kinematics. The Jacobian maps joint rates q̇ to the
%
% Jhex = Jacobian for a single 
%  - Top 3 rows (v) are linear velocity (m/s per rad/s)
%  - Bottom 3 rows (ω) are angular velocity (rad/s per rad/s)
%
%
% theta_list = List of theta parameters for the link robots 
% a_list = List of the link lengths
%
% Raul Medina-Estrada
% 10871796
% MEGN
% 11/7/25

function [Jhex] = hexJacobian(theta_list, a_list)   
    if size(theta_list) ~= 4
        error('theta_list is not 4 elements')
    end
    if size(a_list) ~= 4
        error('a_list is not 4 elements')
    end
    % theta 1 is set to a specific leg
    th1 = theta_list(1);
    % other thetas being chosen
    th2 = theta_list(2);
    th3 = theta_list(3);
    th4 = theta_list(4);
    
    % leg lengths
    a1 = a_list(1);
    a2 = a_list(2);
    a3 = a_list(3);
    a4 = a_list(4);

    %% Initialization
    % Creating Links
    link1 = createLink(a1, 0, 0, th1, 0, a1/2, 1, 1, 0);
    link2 = createLink(a2, 0, pi/2, th2, 0, a2/2, 1, 1, 1);
    link3 = createLink(a3, 0, 0, th3, 0, a3/2, 1, 1, 1);
    link4 = createLink(a4, 0, 0, th4, 0, a4/2, 1, 1, 1);
    
    % Transformation Matrices (Individual)
    T_01 = dhTransform(link1.a,link1.d,link1.alpha,link1.theta);
    T_12 = dhTransform(link2.a,link2.d,link2.alpha,link2.theta);
    T_23 = dhTransform(link3.a,link3.d,link3.alpha,link3.theta);
    T_34 = dhTransform(link4.a,link4.d,link4.alpha,link4.theta);
    
    % Transformation Matrices (Building up)
    T_02 = T_01*T_12;
    T_03 = T_01*T_12*T_23;
    
    % Transformation Matrices (Total)
    T_04 = T_01*T_12*T_23*T_34;
    T_14 = T_12*T_23*T_34;
    
    %% Creating Jacobians

    % z lists
    z0 = [0;0;1];
    z1 = T_01(1:3,3);
    z2 = T_02(1:3,3);
    z3 = T_03(1:3,3);
    
    % o list
    o0 = [0;0;0];
    o1 = T_01(1:3,4);
    o2 = T_02(1:3,4);
    o3 = T_03(1:3,4);
    o4 = T_04(1:3,4);
    
    
    Jhex = [cross(z0,(o4-o0)) cross(z1,(o4-o1)) cross(z2,(o4-o2)) cross(z3,(o4-o3));
        z0                  z1               z2                 z3];
end

% V = J*qdot;