function [theta] = Hexapod_pos_sol(x, y, z, alpha, beta, gamma)
    close all;
    theta = zeros(6,3);

    % give x,y,z in mm and alpha, beta, gamma in degrees
    leg_mounts = [135*cos(deg2rad(30)), 67.5, 0;
                  0, 80, 0;
                  135*cos(deg2rad(150)), 67.5, 0;
                  135*cos(deg2rad(210)), -67.5, 0;
                  0, -80, 0;
                  135*cos(deg2rad(330)), -67.5, 0]; %mm
   
    % Position of body in world frame
    body_pos_w = HTrans(x, y, z, alpha, beta, gamma);

    % Find the position of each leg in the world frame
    leg_mounts_h = [leg_mounts, ones(6,1)].';
    
    leg_world_h = body_pos_w * leg_mounts_h;   
    leg_mounts_w = leg_world_h(1:3,:).';

    % Foot Positions
    foot_pos_w = [135*cos(deg2rad(30)), 67.5, 0;
                  0, 90, 0;
                  135*cos(deg2rad(150)), 67.5, 0;
                  135*cos(deg2rad(210)), -67.5, 0;
                  0, -90, 0;
                  135*cos(deg2rad(330)), -67.5, 0];
    
    rel_pos_w = foot_pos_w - body_pos_w(1:3,4)';

    % Solve for the value of each joint for each leg
    for i=1:6
        theta_temp = Hexapod_IK(rel_pos_w(i,1), rel_pos_w(i,2), rel_pos_w(i,3), i);
        theta(i,1) = theta_temp(1,1);
        theta(i,2) = theta_temp(2,2);
        theta(i,3) = theta_temp(1,3);
        % theta(i,1) = 0;
        % theta(i,2) = 0;
        % theta(i,3) = 0;
    end
       
    % Plot Hexapod Position
    
    figure; hold on; axis equal; grid on; view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
   
    % Plot leg mounts and foot end positions
    plot3(leg_mounts_w(:,1), leg_mounts_w(:,2), leg_mounts_w(:,3), 'k-o', 'LineWidth',1.5);
    plot3(foot_pos_w(:,1),    foot_pos_w(:,2),    foot_pos_w(:,3),    'r.', 'MarkerSize',20);
    
    % For each leg compute FK and plot segments
    for i = 1:6
        jp = Hexapod_FK(leg_mounts_w(i,:), theta(i,:), i);
        
        % Plot legs
        plot3([jp(1,1), jp(2,1)], [jp(1,2), jp(2,2)], [jp(1,3), jp(2,3)], 'b-', 'LineWidth',2);
        plot3([jp(2,1), jp(3,1)], [jp(2,2), jp(3,2)], [jp(2,3), jp(3,3)], 'b-', 'LineWidth',2);
        plot3([jp(3,1), jp(4,1)], [jp(3,2), jp(4,2)], [jp(3,3), jp(4,3)], 'b-', 'LineWidth',2);
        
        % Plot joints
        plot3(jp(2:3,1), jp(2:3,2), jp(2:3,3), 'mo', 'MarkerFaceColor','m');
    end

    disp("Joint Values (Rotate, Hip, Knee): ")
    disp(rad2deg(theta))
end


function theta = Hexapod_IK(x,y,z, joint)
    % Hexapod Parameters
    t1 = deg2rad([30, 90, 150, 210, 270, 330]); % rad
    a1 = [135, 80, 135, 135, 80, 135]; % mm
    a2 = 40;  % mm
    a3 = 80;  % mm
    a4 = 125; % mm
      
    % Solve for theta 2
    p = [x; y];
    a1_vec = [cos(t1(joint)) -sin(t1(joint)); sin(t1(joint)) cos(t1(joint))]*[a1(joint); 0];
    cross = a1_vec(1)*y - a1_vec(2)*x; 
    A = norm(p - a1_vec);
    r = norm(p);
    cost2 = (r^2 - a1(joint)^2 - A^2)/(2*a1(joint)*A); 
    t2 = 2*atan2(sqrt(1-cost2^2), 1+cost2)*sign(cross);

    % Solve for theta 4
    s = norm([z, A-a2]);
    cost4 = (s^2-a3^2-a4^2)/(2*a3*a4);
    t4 = 2*atan2(sqrt(1-cost4^2), 1+cost4);
    t4 = [t4, -t4];

    % Solve for theta 3
    t3 = atan2(z, A-a2) - atan2(a4.*sin(t4), a3+a4.*cos(t4));

    theta = [t2, t3(1) , t4(1); t2, t3(2), t4(2)];
end


% Solves the FK given the joint mount, the joint angels and which leg it is
function joint_pos = Hexapod_FK(mount_pos, theta, joint)
    % Hexapod Parameters
    t1 = deg2rad([30, 90, 150, 210, 270, 330]); % rad
    a2 = 40;  % mm
    a3 = 80;  % mm
    a4 = 125; % mm

    % Find leg mount position
    T = eye(4);
    T = T * Trans(mount_pos(1), mount_pos(2), mount_pos(3));
    P0 = T(1:3,4).';
    T = T * RotZ(t1(joint));
    
    % Find first leg point
    T = T * RotZ(theta(1));
    T = T * Trans(a2, 0, 0);
    P2 = T(1:3,4).';

    % Rotate by alpha2
    T = T * RotY(pi/2);

    % Find second leg point
    T = T * RotY(theta(2));
    T = T * Trans(a3, 0, 0);
    P3 = T(1:3,4).';

    % Find third leg point
    T = T * RotY(theta(3));
    T = T * Trans(a4, 0, 0);
    P4 = T(1:3,4).';

    joint_pos = [P0; P2; P3; P4];
end


% Helper functions
function T = HTrans(x , y, z, a, b, c)
    % Given x,y,z yaw, pitch, roll, this given the homogenous
    % transfromation matrix T
    a = deg2rad(a);
    b = deg2rad(b);
    c = deg2rad(c);

    R = [cos(a)*cos(b), cos(a)*sin(b)*sin(c)-sin(a)*cos(c), cos(a)*sin(b)*cos(c)+sin(a)*sin(c);
         sin(a)*cos(b), sin(a)*sin(b)*sin(c)+cos(a)*cos(c), sin(a)*sin(b)*cos(c)-cos(a)*sin(c);
         -sin(b), cos(b)*sin(c), cos(b)*cos(c)]; %seems to have done yaw, pitch, roll (z,y,x) keep this in mind
    p = [x;y;z];
    T = [R, p; 
        0,0,0,1];
end

function T = Trans(x,y,z)
    T = eye(4);
    T(1:3,4) = [x;y;z];
end

function R = RotX(theta)
    R = [1, 0, 0, 0;
         0, cos(theta), -sin(theta), 0;
         0, sin(theta), cos(theta), 0;
         0, 0, 0, 1];
end

function R = RotY(theta)
    R = [cos(theta), 0, sin(theta), 0;
         0, 1, 0, 0;
         -sin(theta), 0, cos(theta), 0;
         0, 0, 0, 1];
end

function R = RotZ(theta)
    R = [cos(theta), -sin(theta), 0, 0;
         sin(theta), cos(theta), 0, 0;
         0, 0, 1, 0;
         0, 0, 0, 1];
end