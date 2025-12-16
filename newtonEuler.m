% NEWTONEULER computes the inverse dynamics of a serial link mainpulator
%
% [jointTorques] = newtonEuler( linkList, paramList,
%                      paramListDot, paramListDDot,  boundry_conditions )
%
% Outputs:
% jointTorques - The calculated torques of each joint.
%
% Inputs:
% linkList - list of the joint parameters created by createLink
% 
% paramList - the joint angles/distances of the manipulator
% 
% paramListDot - the joint angular/linear velocities of the manipulator
% 
% paramListDDot - the joint angular\linear accelerations of the manipulator
% 
% boundry_conditions - structure containing:
%   base_angular_velocity
%   base_angular_acceleration
%   base_linear_acceleration (add gravity in here)
%   distal_force (in tool frame)
%   distal_torque (in tool frame)
%
% Luke Boyd
% 10900074
% MEGN544
% 12/7/2025


function [jointTorques] = newtonEuler(linkList, paramList, paramListDot, paramListDDot, boundry_conditions)

    % Number of Joints
    numJoints = length(linkList);
    
    z_list  = zeros(3,numJoints); % Z vector
    T_list  = zeros(4,4,numJoints); % Homogenous Transformation
    d_list  = zeros(3,numJoints); % Displacment Vector 
    vd_list = zeros(3,numJoints); % Acceleration Vector
    w_list  = zeros(3,numJoints); % Angular Velocity Vector
    wd_list = zeros(3,numJoints); % Angular Acceleration Vector
    f_list  = zeros(3,numJoints); % Force Vector
    n_list  = zeros(3,numJoints); % Internal Torque Vector
    
    
    
    % Initialize link variables that get propagated forward
    T = eye(4); % Transform from 0 to joint i
    w_list(:, 1)  = boundry_conditions.base_angular_velocity;     % Angluar Velocity in joint frame
    wd_list(:, 1) = boundry_conditions.base_angular_acceleration; % Angular Acceleration in joint frame
    vd_list(:, 1) = boundry_conditions.base_linear_acceleration;  % Linear acceleration in joint frame
    
    num_static = 0; % number of static links encountered
    
    for i=1:numJoints % begin forward iteration from base to tool
        
        % Calculate link transform from i-1 to i
        if linkList(i).isRotary == 1
            % hint i-num_static is the param index to be on
            Ti = dhTransform(linkList(i).a, linkList(i).d, linkList(i).alpha, paramList(i-num_static) - linkList(i).offset);
        elseif linkList(i).isRotary == 0
            Ti = dhTransform(linkList(i).a, paramList(i-num_static) - linkList(i).offset, linkList(i).alpha, linkList(i).theta);
        else
            Ti = dhTransform(linkList(i).a, linkList(i).d, linkList(i).alpha, linkList(i).theta);
            num_static = num_static+1;
        end
    
        % Extract values
        w_prev  = w_list(:, i);
        wd_prev = wd_list(:, i);
        vd_prev = vd_list(:, i);
        z_prev  = T(1:3,3);
        di = T(1:3,1:3)*Ti(1:3,4); 
        
        % Update joint frame acceleartion, angular Acceleration, and
        % angular velocity.  In the i-i frame (so z is [0;0;1])
        if linkList(i).isRotary == 1  
            w = w_prev + paramListDot(i)*z_prev; % update ang vel
            wd = wd_prev + paramListDDot(i)*z_prev + paramListDot(i)*cross(w_prev, z_prev); % update ang accel in joint frame
        else
            w = w_prev;
            wd = wd_prev;
        end
    
        vd = vd_prev + cross(wd, di) + cross(w, cross(w, di));
        if linkList(i).isRotary == 0 % Prismatic
            vd = vd + 2*paramListDot(i)*cross(w_prev, z_prev) + paramListDDot(i)*z_prev; % update accel in joint frame
        end
        
        % Update lists
        w_list(:,i+1) = w;
        wd_list(:,i+1) = wd;
        vd_list(:,i+1) = vd;
    
        T = T*Ti;
        T_list(:,:,i) = T;
        d_list(:,i) = T(1:3, 4);
        R = T(1:3, 1:3);
        z_list(:,i) = z_prev;
        
        % Calculate the displacement from the i-1 frame to the i'th com in
        % the i'th frame
        ri1_i =  R*linkList(i).centOfMass;
        
        % Calculate the Acceleration of the Center of Mass of the link in
        % the ith frame
        vdCOM = vd + cross(wd, ri1_i) + cross(w, cross(w, ri1_i));
        
        % Calculate and Save Inertial Force and Torque in the i'th frame
        m = linkList(i).mass;
        I = linkList(i).inertia;
        J = R*I*R'; % Rotate into world frame
        f_list(:, i) = m * vdCOM; % Newton's Equation
        n_list(:, i) = J * wd + cross(w, J*w); % Euler's Equation    
        
    end % End Forward Iterations
    
    % Initialize variables for force/torque propagation
    f = boundry_conditions.distal_force; % Initialize force to external force on the tool in the tool frame
    n = boundry_conditions.distal_torque; % Initialize torque to external torque on the tool in the tool frame
    
    % Rotate f & n to base frame
    f = R * f;
    n = R * n;
    
    % preallocate joint torque vector
    jointTorques = zeros(numJoints,1); % preallocate for speed
    
    for i = numJoints:-1:1 % From Last joint to base
        
        % displacement from origin i-1 to i in base. Hint: use list(i).doi to help...
        if( i> 1)
            d = d_list(:,i) - d_list(:, i-1);
        else
            d = d_list(:,i);
        end

        r_i = T_list(1:3, 1:3, i) * linkList(i).centOfMass;
        
        % Update Force on joint i in base frame with inertial force from before
        f = f + f_list(:, i);

        % Update Torque all at once with f
        n = n + n_list(:, i) + cross(r_i, f_list(:, i)) + cross(d, f);
            
        if linkList(i).isRotary == 1 % Rotational Joint
            % joint i torque is the Z component
            jointTorques(i-num_static) = z_list(:, i)' * n;
            
        elseif linkList(i).isRotary == 0 % Prysmatic
            % joint i force is the Z component
            jointTorques(i-num_static) = z_list(:, i)' * f;
            
        else
            num_static = num_static-1;
        end
        
    end % End Backword Iterations
end % end newtonEular Function definition