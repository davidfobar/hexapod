% dhFwdKine - Compute the forward kinematics using DH parameters
%
% H = dhFwdKine(linkList, paramList)
% This function computes the homogeneous transformation matrix for a robotic manipulator
% using the DH parameters.
%
% H = the resulting homogeneous transformation matrix representing the pose of the end-effector
%     with respect to the base frame.
%
% linkList = an array of link structures created using createLink.m, containing DH parameters
% paramList = a vector of joint parameters (angles for rotary joints, displacements for prismatic joints)
%
% David Fobar
% CWID - 10950737
% ROBO554
% 03NOV2025
function H = dhFwdKine(linkList, paramList)
    numLinks = length(paramList);

    % initialize the overall transformation matrix as identity
    H = eye(4);

    % iterate through each link to compute the overall transformation
    for i = 1:numLinks
        if linkList(i).isRotary == 1     % rotatry
            h = dhTransform(linkList(i).a, linkList(i).d, linkList(i).alpha, paramList(i) - linkList(i).offset);
        elseif linkList(i).isRotary == 0 % prismatic
            h = dhTransform(linkList(i).a, paramList(i) - linkList(i).offset, linkList(i).alpha, linkList(i).theta);
        else                             % static
            h = dhTransform(linkList(i).a, linkList(i).d, linkList(i).alpha, linkList(i).theta);
        end

        % update the overall transformation matrix
        H = H*h;
    end
end