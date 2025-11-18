% dhTransform: Returns the homogeneous transform corresponding to the provided DH parameters.
%
% H = dhTransform(a, d, alpha, theta)
% Computes the homogeneous transformation matrix using the Denavit-Hartenberg parameters.
%
% H = 4x4 homogeneous transformation matrix
% a = Link length (distance along the x-axis)
% d = Link offset (distance along the z-axis)
% alpha = Link twist (rotation around the x-axis)
% theta = Joint angle (rotation around the z-axis)
%
% Raul Medina-Estrada
% 10871796
% MEGN 544
% 10/6/25

function H = dhTransform(a, d, alpha, theta)
    % Compute the homogeneous transformation matrix using DH parameters
    H = [cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
         sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
         0, sin(alpha), cos(alpha), d;
         0, 0, 0, 1];
end
