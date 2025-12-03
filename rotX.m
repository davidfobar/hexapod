% rotX: Gives a rotation matrix for an angle about the x axis.
%
% R = rotX(theta)
% Returns a rotation matrix of the rotation about the y axis of given angle 
% provided in radians.
%
% R = the rotation matrix in the x direction of 
%
% theta = the angle of rotation in the x direction in radians
%
% Luke Boyd
% 10900074
% MEGN544
% 9/18/2025
function R = rotX(theta)
R = [ 1,          0,           0;
      0, cos(theta), -sin(theta);
      0, sin(theta), cos(theta)];
end