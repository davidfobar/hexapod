% rotY: Gives a rotation matrix for an angle about the y axis.
%
% R = rotY(theta)
% Returns a rotation matrix of the rotation about the y axis of given angle 
% provided in radians.
%
% R = the rotation matrix in the y direction of 
%
% theta = the angle of rotation in the y direction in radians
%
% Luke Boyd
% 10900074
% MEGN544
% 9/18/2025
function R = rotY(theta)
R = [ cos(theta), 0, sin(theta);
               0, 1          0 ;
     -sin(theta), 0, cos(theta)];
end