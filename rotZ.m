% rotZ: Gives a rotation matrix for an angle about the z axis.
%
% R = rotZ(theta)
% Returns a rotation matrix of the rotation about the z axis of given angle 
% provided in radians.
%
% R = the rotation matrix in the z direction of 
%
% theta = the angle of rotation in the z direction in radians
%
% Luke Boyd
% 10900074
% MEGN544
% 9/18/2025
function R = rotZ(theta)
R = [ cos(theta), -sin(theta), 0;
      sin(theta),  cos(theta), 0;
               0,           0, 1];
end