% quat2Rot: Returns a rotation matrix given a quaternion.
%
% Q = rot2Quat(R): Returns the rotation matrix that corresponds to he
% quaternion, stacked [q0;q_vec]. 
% 
% R is the rotation matrix
% 
% Q is the quaternion
% 
% Luke Boyd
% 10900074
% MEGN544
% 9/18/2025

function R = quat2Rot(Q)
    q0 = Q(1); 
    q1 = Q(2); 
    q2 = Q(3); 
    q3 = Q(4);

    normQ = norm(Q);

    q0 = q0 / normQ; 
    q1 = q1 / normQ; 
    q2 = q2 / normQ; 
    q3 = q3 / normQ;

    % Products for ease of typing
    q0q0 = q0*q0; q1q1 = q1*q1; q2q2 = q2*q2; q3q3 = q3*q3;
    q0q1 = q0*q1; q0q2 = q0*q2; q0q3 = q0*q3;
    q1q2 = q1*q2; q1q3 = q1*q3; q2q3 = q2*q3;


    R = [ q0q0 + q1q1 - q2q2 - q3q3, 2*(q1q2 - q0q3), 2*(q1q3 + q0q2);
          2*(q1q2 + q0q3), q0q0 - q1q1 + q2q2 - q3q3, 2*(q2q3 - q0q1);
          2*(q1q3 - q0q2), 2*(q2q3 + q0q1), q0q0 - q1q1 - q2q2 + q3q3 ];
end