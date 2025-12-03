% rot2Quat: Returns a quaternion given an rotation matrix.
%
% Q = rot2Quat(R): Returns the quaternion [qo;q_vec] that corresponds to
% the rotation matrix.
% 
% Q is the quaternion
% 
% R is the rotation matrix
%
% Luke Boyd
% 10900074
% MEGN544
% 9/18/2025

function Q = rot2Quat(R)
    
    q0sq = 0.25*(1+R(1,1)+R(2,2)+R(3,3));
    q1sq = 0.25*(1+R(1,1)-R(2,2)-R(3,3));
    q2sq = 0.25*(1-R(1,1)+R(2,2)-R(3,3));
    q3sq = 0.25*(1-R(1,1)-R(2,2)+R(3,3));
    q = [q0sq, q1sq, q2sq, q3sq];

    [~, index] = max(q);

    if(index == 1)
        q0 = sqrt(max(q0sq,0));
        q1 = (R(3,2)-R(2,3)) / (4*q0);
        q2 = (R(1,3)-R(3,1)) / (4*q0);
        q3 = (R(2,1)-R(1,2)) / (4*q0);

    elseif(index == 2)
        q1 = sqrt(max(q1sq,0));
        q0 = (R(3,2)-R(2,3)) / (4*q1);
        q2 = (R(1,2)+R(2,1)) / (4*q1);
        q3 = (R(3,1)+R(1,3)) / (4*q1);

    elseif(index == 3)
        q2 = sqrt(max(q2sq,0));
        q0 = (R(1,3)-R(3,1)) / (4*q2);
        q1 = (R(1,2)+R(2,1)) / (4*q2);
        q3 = (R(2,3)+R(3,2)) / (4*q2);

    else
        q3 = sqrt(max(q3sq,0));
        q0 = (R(2,1)-R(1,2)) / (4*q3);
        q1 = (R(3,1)+R(1,3)) / (4*q3);
        q2 = (R(2,3)+R(3,2)) / (4*q3);

    end
    Q = [q0; q1 ; q2; q3];

end