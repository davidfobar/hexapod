function HT = interpolateHT(HT0, HT1, steps)

    R0 = HT0(1:3,1:3);
    R1 = HT1(1:3,1:3);

    p0 = HT0(1:3,4);
    p1 = HT1(1:3,4);

    q0 = rot2Quat(R0);
    q1 = rot2Quat(R1);
    q  = quat_slerp(q0, q1, steps);

    R = quat2Rot(q);
    p = (1-steps)*p0 + steps*p1;

    HT = [R p; 0 0 0 1];
end

function q = quat_slerp(q0, q1, t)

    dotp = dot(q0, q1);
    if dotp < 0
        q1 = -q1;
        dotp = -dotp;
    end

    if dotp > 0.9995
        q = (1-t)*q0 + t*q1;
        q = q / norm(q);
        return;
    end

    theta_0 = acos(dotp);
    sin_theta_0 = sin(theta_0);

    theta = theta_0 * t;
    sin_theta = sin(theta);

    s0 = cos(theta) - dotp * sin_theta / sin_theta_0;
    s1 = sin_theta / sin_theta_0;

    q = s0*q0 + s1*q1;
    q = q / norm(q);
end
