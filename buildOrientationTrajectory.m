function trajectory = buildOrientationTrajectory(wayPointTimes, HT0, HT1, footDistFromCoxa, linkParams)
    wayPointTimes = wayPointTimes(:);
    if length(wayPointTimes) ~= 5
        error('This version of buildGaitTrajectory is hard-coded for 5 waypoints.');
    end    

    % extract link parameters and determine coxa position
    a1 = linkParams(1).a;
    t1 = linkParams(1).theta;

    coxaPos = [a1*cos(t1);
               a1*sin(t1);
               0];

    % initial foot position in body frame
    x_init    = coxaPos(1);
    y_init = coxaPos(2) + footDistFromCoxa*sign(coxaPos(2));
    z_init  = -HT0(3,4);

    p_body_init = [x_init; y_init; z_init];

    % compute the foot position in WORLD coordinates (it stays fixed in world)
    R0 = HT0(1:3,1:3);
    t0 = HT0(1:3,4);
    foot_world = R0 * p_body_init + t0;

    % prepare interpolation parameters (s goes from 0..1 across the waypoint set)
    s_vals = (wayPointTimes - wayPointTimes(1)) / (wayPointTimes(end) - wayPointTimes(1));
    N = numel(s_vals);

    % rotation endpoints and quaternions
    R1 = HT1(1:3,1:3);
    t1_vec = HT1(1:3,4);

    q0 = rotm2quat(R0);   % [w x y z]
    q1 = rotm2quat(R1);

    % Build waypoint positions in body frame by inverting interpolated body HT
    x_wps = zeros(N,1);
    y_wps = zeros(N,1);
    z_wps = zeros(N,1);

    for k=1:N
        s = s_vals(k);

        % quaternion slerp for rotation
        q_s = quat_slerp(q0, q1, s);
        R_s = quat2rotm(q_s);

        % linear interpolation of translation (acceptable for body translation)
        t_s = (1-s)*t0 + s*t1_vec;

        % invert transform: p_body = R_s' * (foot_world - t_s)
        p_body_s = R_s' * (foot_world - t_s);

        x_wps(k) = p_body_s(1);
        y_wps(k) = p_body_s(2);
        z_wps(k) = p_body_s(3);
    end

    trajectory = [wayPointTimes(:), x_wps, y_wps, z_wps];
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