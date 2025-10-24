function t_IK = Solve_IK(linkParams, T)
    [t1, a1, a2, a3, a4] = deal(linkParams{:});

    % extract the end effector position
    p_0_4 = T(1:3,4); 
    px = p_0_4(1); py = p_0_4(2); pz = p_0_4(3);
    
    % Solve for theta 2
    p      = [px; py];
    a1_vec = [cos(t1) -sin(t1); sin(t1) cos(t1)]*[a1; 0];
    cross  = a1_vec(1)*p(2) - a1_vec(2)*p(1); 
    A      = norm(p - a1_vec);
    r      = norm(p);
    cos_t2 = (r^2-a1^2-A^2)/(2*a1*A);
    t2     = 2*atan2(sqrt(1-cos_t2^2), 1+cos_t2)*sign(cross);

    % Solve for theta 4
    s      = norm([pz, A-a2]);
    cos_t4 = (s^2-a3^2-a4^2)/(2*a3*a4);
    t4     = 2*atan2(sqrt(1-cos_t4^2), 1+cos_t4);
    t4     = [t4, -t4];

    % Solve for theta 3
    t3     = atan2(pz, A-a2) - atan2(a4.*sin(t4), a3+a4.*cos(t4));

    t_IK = [t1 t2 t3(1) t4(1);
            t1 t2 t3(2) t4(2);];
end