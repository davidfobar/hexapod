function [t2, t3, t4] = Solve_IK(linkList, P, th_last)
    % extract link parameters
    t1 = linkList(1).theta;
    a1 = linkList(1).a;
    a2 = linkList(2).a;
    a3 = linkList(3).a;
    a4 = linkList(4).a;

    % extract the end effector position
    px = P(1); py = P(2); pz = P(3);
    
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

    % if th_last is provided, only return the results for that solution
    if exist('th_last', 'var')
        errs = vecnorm(th_last - t_IK);
        [~, idx]  = min(errs);
        t3 = t3(idx);
        t4 = t4(idx);
    else
        t3 = t3(2);
        t4 = t4(2);
    end
end