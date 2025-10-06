% Homework 5 - 
clear all; clc;
addpath('../../functions');

function t_IK = Hexapod_IK(linkParams, H_0_4_test)
    [t1, a1, a2, a3, a4] = deal(linkParams{:});

    % extract the end effector position
    p_0_4 = H_0_4_test(1:3,4); 
    px = p_0_4(1); py = p_0_4(2); pz = p_0_4(3);
    
    % Solve for theta 2
    p      = [px;py];
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

N_tests = 100;

% Link parameters
a1 = 135; a2 = 40; a3 = 80; a4 = 125;
t1 = pi/6;

linkParams = {t1, a1, a2, a3, a4};

% Joint limits [deg]
JL_deg = [-pi/4     pi/4;    % t2
          -pi/3   2*pi/3;    % t3
          -5*pi/6   pi/10;]; % t4
JL_min = JL_deg(:,1).';
JL_max = JL_deg(:,2).';
theta(1:N_tests,1) = t1;
theta(:,2:4)       = JL_min + (JL_max - JL_min).*rand(N_tests,3);

total_Error = 0;
t0 = tic;
for i=1:N_tests

    % DH table with indices of a, d, alpha, theta
    DH = [a1 0  0    theta(i,1);
          a2 0  pi/2 theta(i,2);
          a3 0  0    theta(i,3);
          a4 0  0    theta(i,4);];

    H_0_1 = dhTransform(DH(1,1), DH(1,2), DH(1,3), DH(1,4));
    H_1_2 = dhTransform(DH(2,1), DH(2,2), DH(2,3), DH(2,4));
    H_2_3 = dhTransform(DH(3,1), DH(3,2), DH(3,3), DH(3,4));
    H_3_4 = dhTransform(DH(4,1), DH(4,2), DH(4,3), DH(4,4));

    H_0_4_test = H_0_1 * H_1_2 * H_2_3 * H_3_4;
    
    t_IK =  Hexapod_IK(linkParams, H_0_4_test);

    for j=1:2 %two solutions
        H_0_1_sol = dhTransform(DH(1,1), DH(1,2), DH(1,3), t_IK(j,1));
        H_1_2_sol = dhTransform(DH(2,1), DH(2,2), DH(2,3), t_IK(j,2));
        H_2_3_sol = dhTransform(DH(3,1), DH(3,2), DH(3,3), t_IK(j,3));
        H_3_4_sol = dhTransform(DH(4,1), DH(4,2), DH(4,3), t_IK(j,4));

        H_0_4_sol = H_0_1_sol * H_1_2_sol * H_2_3_sol * H_3_4_sol;
        
        % compare the end effector position
        pos_sol  = H_0_4_sol(1:3, 4);
        pos_test = H_0_4_test(1:3, 4);
        error    = norm(pos_test - pos_sol);
        if error > 1e-6
            %print diagnostic information
            fprintf('Test %d, Solution %d: Position error = %.6f\n', i, j, error);
            fprintf('theta [%0.3f %0.3f %0.3f %0.3f]\n', theta(:,i))
            fprintf('sol   [%0.3f %0.3f %0.3f %0.3f]\n\n', t_IK(j,:))
        end
        total_Error = total_Error + error;
    end

end

t1 = toc(t0);
fprintf("%d tests conducted over %.03f seconds \n", N_tests, t1);
% Calculate the average error over all tests
average_Error = total_Error / (N_tests * 2);
fprintf("Average error: %.4f\n", average_Error);
