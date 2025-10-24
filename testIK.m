clear; clc;
addpath('../functions');

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
    T_test = Solve_FK(linkParams, theta);

    t_IK   = Solve_IK(linkParams, T_test);

    for j=1:2 %two solutions
        DH = [a1 0  0    t_IK(j,1);
              a2 0  pi/2 t_IK(j,2);
              a3 0  0    t_IK(j,3);
              a4 0  0    t_IK(j,4);];

        H_0_1_sol = dhTransform(DH(1,1), DH(1,2), DH(1,3), DH(1,4));
        H_1_2_sol = dhTransform(DH(2,1), DH(2,2), DH(2,3), DH(2,4));
        H_2_3_sol = dhTransform(DH(3,1), DH(3,2), DH(3,3), DH(3,4));
        H_3_4_sol = dhTransform(DH(4,1), DH(4,2), DH(4,3), DH(4,4));

        T_sol = H_0_1_sol * H_1_2_sol * H_2_3_sol * H_3_4_sol;
        
        % compare the end effector position
        pos_sol  = T_sol(1:3, 4);
        pos_test = T_test(1:3, 4);
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