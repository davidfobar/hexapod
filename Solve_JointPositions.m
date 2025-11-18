function JointPos = Solve_JointPositions(linkParams, theta)
    % extract link parameters
    t1 = linkParams(1).theta;
    a1 = linkParams(1).a;
    a2 = linkParams(2).a;
    a3 = linkParams(3).a;
    a4 = linkParams(4).a;
    
    % DH table with indices of a, d, alpha, theta
    DH = [a1 0  0    t1;
          a2 0  pi/2 theta(1);
          a3 0  0    theta(2);
          a4 0  0    theta(3);];

    H_0_1 = dhTransform(DH(1,1), DH(1,2), DH(1,3), DH(1,4));
    H_1_2 = dhTransform(DH(2,1), DH(2,2), DH(2,3), DH(2,4));
    H_2_3 = dhTransform(DH(3,1), DH(3,2), DH(3,3), DH(3,4));
    H_3_4 = dhTransform(DH(4,1), DH(4,2), DH(4,3), DH(4,4));

    j1 = H_0_1(1:3,4);

    T2 = H_0_1*H_1_2;
    j2 = T2(1:3,4);

    T3 = T2 * H_2_3;
    j3 = T3(1:3,4);

    T4 = T3 * H_3_4;
    j4 = T4(1:3,4);
    
    JointPos = [j1 j2 j3 j4]; % Extract the joint positions from the transformation matrices
end