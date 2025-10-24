function T = Solve_FK(linkParams, theta)
    [t1, a1, a2, a3, a4] = deal(linkParams{:});
    
    % DH table with indices of a, d, alpha, theta
    DH = [a1 0  0    t1;
          a2 0  pi/2 theta(2);
          a3 0  0    theta(3);
          a4 0  0    theta(4);];

    H_0_1 = dhTransform(DH(1,1), DH(1,2), DH(1,3), DH(1,4));
    H_1_2 = dhTransform(DH(2,1), DH(2,2), DH(2,3), DH(2,4));
    H_2_3 = dhTransform(DH(3,1), DH(3,2), DH(3,3), DH(3,4));
    H_3_4 = dhTransform(DH(4,1), DH(4,2), DH(4,3), DH(4,4));

    T = H_0_1 * H_1_2 * H_2_3 * H_3_4;
end