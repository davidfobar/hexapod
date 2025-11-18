function trajectory = buildGaitTrajectory(wayPointTimes, bodyHeight, footLiftHeight, stepLength, footSweepDist, footDistFromCoxa, linkParams)
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

    x_mid    = coxaPos(1);
    x_rear   = x_mid - stepLength/2;
    x_front  = x_mid + stepLength/2;
    x_wps    = [x_mid; x_rear; x_mid; x_front; x_mid];

    y_stance = coxaPos(2) + footDistFromCoxa*sign(coxaPos(2));
    y_swing  = y_stance + footSweepDist*sign(coxaPos(2));
    y_wps    = [y_stance; y_stance; y_swing; y_stance; y_stance];

    z_stance = -bodyHeight;
    z_swing  = z_stance + footLiftHeight;
    z_wps    = [z_stance; z_stance; z_swing; z_stance; z_stance];

    trajectory = [wayPointTimes, x_wps, y_wps, z_wps];
end