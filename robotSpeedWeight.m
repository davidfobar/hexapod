% Returns the speed that the robot can walk at in m/s and the total weight
% it can lift in kg

function [speed, weight] = robotSpeedWeight(torques, totalTime, stepLength, mass, servoTorque)
    % stepLength in mm
    % totalTime in s
    % torques in kg*cm

    maxTorque_J1 = zeros(1,6);
    maxTorque_J2 = zeros(1,6);
    maxTorque_J3 = zeros(1,6);
    biggestLegVal = 0;

    speed = stepLength/1000 / totalTime; %m/s
    
    % Find the maximum torques for each joint in each leg
    for i = 1:size(torques, 3)
        torques_leg = torques(:,:, i);
        for j = 1:size(torques, 1)
            torques_joint = torques_leg(j, :);
            if j == 1 && max(torques_joint) > maxTorque_J1(i)
                maxTorque_J1(i) = max(torques_joint);
            elseif j == 2 && max(torques_joint) > maxTorque_J2(i)
                maxTorque_J2(i) = max(torques_joint);
            elseif j == 3 && max(torques_joint) > maxTorque_J3(i)
                maxTorque_J3(i) = max(torques_joint);
            end
        end
        
        if biggestLegVal < (maxTorque_J1(i) + maxTorque_J2(i) + maxTorque_J3(i))
            biggestLegVal = (maxTorque_J1(i) + maxTorque_J2(i) + maxTorque_J3(i));
        end
    end

    extraTorque = 3*servoTorque - biggestLegVal;
    weightPerTorque = mass/biggestLegVal;
    weight = weightPerTorque * extraTorque;