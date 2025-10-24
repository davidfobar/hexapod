% cubic_Interpolation_Eval: Evaluates cubic polynomial trajectory at specified time points.
%
% xt = cubic_Interpolation_Eval(a, t0, t, tf)
% This function evaluates a cubic polynomial trajectory at given time points
% using pre-computed coefficients. The cubic polynomial is of the form:
% x(s) = a0 + a1*s + a2*s^2 + a3*s^3 where s is normalized time in [0,1]
%
% xt = evaluated trajectory positions at time points t (nTimestep x nJoints)
%
% a = matrix of coefficients where each row corresponds to a joint and 
%     columns represent [a0, a1, a2, a3] coefficients
% t0 = initial time
% t = vector of time points to evaluate the trajectory
% tf = final time
function xt = cubic_Interpolation_Eval(a, t0, t, tf)
    xt = zeros(size(t, 2), size(a, 1)); %(nTimestep, number output)

    % for each timestep
    for i=1:numel(t)
        % determine the interpolation [0,1]
        s = (t(i)-t0)/(tf-t0);

        % compute the requried intermediate positions using a cubic function
        xt(i,:) = a(:,1) + a(:,2)*s^1 + a(:,3)*s^2 + a(:,4)*s^3;
    end
end