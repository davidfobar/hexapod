% compute_Cubic_Interpolation_Coef_for_Zero_Velocity: Computes cubic polynomial coefficients for 
% trajectory interpolation with zero boundary velocities.
%
% a = compute_Cubic_Interpolation_Coef_for_Zero_Velocity(xi, xf)
% This function calculates the coefficients for a cubic polynomial trajectory
% that interpolates between initial and final positions with zero velocity
% constraints at both endpoints. The cubic polynomial is of the form:
% x(t) = a0 + a1*t + a2*t^2 + a3*t^3 (normalized time t in [0,1])
%
% a = matrix of coefficients where each row corresponds to a joint and 
%     columns represent [a0, a1, a2, a3] coefficients
%
% xi = initial position(s) for each joint
% xf = final position(s) for each joint
function a = compute_Cubic_Interpolation_Coef_for_Zero_Velocity(xi, xf)
    % compute the coefficients requried for the displacement, assuming
    % zero velocity at start and end points
    a0 = xi;
    a1 = zeros(size(xi));
    a2 = 3*(xf - xi);
    a3 = -2*(xf - xi);
    a  = [a0; a1; a2; a3].'; % (joint #, a(#))
end