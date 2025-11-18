% constAccelInterp - Constant Acceleration Interpolation
%
% [p,v,a] = constAccelInterp(t, trajectory, transPercent)
%  This function performs constant acceleration interpolation along a trajectory and
%  returns the position, velocity, and acceleration at time t.
%
% t = the time at which to evaluate the trajectory
% trajectory = an Nx(M+1) matrix where the first column contains time waypoints and the
%              remaining M columns contain position waypoints in M-dimensional space
% transPercent = the percentage of each segment's duration to allocate for acceleration/deceleration
%
% p = the interpolated position at time t (1xM vector)
% v = the interpolated velocity at time t (1xM vector)
% a = the interpolated acceleration at time t (1xM vector)
%
% David Fobar
% CWID - 10950737
% ROBO554
% 04NOV2025
function [p,v,a] = constAccelInterp(t, trajectory, transPercent)
    T = trajectory(:,1);     % times
    P = trajectory(:,2:end); % paths
    [N,M] = size(P);         % numWaypoints, numDimensions

    % default outputs
    p = zeros(1, M); % position vector
    v = zeros(1, M); % velocity vector
    a = zeros(1, M); % acceleration vector

    % confirm t is in range, return extreme position if not, assume 
    % velocity and acceleration are zero
    if t <= T(1)
        p = P(1,:);
        return;
    elseif t >= T(end)
        p = P(end,:);
        return;
    end

    % compute constant velocities
    dt = diff(T);
    vc = diff(P)./dt;

    % compute transition time periods (tau)
    tau = zeros(N-2,1);
    for i = 1:(N-2)
        tau(i) = transPercent*min(dt(i), dt(i+1));
    end

    % starting from the initial time/position, determine where in the path
    % we are at for time t
    p_cur = P(1,:);

    % iterate through each waypoint
    for i = 1:(N-1)
        % get transition time lengths, which are equal to zero for the
        % first and last waypoints
        tau_prev = 0;
        if (i>1) 
            tau_prev = tau(i-1); 
        end
        tau_next = 0;
        if (i<=N-2)
            tau_next = tau(i);
        end

        % determine times for this segment
        t_lin_start  = T(i)   + tau_prev;
        t_tran_start = T(i+1) - tau_next;
        t_lin_end    = t_tran_start;
        t_tran_end   = T(i+1) + tau_next;

        % velocity for this segment
        v_i = vc(i,:);

        % assume that the segment motion starts with linear motion, if t
        % falls after the waypoint, but before linear motion begins, it
        % will get caught in the constant acceleration portion
        if t <= t_lin_end
            dt_seg = t - t_lin_start;
            p = p_cur + v_i*dt_seg;
            v = v_i;
            a = zeros(1,M);
            return;
        else
            % we must be beyond the linear segment, advance our position
            % accordingly
            p_cur = p_cur + v_i*max(0,t_lin_end-t_lin_start);
        end
    
        % check to see if t falls within this segments constant 
        % acceleration portion
        a_i = (vc(i+1,:) - v_i)/(2*tau_next);
        if t <= t_tran_end
            dt_seg = t - t_tran_start;
            p = p_cur + v_i*dt_seg + 0.5*a_i*(dt_seg.^2);
            v = v_i + a_i*dt_seg;
            a = a_i;
            return;
        else
            % update the position for the end of the constant accel segment
            dt_seg = 2*tau_next;
            p_cur = p_cur + v_i*dt_seg + 0.5*a_i*(dt_seg.^2);
        end
    end
end