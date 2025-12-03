function generateOrientationVideo(Q, legs, bodyHTs)

    nTimeSteps = size(Q, 2);
    nLegs      = length(legs);
    nJoints    = length(legs(1).linkParams);

    f1 = figure( ...
        'Position', [200 200 1200 800], ...
        'Color', 'w', ...
        'Renderer', 'opengl', ...
        'Visible', 'on');

    ax = axes('Parent', f1);
    hold(ax, 'on');
    grid(ax, 'on');

    view(3);
    axis(ax, [-300, 300, -300, 300, -100, 300]);
    axis(ax, 'manual');
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Hexapod Body Orientation with Fixed Feet');


    v = VideoWriter('orientation.mp4', 'MPEG-4');
    v.FrameRate = 30;
    open(v);

    colors = lines(nJoints); 
    fixedSize = [];

    % Ground plane
    groundSize = 400;
    [Xg, Yg] = meshgrid([-groundSize groundSize], [-groundSize groundSize]);
    Zg = zeros(size(Xg));

    for i = 1:nTimeSteps

        cla(ax);
        hold(ax, 'on');

        % ground
        surf(ax, Xg, Yg, Zg, ...
             'FaceColor', [0.85 0.85 0.85], ...
             'FaceAlpha', 0.4, ...
             'EdgeColor', 'none');

        % Body transform
        HT = bodyHTs(:,:,i);
        R  = HT(1:3,1:3);
        p0 = HT(1:3,4);


        L = 40;
        O = p0;
        X = p0 + R(:,1)*L;
        Y = p0 + R(:,2)*L;
        Z = p0 + R(:,3)*L;

        plot3([O(1) X(1)], [O(2) X(2)], [O(3) X(3)], 'r', 'LineWidth',2);
        plot3([O(1) Y(1)], [O(2) Y(2)], [O(3) Y(3)], 'g', 'LineWidth',2);
        plot3([O(1) Z(1)], [O(2) Z(2)], [O(3) Z(3)], 'b', 'LineWidth',2);

        % draw legs
        for j = 1:nLegs
            P = Solve_JointPositions(legs(j).linkParams, Q(:,i,j));

            % Body frame origin
            P = [[0;0;0], P];

            % World transform
            P_world = R * P + p0;

            for k = 1:nJoints
                plot3([P_world(1,k) P_world(1,k+1)], ...
                      [P_world(2,k) P_world(2,k+1)], ...
                      [P_world(3,k) P_world(3,k+1)], ...
                      '-', 'Color', colors(k,:), 'LineWidth', 2);
            end
        end


        drawnow;

        frame = getframe(f1);

        if isempty(fixedSize)
            fixedSize = size(frame.cdata);
        else
            [h,w,~] = size(frame.cdata);
            if h ~= fixedSize(1) || w ~= fixedSize(2)
                frame.cdata = imresize(frame.cdata, fixedSize(1:2));
            end
        end

        writeVideo(v, frame);
    end

    close(v);

end
