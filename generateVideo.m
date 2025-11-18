function generateVideo(Q, legs)
    nTimeSteps = size(Q, 2);
    nLegs      = length(legs);
    nJoints    = length(legs(1).linkParams);

    % generate a 3d plot that shows the motion of the leg
    f1 = figure('Position', [200 200 1200 800], 'Color','w', 'Resize','off'); 
    ax = axes('Parent',f1);
    hold(ax, 'on');
    grid(ax, 'on');
    ax.Color  = 'w';  
    ax.XColor = 'k';   
    ax.YColor = 'k';
    ax.ZColor = 'k';

    view(3);
    axis(ax, [-300, 300, -300, 300, -300, 300]);
    axis(ax, 'manual');
    xlabel('x', Color='k'); ylabel('y', Color='k'); zlabel('z', Color='k');
    title('3x Tripod Gait Cycles', Color='k');

    v = VideoWriter('./tripodGait.avi','Motion JPEG AVI');
    open(v);
    colors = lines(nJoints); 
    fixedSize = [];
    for i=1:nTimeSteps
        cla(ax);        
        hold(ax,'on');

        for j = 1:nLegs
            P = Solve_JointPositions(legs(j).linkParams, Q(:,i,j));
            % append the origin
            P = [ [0;0;0], P];
            
            for k = 1:nJoints
                plot3([P(1,k) P(1,k+1)], ...
                      [P(2,k) P(2,k+1)], ...
                      [P(3,k) P(3,k+1)], ...
                      '-', 'Color', colors(k,:), 'LineWidth', 2);
            end
        end
        %drawnow;
        frame = getframe(f1);

        if isempty(fixedSize)
            fixedSize = size(frame.cdata);
        else
            [h,w,~] = size(frame.cdata);
            if h ~= fixedSize(1) || w ~= fixedSize(2)
                frame.cdata = imresize(frame.cdata, fixedSize(1:2));
            end
        end
        writeVideo(v,frame);
    end
    close(v);
end