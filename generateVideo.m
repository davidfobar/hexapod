function generateVideo(t, linkParams,  title_text, filename)
    nTimeSteps = size(t{1}, 1);
    nLegs = size(t,2);
    nJoints = size(t{1},2);

    % generate a 3d plot that shows the motion of the leg
    f1 = figure('Position', [200 200 1200 800], 'Color','w'); 
    ax = axes('Parent',f1);
    hold(ax, 'on');
    grid(ax, 'on');
    ax.Color  = 'w';  
    ax.XColor = 'k';   
    ax.YColor = 'k';
    ax.ZColor = 'k';

    view(3);
    axis([-300,300, -300,300, -300,300]);
    xlabel('x', Color='k'); ylabel('y', Color='k'); zlabel('z', Color='k');
    title(sprintf('End-Effector Cartesian Path with Orientation for %s', title_text), Color='k');

    v = VideoWriter(sprintf('./%s_arm.avi', filename),'Motion JPEG AVI');
    open(v);
    colors = lines(4); 
    for i=1:nTimeSteps
                
        hold on;

        for j = 1:nLegs
            t1_j = linkParams{1}(j);
            a1_j = linkParams{2}(j);
        
            % build leg-specific parameters
            linkParams_j = {t1_j, a1_j, linkParams{3}, linkParams{4}, linkParams{5}};

            t_thisLeg = t{j};
            P = Solve_JointPositions(linkParams_j, t_thisLeg(i,:));
            % append the origin
            P = [ [0;0;0], P];
            
            for k = 1:nJoints
                plot3([P(1,k) P(1,k+1)], ...
                      [P(2,k) P(2,k+1)], ...
                      [P(3,k) P(3,k+1)], ...
                      '-', 'Color', colors(k,:), 'LineWidth', 2);
            end
        end
        hold off;

        frame = getframe(f1);
        writeVideo(v,frame);
        pause(1/30);

        ax = gca;
        cla(ax);
    end
    close(v);

end