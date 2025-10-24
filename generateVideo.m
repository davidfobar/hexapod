function generateVideo(t, linkParams,  title_text, filename)
    nTimeSteps = size(t, 1);

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
   
    for i=1:nTimeSteps
                
        hold on;
        JointPos = Solve_JointPositions(linkParams, t(i,:));
        plot3(JointPos(:,1),JointPos(:,2),JointPos(:,3),'Color','k','LineWidth',2);
        hold off;

        frame = getframe(f1);
        writeVideo(v,frame);
        pause(1/30);

        ax = gca;
        cla(ax);
    end
    close(v);

end