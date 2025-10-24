function generatePlot(t, tEvalList, title_text, filename)
    f1 = figure('Position', [200 200 1200 800], 'Color','w');
    ax = axes('Parent',f1);
    hold(ax, 'on');
    grid(ax, 'on');
    ax.Color  = 'w';  
    ax.XColor = 'k';   
    ax.YColor = 'k';

    nJoints           = size(t,2);
    ax.GridColor      = [0.7 0.7 0.7];
    ax.MinorGridColor = [0.85 0.85 0.85];
    colors            = lines(nJoints);

    for i = 1:nJoints 
        plot(ax, tEvalList, t(:,i), 'LineWidth', 1.4, 'Color', colors(i,:));
    end 
    
    xlabel('Time [s]', Color='k'); 
    ylabel('Joint variable (θ)', Color='k'); 
    title(title_text, Color='k'); 
    legend(arrayfun(@(j)sprintf('Joint %d',j),1:nJoints,'UniformOutput',false),'Location', 'southeast','Color','w','TextColor','k');

    exportgraphics(f1, sprintf('./%s_joint.png', filename),...
    'Resolution', 300,'BackgroundColor','white');
end