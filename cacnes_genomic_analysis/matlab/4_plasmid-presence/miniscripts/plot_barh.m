function plot_barh( bar_data, reorder, x_limits, x_label, plot_title )

% Colors
color_for_bars = 0.66*[ 1 1 1 ];

% Bar chart
barh( bar_data(fliplr(reorder)), 'FaceColor', color_for_bars )
xlim([0 1])

% Labels
xlabel(x_label)
xlim(x_limits)
if isequal(plot_title,'coverage (norm)')
    xticks([ x_limits(1) 1 x_limits(2) ])
    xticklabels({ num2str(x_limits(1)), '1', [num2str(x_limits(2)) '+'] } )
end
yticks([])
title(plot_title)

end
