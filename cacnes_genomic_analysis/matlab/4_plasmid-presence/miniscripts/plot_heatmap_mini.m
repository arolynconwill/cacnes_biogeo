function plot_heatmap_mini( coverage_scaffold_norm_binned, cov_norm_max_for_colormap, reorder, plot_title, num_bins, bin_size, this_set_names )

% Plot
imagesc( coverage_scaffold_norm_binned( reorder,:), [0 cov_norm_max_for_colormap] )

% Colorbar
h=colorbar;
h.Ticks = 0:0.5:cov_norm_max_for_colormap;

% Axes and labels
temp_tick_labels = arrayfun(@(x) num2str(x), 0:0.5:cov_norm_max_for_colormap, 'UniformOutput', false);
temp_tick_labels{end} = [ temp_tick_labels{end} '+' ];
h.TickLabels = temp_tick_labels;
xticks( [1 num_bins] )
xticklabels( [1 num_bins*bin_size] )
yticks( 1:1:numel(this_set_names) )
yticklabels( this_set_names(reorder) )
title( plot_title, 'Interpreter', 'none' )
set(gca,'TickLabelInterpreter','none')

end