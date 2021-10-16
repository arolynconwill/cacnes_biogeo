function plot_heatmap( coverage_scaffold_norm_binned, cov_norm_max_for_colormap, reorder, plot_title, num_bins, bin_size )

% Plot
imagesc( coverage_scaffold_norm_binned( reorder,:), [0 cov_norm_max_for_colormap] )

% Colorbar
h=colorbar;
h.Ticks = 0:0.5:cov_norm_max_for_colormap;

% Axes and labels
temp_tick_labels = arrayfun(@(x) num2str(x), 0:0.5:cov_norm_max_for_colormap, 'UniformOutput', false);
temp_tick_labels{end} = [ temp_tick_labels{end} '+' ];
h.TickLabels = temp_tick_labels;
h.Label.String = 'coverage normalized to chromosomal coverage';
ylabel('samples')
xticks( 1:25:num_bins )
xticklabels( 0:25*bin_size:num_bins*bin_size )
xtickangle( 45 )
xlabel([ 'position on mobile element scaffold (bin size = ' num2str(bin_size) ' bp)' ])
title( plot_title, 'Interpreter', 'none' )

end