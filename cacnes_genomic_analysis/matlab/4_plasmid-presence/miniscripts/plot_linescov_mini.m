function plot_linescov( coverage_scaffold_norm_binned, cov_norm_max_for_colormap, coverage_cutoff, plot_title, num_samples, num_bins, bin_size )

hold on

box on

% Colors to cycle through
my_colors = lines(7);
my_color_alpha = 0.2;

% Line for each sample
for n=1:num_samples
    plot( coverage_scaffold_norm_binned(n,:), 'Color', [my_colors(mod(n,7)+1,:), my_color_alpha] );
end

% Labels
xticks( [1 num_bins] )
xticklabels( [1 num_bins*bin_size] )
xlabel([ 'position' ])
xlim([ 1 num_bins ])
ylabel('coverage (norm)')
ylim([0 3*cov_norm_max_for_colormap])
line( xlim, [coverage_cutoff coverage_cutoff], 'Color', rgb('DeepPink') )
title( plot_title, 'Interpreter', 'none' )

hold off

end