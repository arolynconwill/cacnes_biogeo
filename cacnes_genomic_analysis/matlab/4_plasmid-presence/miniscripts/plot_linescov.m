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
xticks( 1:10:num_bins )
xticklabels( 0:10*bin_size:num_bins*bin_size )
xtickangle( 45 )
xlabel([ 'position on mobile element scaffold (bin size = ' num2str(bin_size) ' bp)' ])
xlim([ 1 num_bins ])
ylabel('coverage normalized to chromosomal coverage')
ylim([0 2*cov_norm_max_for_colormap])
line( xlim, [coverage_cutoff coverage_cutoff], 'Color', rgb('DeepPink') )
title( plot_title, 'Interpreter', 'none' )

hold off

end