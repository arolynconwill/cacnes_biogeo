function plot_parallelevo_hist_exp_line_obs_poisson( pval_cutoff_from_bh, obs_pvals_of_mutated_genes, genes_of_interest_bool, ...
    plot_title, x_axis_label, y_axis_label, line_label, dir_save, plot_file_name )

% Colors
colors_list = lines(2);
color_obs_others = colors_list(1,:);
color_obs_genes_of_interest = colors_list(2,:);
color_cutoff = 0.33*[ 1 1 1 ];

% Quick calculations
max_pval = ceil( max( pval_cutoff_from_bh, max(obs_pvals_of_mutated_genes) ) );
min_pval = floor( min( pval_cutoff_from_bh, min(obs_pvals_of_mutated_genes) ) );
num_sig_genes = sum( genes_of_interest_bool );

% Figure
figure(3)
clf(3)
hold on
box on
% bar chart for expected number of genes
histogram( obs_pvals_of_mutated_genes( ~genes_of_interest_bool ), min_pval:0.25:max_pval, 'FaceColor', color_obs_others );
histogram( obs_pvals_of_mutated_genes( genes_of_interest_bool ), min_pval:0.25:max_pval, 'FaceColor', color_obs_genes_of_interest );
% line for observed number of genes
line( [pval_cutoff_from_bh pval_cutoff_from_bh], ylim, 'Color', color_cutoff, 'LineWidth', 2 )
% x axis
xlabel(x_axis_label)
xlim([min_pval-1 max_pval+1])
xticks( min_pval:1:max_pval )
% y axis
ylabel(y_axis_label)
% title
title( plot_title )
set(gca, 'FontSize', 20, 'FontName', 'Helvetica')
% legend
l=legend({'other genes','genes of interest','pval cutoff (BH)'},'Location','northwest');
set(l, 'Interpreter', 'none')
hold off
% text
text( pval_cutoff_from_bh-0.25, 0.67*max(ylim), [ num2str(num_sig_genes) ' genes' ], 'FontSize', 14, 'HorizontalAlignment', 'right' )

% Save image
if ~exist( dir_save, 'dir')
    mkdir( dir_save )
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
print([ dir_save '/' plot_file_name ],'-dpng')

end