function pval = plot_parallelevo_hist_exp_line_obs( obs_num_genes_mutated_a_lot, sim_num_genes_mutated_a_lot, ...
    plot_title, x_axis_label, y_axis_label, dir_save, plot_file_name )

% Colors
color_obs = 0.33*[ 1 1 1 ];
color_exp = 0.67*[ 1 1 1 ];

% Quick calculations
max_num_genes = max( obs_num_genes_mutated_a_lot, max(sim_num_genes_mutated_a_lot) )+1;
pval = sum(sim_num_genes_mutated_a_lot>=obs_num_genes_mutated_a_lot)/numel(sim_num_genes_mutated_a_lot);

% Figure
figure(2)
clf(2)
hold on
box on
% bar chart for expected number of genes
histogram( sim_num_genes_mutated_a_lot, -0.5:1:max_num_genes+0.5, 'FaceColor', color_exp );
% line for observed number of genes
line( [obs_num_genes_mutated_a_lot obs_num_genes_mutated_a_lot], ylim, 'Color', color_obs, 'LineWidth', 2 )
% x axis
xlabel(x_axis_label)
xlim([-1 max_num_genes+1])
xticks( 0:1:max_num_genes)
% y axis
ylabel(y_axis_label)
% title
title( plot_title )
set(gca, 'FontSize', 20, 'FontName', 'Helvetica')
% legend
legend({'expected','observed'},'Location','northeast')
hold off
% text
text( obs_num_genes_mutated_a_lot+0.25, 0.67*max(ylim), [ 'p=' num2str(pval) ], 'FontSize', 12, 'HorizontalAlignment', 'left' )

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