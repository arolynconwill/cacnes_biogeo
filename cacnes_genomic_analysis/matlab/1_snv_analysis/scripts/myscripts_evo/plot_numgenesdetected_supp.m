function plot_numgenesdetected_supp( fdr_to_test, num_genes_detected, series_labels, plot_title, ...
    dir_save, plot_file_name )

% Colors
color_lines = 'k';
% Fonts
fs = 26;

% Figure
figure(11)
clf(11)
hold on
set(gca, 'FontSize', fs, 'FontName', 'Helvetica')
% plot
box on
for i=1:numel(series_labels)
    plot( fdr_to_test, num_genes_detected(i,:), 'LineWidth', 2, 'Color', color_lines )
    scatter( fdr_to_test, num_genes_detected(i,:), 100, color_lines, 'filled', 'MarkerFaceAlpha', 1/numel(series_labels) )
end
% x axis
xticks(fdr_to_test)
xlabel( 'false discovery rate (%)' )
% y axis
ylim([0 5])
ylabel({'number of genes','enriched for mutations'}, 'FontSize', fs) % y-axis label
% title
title( plot_title )
hold off

% Save image
if ~exist( dir_save, 'dir')
    mkdir( dir_save )
end
print([ dir_save '/' plot_file_name '.png'],'-dpng')



end