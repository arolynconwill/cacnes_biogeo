function make_heatmap_fig( prob_in_pore_mat, min_prob_to_show, ...
    parameter_test_multiplier_vD, drift_speeds_to_test, diffusion_coefficients_to_test, ...
    my_colors, dir_save )


figure(1); clf(1)
% Probabilities in log10
prob_in_pore_mat_log10 = log10(prob_in_pore_mat);
prob_in_pore_mat_log10( prob_in_pore_mat_log10<min_prob_to_show ) = min_prob_to_show;
% Make heatmap
imagesc(flipud(prob_in_pore_mat_log10))
colormap( flipud(my_colors) )
c=colorbar;
c.Limits = [ min_prob_to_show 0 ];
c.Ticks = min_prob_to_show:1:0;
c_tick_labels = arrayfun(@(x) num2str(x), min_prob_to_show:1:0, 'UniformOutput', false );
c_tick_labels{1} = [ '<' num2str(min_prob_to_show) ];
c.TickLabels = c_tick_labels;
c.Label.String = 'probability in pore (log10)';
% X axis
xticks(1:1:numel(drift_speeds_to_test))
%xticklabels(drift_speeds_to_test)
xticklabels( arrayfun(@(x) [ '2^{' num2str(x) '}v'], -parameter_test_multiplier_vD:1:parameter_test_multiplier_vD, 'UniformOutput', false ) )
xlabel('sebum speed (um/day)')
xtickangle(90)
% Y axis
yticks(1:1:numel(diffusion_coefficients_to_test))
%yticklabels(diffusion_coefficients_to_test)
yticklabels( fliplr( arrayfun(@(x) [ '2^{' num2str(x) '}D'], -parameter_test_multiplier_vD:1:parameter_test_multiplier_vD, 'UniformOutput', false ) ) )
ylabel('diffusion coefficient (um^2/day)')
% Title
title('probability in pore after one generation time')
% Appearance
set(gca,'FontSize',20)

% Save
print( [ dir_save '/' 'fig_prob-grid.png'],'-dpng')


end