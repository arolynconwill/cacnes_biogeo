function make_prob_plots( prob_in_pore_mat, ...
    t_vals, t_generation, tgen_index, ...
    drift_speed, drift_speeds_to_test, ...
    diffusion_coefficient, diffusion_coefficients_to_test, ...
    parameter_test_multiplier_vD, ...
    min_prob_to_show, ...
    my_colors, dir_figs )

%% Rail probabilities

prob_in_pore_mat_rail = prob_in_pore_mat;
for d=1:numel(diffusion_coefficients_to_test)
    for v=1:numel(drift_speeds_to_test)
        this_prob_in_pore = prob_in_pore_mat{ d,v };
        this_prob_in_pore( log10(this_prob_in_pore)<2*min_prob_to_show ) = 10^(2*min_prob_to_show);
        prob_in_pore_mat_rail{ d,v} = this_prob_in_pore;
    end
end


%% Fixed drift speed; vary diffusion coefficient and generation time

figure(4)
clf(4)
hold on
box on
% Plot initial condition
% Plot probabilities
for i=1:numel(diffusion_coefficients_to_test)
    next_series = prob_in_pore_mat_rail{ i, drift_speed==drift_speeds_to_test };
    plot( t_vals, next_series, 'Color',my_colors(end-i+1,:),'LineWidth',2 );
end
set(gca, 'YScale', 'log')
% X axis
xlabel('time (days)')
% Y axis
ylabel('probability in pore')
yticks( arrayfun(@(x) 10^x, min_prob_to_show:1:0 ))
yticklabels( arrayfun(@(x) [ '10^{' num2str(x) '}' ], min_prob_to_show:1:0, 'UniformOutput', false) )
ylim([(10^min_prob_to_show) 1])
% title
title('constant sebum flow')
% Appearance
set(gca,'FontSize',20)
% Legend
leg_str = arrayfun(@(x) [ num2str(round(diffusion_coefficients_to_test(x+parameter_test_multiplier_vD+1))) ' um^2/day (2^{' num2str(x) '}D)'], -parameter_test_multiplier_vD:1:parameter_test_multiplier_vD, 'UniformOutput', false );
leg=legend(leg_str,'Location','northeastoutside');
leg.Title.String = 'diffusion coefficient';
leg.FontSize = 12;
% Line for estimated generation time
ylim_temp = ylim;
line( [t_generation t_generation], ylim_temp, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'HandleVisibility', 'off' )

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
print( [ dir_figs '/' 'fig_probs_vary-D-tgen.png' ],'-dpng')


%% Fixed diffusion coefficient; vary sebum speed and generation time

figure(5)
clf(5)
hold on
box on
% Plot initial condition
% Plot probabilities
for i=1:numel(drift_speeds_to_test)
    next_series = prob_in_pore_mat_rail{ diffusion_coefficient==diffusion_coefficients_to_test, i };
    next_series(end)
    plot( t_vals, next_series, 'Color', my_colors(i,:),'LineWidth',2 );
end
set(gca, 'YScale', 'log')
% X axis
xlabel('time (days)')
% Y axis
ylabel('probability in pore')
yticks( arrayfun(@(x) 10^x, min_prob_to_show:1:0 ))
yticklabels( arrayfun(@(x) [ '10^{' num2str(x) '}' ], min_prob_to_show:1:0, 'UniformOutput', false) )
ylim([(10^min_prob_to_show) 1])
% title
title('constant diffusion coefficient')
% Appearance
set(gca,'FontSize',20)
% Legend
leg_str = arrayfun(@(x) [ num2str(round(drift_speeds_to_test(x+parameter_test_multiplier_vD+1))) ' um/day (2^{' num2str(x) '} v)'], -parameter_test_multiplier_vD:1:parameter_test_multiplier_vD, 'UniformOutput', false );
leg=legend(leg_str,'Location','northeastoutside');
leg.Title.String = 'sebum speed';
leg.FontSize = 12;
% Line for estimated generation time
ylim_temp = ylim;
line( [t_generation t_generation], ylim_temp, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'HandleVisibility', 'off' )

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
print( [ dir_figs '/' 'fig_probs_vary-v-tgen.png' ],'-dpng')


end