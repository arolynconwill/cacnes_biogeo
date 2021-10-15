function make_prob_plots( prob_in_pore_mat, ...
    t_generation, t_generations_to_test, parameter_test_multiplier_tgen, ...
    drift_speed, drift_speeds_to_test, ...
    diffusion_coefficient, diffusion_coefficients_to_test, ...
    parameter_test_multiplier_vD, ...
    min_prob_to_show, ...
    my_colors, dir_figs )

%% Rail probabilities

prob_in_pore_mat_rail = prob_in_pore_mat;
prob_in_pore_mat_rail( log10(prob_in_pore_mat_rail)<2*min_prob_to_show ) = 10^(2*min_prob_to_show);


%% Fixed drift speed; vary diffusion coefficient and generation time

figure(4)
clf(4)
hold on
box on
% Plot initial condition
% Plot probabilities
for i=1:numel(diffusion_coefficients_to_test)
    next_series = squeeze(prob_in_pore_mat_rail(:,i,drift_speed==drift_speeds_to_test));
    plot(next_series,'Color',my_colors(end-i+1,:),'LineWidth',2);
end
set(gca, 'YScale', 'log')
ylim([10^(2*min_prob_to_show) 1])
% X axis
xlabel('generation time (t_{gen})')
xticks(1:1:numel(t_generations_to_test))
xticklabels( arrayfun(@(x) [ '2^{' num2str(x) '}t_{g}'], -parameter_test_multiplier_tgen:1:parameter_test_multiplier_tgen, 'UniformOutput', false ) )
x_buffer = 0;
xlim([ 1-x_buffer numel(t_generations_to_test)+x_buffer ])
% Y axis
ylabel('probability in pore after t_{gen}')
yticks( arrayfun(@(x) 10^x, min_prob_to_show:1:0 ))
yticklabels( arrayfun(@(x) [ '10^{' num2str(x) '}' ], min_prob_to_show:1:0, 'UniformOutput', false) )
ylim([(10^min_prob_to_show) 1])
% title
title('constant sebum flow')
% Appearance
set(gca,'FontSize',20)
% Legend
leg_str = arrayfun(@(x) [ '2^{' num2str(x) '} D'], -parameter_test_multiplier_vD:1:parameter_test_multiplier_vD, 'UniformOutput', false );
leg=legend(leg_str,'Location','northeastoutside');
leg.Title.String = 'diffusion coefficient (D)';
leg.FontSize = 12;

print( [ dir_figs '/' 'fig_probs_vary-D-tgen.png' ],'-dpng')


%% Fixed diffusion coefficient; vary sebum speed and generation time

figure(5)
clf(5)
hold on
box on
% Plot initial condition
% Plot probabilities
for i=1:numel(drift_speeds_to_test)
    next_series = squeeze(prob_in_pore_mat_rail(:,diffusion_coefficient==diffusion_coefficients_to_test,i));
    plot(next_series,'Color',my_colors(end-numel(drift_speeds_to_test)+i-1,:),'LineWidth',2);
end
set(gca, 'YScale', 'log')
ylim([10^(2*min_prob_to_show) 1])
% X axis
xlabel('generation time (t_{gen})')
xticks(1:1:numel(t_generations_to_test))
xticklabels( arrayfun(@(x) [ '2^{' num2str(x) '}t_{g}'], -parameter_test_multiplier_tgen:1:parameter_test_multiplier_tgen, 'UniformOutput', false ) )
x_buffer = 0;
xlim([ 1-x_buffer numel(t_generations_to_test)+x_buffer ])
% Y axis
ylabel('probability in pore after t_{gen}')
yticks( arrayfun(@(x) 10^x, min_prob_to_show:1:0 ))
yticklabels( arrayfun(@(x) [ '10^{' num2str(x) '}' ], min_prob_to_show:1:0, 'UniformOutput', false) )
ylim([(10^min_prob_to_show) 1])
% title
title('constant diffusion coefficient')
% Appearance
set(gca,'FontSize',20)
% Legend
leg_str = arrayfun(@(x) [ '2^{' num2str(x) '} v'], -parameter_test_multiplier_vD:1:parameter_test_multiplier_vD, 'UniformOutput', false );
leg=legend(leg_str,'Location','northeastoutside');
leg.Title.String = 'sebum speed (v)';
leg.FontSize = 12;

print( [ dir_figs '/' 'fig_probs_vary-v-tgen.png' ],'-dpng')


end