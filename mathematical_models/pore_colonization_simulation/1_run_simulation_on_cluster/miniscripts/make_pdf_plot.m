function make_pdf_plots( diff_eqn_solutions, min_prob_to_show, ...
    drift_speed, drift_speeds_to_test, diffusion_coefficient, diffusion_coefficients_to_test, ...
    parameter_test_multiplier_vD, ...
    my_colors, dir_save )

% Fixed drift speed; vary diffusion coefficient

figure(2)
clf(2)
hold on
box on
% Plot initial condition
% Plot pdfs 
for i=1:numel(diffusion_coefficients_to_test)
    next_series = diff_eqn_solutions{ i, drift_speed==drift_speeds_to_test };
%     if i==1 % plot initial condition
%         next_vec = next_series(1,:); 
%         plot(next_vec,'Color','k','LineWidth',2);
%     end
    next_vec =  next_series(end,:); 
    plot(next_vec,'Color',my_colors(end-i+1,:),'LineWidth',2);
end
set(gca, 'YScale', 'log')
%ylim([10^(2*min_prob_to_show) 1])
% X axis
xlabel('depth into pore (um) after one generation time')
xlim([0 1000])
% Y axis
ylabel('probability')
%ylim([0 5*10^-3])
ylim([(10^min_prob_to_show)^2 1])
% title
title('constant sebum flow')
% Appearance
set(gca,'FontSize',20)
% Legend
leg_str = arrayfun(@(x) [ '2^{' num2str(x) '}D'], -parameter_test_multiplier_vD:1:parameter_test_multiplier_vD, 'UniformOutput', false );
leg=legend(leg_str);
leg.Title.String = 'diffusion coefficient (D)';

print( [ dir_save '/' 'fig_series-vary-D_log.png' ],'-dpng')

% Fixed drift speed; vary diffusion coefficient

figure(3)
clf(3)
hold on
box on
% Plot initial condition
% Plot pdfs 
for i=1:numel(drift_speeds_to_test)
    next_series = diff_eqn_solutions{ diffusion_coefficient==diffusion_coefficients_to_test, i };
%     if i==1 % plot initial condition
%         next_vec = next_series(1,:); 
%         plot(next_vec,'Color','k','LineWidth',2);
%     end
    next_vec =  next_series(end,:); 
    plot(next_vec,'Color',my_colors(end+i-numel(drift_speeds_to_test),:),'LineWidth',2);
end
set(gca, 'YScale', 'log')
%ylim([10^(2*min_prob_to_show) 1])
% X axis
xlabel('depth into pore (um) after one generation time')
xlim([0 1000])
% Y axis
ylabel('probability')
%ylim([0 5*10^-3])
ylim([(10^min_prob_to_show)^2 1])
% title
title('constant diffusion coefficient')
% Appearance
set(gca,'FontSize',20)
% Legend
leg_str = arrayfun(@(x) [ '2^{' num2str(x) '}v'], -parameter_test_multiplier_vD:1:parameter_test_multiplier_vD, 'UniformOutput', false );
leg=legend(leg_str);
leg.Title.String = 'sebum speed (v)';

print( [ dir_save '/' 'fig_series-vary-v_log.png' ],'-dpng')

end