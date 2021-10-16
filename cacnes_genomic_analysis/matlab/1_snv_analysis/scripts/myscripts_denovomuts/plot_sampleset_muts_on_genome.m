function plot_sampleset_muts_on_genome( genome_length, mutation_positions, this_cluster_name, directory_name, fig_num, save_plots )

% Make a plot
fig = figure(fig_num);
clf(fig_num);
%
% Line plot
subplot(2,1,1)
hold on
box on
xlim([0,genome_length])
xlabel('Genome Position (bp)') % x-axis label
ylim([0,1])
yticks([])
for i=1:numel(mutation_positions)
    mark_pos = mutation_positions(i);
    next_line = line([mark_pos,mark_pos],ylim,'Color',[0 0 0 .25],'LineStyle','-','LineWidth',1);
end
title([ 'Mutation Positions: ' this_cluster_name ]);
set(gca, 'FontSize', 12, 'FontName', 'Verdana')
hold off
%
% Histogram
subplot(2,1,2)
hold on
box on
xlim([0,genome_length])
xlabel('Genome Position (bp)') % x-axis label
histogram(mutation_positions,0:25000:genome_length,'FaceColor',0.75*[1 1 1])
title([ 'Mutation Positions: ' this_cluster_name ]);
set(gca, 'FontSize', 12, 'FontName', 'Verdana')
hold off

% Save figure
if save_plots
    % Provide dimensions
    printHeight = 4; % inches?
    printWidth = 8;
    % Set dimensions
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [printWidth printHeight]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 printWidth printHeight]);
    % Save
    print([pwd '/' directory_name '/' 'Figure_MutationPositions_' this_cluster_name '.png'],'-dpng')
end

end