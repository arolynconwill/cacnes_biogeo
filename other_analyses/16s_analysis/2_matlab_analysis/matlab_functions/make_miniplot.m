function make_miniplot( these_samples_name, bar_temp, x_labels, graph_labels, legend_bool, fig_filename )

% Colors
my_colormap = flipud( ...
    [168,168,168; 
    141,211,199;
    255,255,179;
    190,186,218;
    251,128,114;
    128,177,211;
    253,180,98;
    179,222,105]/255 );

% blue-greens
% my_colormap = flipud( ...
%     [168,168,168; 
%     255,255,204;
%     199,233,180;
%     127,205,187;
%     65,182,196;
%     29,145,192;
%     34,94,168;
%     12,44,132]/256 );

% Plot

hold on
box on
set(gca,'FontSize',24)

b=bar( bar_temp', 'stacked', 'FaceColor','flat');
for k = 1:numel(graph_labels)
    b(k).CData = my_colormap(k,:);
end

if legend_bool
    legend(graph_labels, 'Location', 'northeastoutside', 'FontSize', 14)
end

xticks(1:1:size(bar_temp,2))
xticklabels(x_labels)
ax = gca;
ax.XAxis.FontSize = 8;
ax.XAxis.TickLabelInterpreter = 'none';
xtickangle(90)
%xlabel('sample')

ylim([0 1])
yticks([0 0.5 1])
ylabel('relative abundance')

title(these_samples_name)

hold off

%% Save figure

if contains( these_samples_name, 'pore strip' )
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [10 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 10 6]);
else
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 6 6]);
end

print( fig_filename, '-dpng')
    
    
end