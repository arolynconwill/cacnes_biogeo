function plot_parallelevo_doublebar_supp( obs_list_cum, exp_list_cum, exp_list_cum_stddev, plot_title, x_axis_label, dir_save, plot_file_name, make_inset, inset_start ) 

% Colors
color_obs = 0.33*[ 1 1 1 ];
color_exp = 0.67*[ 1 1 1 ];

fs = 26;

% Figure
figure(10)
clf(10)
hold on
box on
% bar chart
b = bar( [obs_list_cum;exp_list_cum]', 'FaceColor', 'flat' );
b(1).CData = color_obs;
b(2).CData = color_exp;
% ebars
ebar_shift = 0.15;
for i=1:numel(exp_list_cum)
    if i==1
        line( [i+ebar_shift,i+ebar_shift], [exp_list_cum(i)-exp_list_cum_stddev(i) exp_list_cum(i)+exp_list_cum_stddev(i)], 'Color', 'k', 'LineWidth', 1.5 )
    else
        line( [i+ebar_shift,i+ebar_shift], [exp_list_cum(i)-exp_list_cum_stddev(i) exp_list_cum(i)+exp_list_cum_stddev(i)], 'Color', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off' )
    end
end
% x axis
xlabel(x_axis_label);
if contains( plot_file_name, 'short' )
    fprintf(1,'hello')
    xlim( [ 1.5 numel(obs_list_cum)+0.5 ] )
    xticks( 2:1:numel(obs_list_cum))    
else
    xticks( 1:1:numel(obs_list_cum))
end
% y axis
ylabel('number of genes', 'FontSize', fs) % y-axis label
if contains( plot_file_name, 'short' )
    ylim( [0 2*max([obs_list_cum(2:end),exp_list_cum(2:end)]) ] )
else
    ylim( [0 1.2*max([obs_list_cum,exp_list_cum]) ] )
end
% title
title( plot_title, 'Interpreter', 'none' )
set(gca, 'FontSize', fs, 'FontName', 'Helvetica')
% legend
if contains( plot_file_name, 'short' )
    legend({'observed','expected'},'Location','northwest')
else
    legend({'observed','expected'},'Location','northeast')
end
% inset
if make_inset
    axes('Position',[.5 .25 .375 .35])
    box on
    b2=bar( [obs_list_cum(inset_start:end);exp_list_cum(inset_start:end)]', 'FaceColor', 'flat' );
    b2(1).CData = color_obs;
    b2(2).CData = color_exp;
    xticklabels(inset_start:1:numel(obs_list_cum))
    for i=inset_start:numel(obs_list_cum)
        line( [i+0.125-inset_start+1,i+0.125-inset_start+1], [exp_list_cum(i)-exp_list_cum_stddev(i) exp_list_cum(i)+exp_list_cum_stddev(i)], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off' )
    end
    ylim( [0 1+max([obs_list_cum(inset_start:end),exp_list_cum(inset_start:end)]) ] )
    yticks( 0:1:1+max([obs_list_cum(inset_start:end),exp_list_cum(inset_start:end)]) )
    set(gca,'FontSize',16)
end
hold off

% Save image
if ~exist( dir_save, 'dir')
    mkdir( dir_save )
end
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [8 6]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 8 6]);
print([ dir_save '/' plot_file_name '.png'],'-dpng')



end