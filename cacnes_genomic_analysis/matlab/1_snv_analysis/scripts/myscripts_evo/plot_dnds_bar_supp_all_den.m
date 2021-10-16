function plot_dnds_bar_supp_all_den( dNdS_for_plot, CI_for_plot, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, dir_save, plot_file_name ) 

% Colors
color_face = .5*[ 1 1 1 ];
color_edge = .1*[ 1 1 1 ];
color_ebars = .1*[ 1 1 1 ];
% Fonts
fs = 20;

% Figure
figure(10)
clf(10)
hold on
set(gca, 'FontSize', fs, 'FontName', 'Helvetica')
% bar chart
b=bar(dNdS_for_plot,'BaseValue',1,'FaceColor',color_face);
b.EdgeColor = color_edge;
b.LineWidth = 1;
box on
% x axis
xticks(1:1:numel(x_tick_labels))
xticklabels(x_tick_labels)
%ax = gca; ax.XAxis.FontSize = 12;
xlabel(x_axis_label, 'FontSize', fs)
% y axis
ylim([1/y_axis_factor y_axis_factor])
yticks([.25 0.5 1 2 4])
ylabel('dN/dS', 'FontSize', fs) % y-axis label
set(gca,'YScale','log')
% title
title( plot_title, 'FontWeight', 'normal' )
% Add line for expected dN/dS if mutations are random (normalized to 1)
l=line( xlim, [1 1] );
l.Color = color_ebars;
l.LineWidth = 1;
% Add error bars
if numel(dNdS_for_plot)>1
    for i=1:numel(x_tick_labels)
        l_temp = line( [i i], [CI_for_plot(1,i), CI_for_plot(2,i)]);
        l_temp.Color = color_ebars;
        l_temp.LineWidth = 1.5;
    end
else
    l_temp = line( [1 1], [CI_for_plot(1), CI_for_plot(2)]);
    l_temp.Color = color_ebars;
    l_temp.LineWidth = 1.5;
end

hold off

% Save image
if ~exist( dir_save, 'dir')
    mkdir( dir_save )
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 4]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 6 4]);
print([ dir_save '/' plot_file_name '.png'],'-dpng')



end