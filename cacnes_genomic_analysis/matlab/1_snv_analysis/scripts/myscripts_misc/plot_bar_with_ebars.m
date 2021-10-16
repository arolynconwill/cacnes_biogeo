function plot_bar_with_ebars( y_values, ebar_values, ...
    plot_title, x_axis_label, x_tick_labels, dir_save, plot_file_name ) 

% Colors
color_face = .75*[ 1 1 1 ];
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
b=bar(y_values,'BaseValue',0);
b.FaceColor = color_face;
b.EdgeColor = color_edge;
b.LineWidth = 1;
box on
% x axis
xticks(1:1:numel(x_tick_labels))
xticklabels(x_tick_labels)
xlabel(x_axis_label, 'FontSize', fs)
% y axis
ylabel('molecular clock (SNVs/mb/yr)', 'FontSize', fs) % y-axis label
ylim([-0.2 0.6])
% title
title( plot_title )
% Add line for expected dN/dS if mutations are random (normalized to 1)
l=line( xlim, [1 1] );
l.Color = color_ebars;
l.LineWidth = 1;
% Add error bars
ebar_ranges = [ y_values-ebar_values, y_values+ebar_values ]
for i=1:numel(y_values)
    l_temp = line( [i i], [ebar_ranges(i,1), ebar_ranges(i,2)]);
    l_temp.Color = color_ebars;
    l_temp.LineWidth = 1.5;
end
hold off

% Save image
if ~exist( dir_save, 'dir')
    mkdir( dir_save )
end
print([ dir_save '/' plot_file_name '.png'],'-dpng')



end