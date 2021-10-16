function plot_dnds_bar_supp_all_den_w_pores( dNdS_for_plot, CI_for_plot, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, dir_save, plot_file_name ) 

% Colors
color_face_intra = [.1 .4 .5];
color_face = .75*[ 1 1 1 ];
color_edge = 0*[ 1 1 1 ];
color_ebars = 0*[ 1 1 1 ];
% Fonts
fs = 20-4;

% Figure
figure(10)
clf(10)
hold on
% bar chart
b=bar(dNdS_for_plot,'BaseValue',1,'FaceColor','flat','LineWidth',1);
for i=1:numel(dNdS_for_plot)
    if i==9
        b.CData(i,:) = color_face_intra;
    else
        b.CData(i,:) = color_face;
    end
%     b(i).EdgeColor = color_edge;
%     b(i).LineWidth = 1;
end
box on
% Axes
set(gca, 'FontSize', fs, 'FontName', 'Helvetica', 'LineWidth', 1)
% x axis
%xticks(1:1:numel(x_tick_labels))
xticks([1 3 4 5 6 7 9])
xticklabels(x_tick_labels([1 3 4 5 6 7 9]))
xlabel(x_axis_label, 'FontSize', fs)
% y axis
ylim([1/y_axis_factor y_axis_factor])
ylabel('dN/dS', 'FontSize', fs) % y-axis label
set(gca,'YScale','log')
y_tick_nums = 0.5:0.1:2;
yticks( y_tick_nums )
y_tick_labels = arrayfun(@(x) {''}, y_tick_nums );
y_tick_labels{1} = '0.5'; y_tick_labels{end} = '2'; y_tick_labels{6} = '1';
yticklabels(y_tick_labels)
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
set(gcf, 'PaperSize', [6.5 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 6.5 5]);
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [9 5.75]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 9 5.75]);
print([ dir_save '/' plot_file_name '.png'],'-dpng')



end