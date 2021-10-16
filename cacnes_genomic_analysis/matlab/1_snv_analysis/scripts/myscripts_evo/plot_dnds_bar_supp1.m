function plot_dnds_bar_supp1( dNdS_for_plot, CI_for_plot, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, dir_save, plot_file_name, ...
    text_labels_1, text_labels_2 ) 

% Colors
color_face = .25*[ 1 1 1 ];
color_edge = .1*[ 1 1 1 ];
color_ebars = .1*[ 1 1 1 ];

% Figure
figure(6)
clf(6)
hold on
fs=24;
set(gca, 'FontSize', fs, 'FontName', 'Helvetica')
set(gca,'TickLabelInterpreter','none')
% bar chart
b=bar(dNdS_for_plot,'BaseValue',1);
b.FaceColor = color_face;
b.EdgeColor = color_edge;
b.FaceAlpha = 0.5;
b.LineWidth = 1;
box on
% x axis
xticks(1:1:numel(x_tick_labels))
xticklabels(x_tick_labels)
if contains( plot_title,'summary' ) || contains( plot_title,'superSLST' )
    xtickangle( 45 )
elseif contains( plot_title,'kegg' )
    xtickangle( 45 )
elseif contains( plot_title,'genes' )
    xtickangle( 45 )
end
ax = gca;
ax.YAxis.FontSize = fs-2;
ax.XAxis.FontSize = fs-2;
xlabel(x_axis_label, 'FontSize', fs )
% y axis
ylim([1/y_axis_factor y_axis_factor])
yticks([.125 .25 0.5 1 2 4 8])
ylabel('dN/dS', 'FontSize', fs) % y-axis label
set(gca,'YScale','log')
% title
title( plot_title )
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
% Add text labels
if nargin>=9
    fs_labels = 16;
    if ~iscell(text_labels_1)
        if numel(dNdS_for_plot)>1
            for i=1:numel(x_tick_labels)
                text( i, y_axis_factor/sqrt(2), num2str(text_labels_1(i)), 'HorizontalAlignment', 'center', 'Color', color_face, 'FontSize', fs_labels )
            end
        else
            text( 1, y_axis_factor/sqrt(2), num2str(text_labels_1(1)), 'HorizontalAlignment', 'center', 'Color', color_face, 'FontSize', fs_labels )
        end
    else
        for i=1:numel(x_tick_labels)
            text( i, y_axis_factor/sqrt(2), text_labels_1{i}, 'HorizontalAlignment', 'center', 'Color', color_face, 'FontSize', fs_labels )
        end
    end
end
if nargin>=10
    if ~iscell(text_labels_2)
        if numel(dNdS_for_plot)>1
            for i=1:numel(x_tick_labels)
                text( i, y_axis_factor/2, num2str(text_labels_2(i)), 'HorizontalAlignment', 'center', 'Color', color_face, 'FontSize', fs_labels )
            end
        else
            text( 1, y_axis_factor/2, num2str(text_labels_2(1)), 'HorizontalAlignment', 'center', 'Color', color_face, 'FontSize', fs_labels )
        end
    else
        for i=1:numel(x_tick_labels)
            text( i, y_axis_factor/2, text_labels_2{i}, 'HorizontalAlignment', 'center', 'Color', color_face, 'FontSize', fs_labels )
        end
    end
end
hold off

% Save image
if ~exist( dir_save, 'dir')
    mkdir( dir_save )
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
print([ dir_save '/' plot_file_name '.png'],'-dpng')


end