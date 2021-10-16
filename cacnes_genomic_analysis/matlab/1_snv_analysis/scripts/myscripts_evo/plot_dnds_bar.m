function plot_dnds_bar( dNdS_for_plot, CI_for_plot, ...
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
set(gca, 'FontSize', 16, 'FontName', 'Helvetica')
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
ax.XAxis.FontSize = 12;
xlabel(x_axis_label)
% y axis
ylim([1/y_axis_factor y_axis_factor])
yticks([.25 0.5 1 2 4])
ylabel('dN/dS', 'FontSize', 16) % y-axis label
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
    if ~iscell(text_labels_1)
        if numel(dNdS_for_plot)>1
            for i=1:numel(x_tick_labels)
                text( i, 0.9*y_axis_factor, num2str(text_labels_1(i)), 'HorizontalAlignment', 'center', 'Color', color_face )
            end
        else
            text( 1, 0.9*y_axis_factor, num2str(text_labels_1(1)), 'HorizontalAlignment', 'center', 'Color', color_face )
        end
    else
        for i=1:numel(x_tick_labels)
            text( i, 0.9*y_axis_factor, text_labels_1{i}, 'HorizontalAlignment', 'center', 'Color', color_face )
        end
    end
end
if nargin>=10
    if ~iscell(text_labels_2)
        if numel(dNdS_for_plot)>1
            for i=1:numel(x_tick_labels)
                text( i, 0.8*y_axis_factor, num2str(text_labels_2(i)), 'HorizontalAlignment', 'center', 'Color', color_face )
            end
        else
            text( 1, 0.8*y_axis_factor, num2str(text_labels_2(1)), 'HorizontalAlignment', 'center', 'Color', color_face )
        end
    else
        for i=1:numel(x_tick_labels)
            text( i, 0.8*y_axis_factor, text_labels_2{i}, 'HorizontalAlignment', 'center', 'Color', color_face )
        end
    end
end
hold off

% Save image
if ~exist( dir_save, 'dir')
    mkdir( dir_save )
end
if contains( plot_title,'kegg-coarse' ) || contains( plot_title,'kegg-fine' ) || contains( plot_title,'gene' )
    alt_width = max(6,ceil(numel(dNdS_for_plot))/2);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [max(10,alt_width) 10]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 max(10,alt_width) 10]);
    print([ dir_save '/' plot_file_name '.png'],'-dpng')
else
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 8 6]);
    print([ dir_save '/' plot_file_name '.png'],'-dpng')
end


end