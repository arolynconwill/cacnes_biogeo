function plot_dnds_bar_supp_multi( dNdS_for_plot, CI_for_plot, top_num_genes_to_test, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, dir_save, plot_file_name, ...
    text_labels_1, text_labels_2 ) 

% Colors
color_face_1 = .5*[ 1 1 1 ];
color_face_2 = .25*[ 1 1 1 ];
color_edge = .1*[ 1 1 1 ];
color_ebars = .1*[ 1 1 1 ];
% Fonts
fs = 26;

% Figure
figure(10)
clf(10)
hold on
set(gca, 'FontSize', fs, 'FontName', 'Helvetica')
% bar chart
b=bar(dNdS_for_plot,'BaseValue',1,'FaceColor','flat','EdgeColor',color_edge,'LineWidth',1);
b(1).FaceColor = color_face_1;
b(2).FaceColor = color_face_2;
box on
% x axis
if contains( plot_file_name, '_all_' )
    xlim( [ 1-.75 numel(x_tick_labels)+0.75 ] )
end
xticks(1:1:numel(x_tick_labels))
xticklabels(x_tick_labels)
if contains( x_axis_label,'KEGG' )
    xtickangle( 45 )
end
ax = gca;
if contains( x_axis_label,'KEGG' )
    ax.XAxis.FontSize = 12;
end
xlabel(x_axis_label, 'FontSize', fs)
% y axis
ylim([1/y_axis_factor y_axis_factor])
yticks([.25 0.5 1 2 4])
ylabel('dN/dS', 'FontSize', fs) % y-axis label
set(gca,'YScale','log')
% title
title( plot_title )
% legend
l=legend( arrayfun(@(x) ['top ' num2str(x)], top_num_genes_to_test, 'UniformOutput', false) );
l.Title.String = 'genes';
l.FontSize = fs/2;
l.Location = 'northwest';
% Add line for expected dN/dS if mutations are random (normalized to 1)
l_axis=line( xlim, [1 1], 'HandleVisibility', 'off' );
l_axis.Color = color_ebars;
l_axis.LineWidth = 1;
% Add error bars
for i=1:numel(x_tick_labels)
    x_buffer = 0.125;
    l_temp_1 = line( [i-x_buffer i-x_buffer], [CI_for_plot(i,1,1), CI_for_plot(i,1,2)], 'HandleVisibility', 'off');
    l_temp_1.Color = color_ebars;
    l_temp_1.LineWidth = 1.5;
    l_temp_2 = line( [i+x_buffer i+x_buffer], [CI_for_plot(i,2,1), CI_for_plot(i,2,2)], 'HandleVisibility', 'off');
    l_temp_2.Color = color_ebars;
    l_temp_2.LineWidth = 1.5;
end

hold off

% Save image
if ~exist( dir_save, 'dir')
    mkdir( dir_save )
end
print([ dir_save '/' plot_file_name '.png'],'-dpng')



end