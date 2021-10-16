function is_significant = plot_nummuts_hists( genes_of_interest_nummuts, genes_of_interest_nummuts_expected, ...
    plot_title, x_axis_label, x_tick_labels, ...
    pval_sig, ...
    dir_save, plot_file_name )

% Compute info for figure
num_genes_of_interest = numel(genes_of_interest_nummuts);
max_num_muts = max( max(genes_of_interest_nummuts), max(max(genes_of_interest_nummuts_expected)) );
% Histogram for expected number of mutations
nummuts_expected_hist = zeros( num_genes_of_interest, max_num_muts+1 );
for n=1:num_genes_of_interest
    nummuts_expected_hist(n,:) = histcounts( genes_of_interest_nummuts_expected(:,n), -0.5:1:max_num_muts+0.5 );
end
% Probability of observed number of mutations
prob_nummuts_expected = zeros( num_genes_of_interest, 1 ); % probability of observing this number of expected mutations
prob_nummuts_expected_text = cell( num_genes_of_interest,1 );
is_significant = zeros( num_genes_of_interest, 1 );
num_trials = size(genes_of_interest_nummuts_expected,1);
for n=1:num_genes_of_interest
    prob_of_obs = sum( genes_of_interest_nummuts_expected(:,n) >= genes_of_interest_nummuts(n) )/num_trials;
    prob_nummuts_expected(n) = prob_of_obs;
    if prob_of_obs == 0
        prob_nummuts_expected_text{n} = [ '<1/' num2str(num_trials) ];
        prob_of_obs = 1/num_trials;
    else
        prob_nummuts_expected_text{n} = num2str( prob_of_obs );
    end
    is_significant(n) = prob_of_obs < pval_sig;
end

% Figure
figure(4)
clf(4)
hold on
set(gca, 'FontSize', 16, 'FontName', 'Helvetica')
% heatmap
imagesc( nummuts_expected_hist'/num_trials )
% Colormap
my_colormap = ones( 100,3 ); my_colormap_scale = fliplr(0.1:(1-0.1)/99:1);
my_colormap = my_colormap_scale.*my_colormap'; my_colormap = my_colormap';
caxis([0 1])
colormap( my_colormap )
colorbar( 'Ticks', 0:0.25:1 )
% y axis
ylim([ 0.5 max_num_muts+0.5 ])
yticks( 1:1:max_num_muts )
ylabel('number of mutations', 'FontSize', 16) % y-axis label
% x axis
xticks( 1:1:num_genes_of_interest )
xticklabels(x_tick_labels)
xtickangle( 45 )
ax = gca;
ax.XAxis.FontSize = 12;
xlabel(x_axis_label)
% title
title( plot_title )
% Add line for expected dN/dS if mutations are random (normalized to 1)
for i=1:num_genes_of_interest
    l=line( [i-0.45,i+0.45], [genes_of_interest_nummuts(i) genes_of_interest_nummuts(i)] );
    l.Color = 'k';
    l.LineWidth = 1.5;
end
% text labels
for i=1:num_genes_of_interest
    text( i, max_num_muts, prob_nummuts_expected_text{i}, 'HorizontalAlignment', 'center', 'Color', 'k' )
end
for i=1:num_genes_of_interest
    text( i, max_num_muts-1, num2str(is_significant(i)), 'HorizontalAlignment', 'center', 'Color', 'k' )
end

% Save image
if ~exist( dir_save, 'dir')
    mkdir( dir_save )
end
if contains( plot_title,'kegg-coarse' ) || contains( plot_title,'kegg-fine' ) || contains( plot_title,'gene' )
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [10 10]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 10 10]);
    print([ dir_save '/' plot_file_name ],'-dpng')
else
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 8 6]);
    print([ dir_save '/' plot_file_name ],'-dpng')
end

end