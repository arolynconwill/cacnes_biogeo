function [ h_benjamini_hochberg_reindexed, pval_max ] = benjamini_hochberg_sig_filter( pvalues, num_tests, fdr_bh, dir_save, plot_file_name )

% Benjamini-Hochberg procedure
num_pvals = numel(pvalues); 
% Sort pvalues
[ pvalues_sorted, pvalues_index ] = sort( pvalues );
% Compute BH threshold for each index
bh_cutoffs = fdr_bh*(1:1:num_pvals)'/num_tests;
% Determine if pvalue is \leq its BH threshold
bh_bool = ( pvalues_sorted <= bh_cutoffs );
if sum( bh_bool )>0
    pval_max = max( pvalues_sorted( bh_bool ) );
else
    pval_max = 0;
end
% Boolean vector for which pvals are significant
h_benjamini_hochberg = zeros( size( pvalues ), 'logical' ); % initialize
if sum( bh_bool )>0
    h_benjamini_hochberg( pvalues_sorted <= pval_max ) = 1; % h=0 accept null; h=1 accept alternative
end
% Reindex
h_benjamini_hochberg_reindexed = zeros( size( pvalues ), 'logical' ); % initialize
h_benjamini_hochberg_reindexed( pvalues_index( h_benjamini_hochberg ) ) = 1;

% Plot prep
index_list = 1:1:num_pvals;
inset_index = min( ceil( 1.5*sum(h_benjamini_hochberg) ), num_pvals );

% Plot
figure(1)
clf(1)
hold on
% visuals
dot_size = 50;
lw = 1.5;
fs = 24;
colors_list = lines(2);
color_dot = colors_list(1,:);
color_dot_sig = colors_list(2,:);
color_line = 0.2*[ 1 1 1 ];
% axes
box on
title(['benjamini-hochberg (fdr=' num2str(fdr_bh) ')'])
xlabel('index')
xlim( [ 0 numel(pvalues) ] )
ylabel('pvalue')
ylim( [ 0 1 ] )
set(gca,'FontSize',fs)
% pvalues
scatter( index_list(~h_benjamini_hochberg), pvalues_sorted(~h_benjamini_hochberg), dot_size, color_dot )
scatter( index_list(h_benjamini_hochberg), pvalues_sorted(h_benjamini_hochberg), dot_size, color_dot_sig )
% benjamini-hochberg threshold
x_bh = 0:1:numel(pvalues);
y_bh = fdr_bh*x_bh/num_tests;
plot( x_bh, y_bh, 'Color', color_line, 'LineWidth', lw )
% inset
if inset_index>0
    axes('Position',[.2 .5 .35 .35])
    hold on
    box on
    plot( x_bh(1:1:inset_index+1), y_bh(1:1:inset_index+1), 'Color', color_line, 'LineWidth', lw )
    scatter( index_list(~h_benjamini_hochberg(1:1:inset_index)), pvalues_sorted(~h_benjamini_hochberg(1:1:inset_index)), dot_size, color_dot )
    scatter( index_list(h_benjamini_hochberg(1:1:inset_index)), pvalues_sorted(h_benjamini_hochberg(1:1:inset_index)), dot_size, color_dot_sig )
    xlim([ 0 inset_index ])
    title([ num2str(sum(h_benjamini_hochberg)) '/' num2str(num_tests) ])
    set(gca,'FontSize',fs-2)
    hold off
else
    inset_index = min( 10, num_pvals );
    axes('Position',[.2 .5 .35 .35])
    hold on
    box on
    plot( x_bh(1:1:inset_index+1), y_bh(1:1:inset_index+1), 'Color', color_line, 'LineWidth', lw )
    scatter( index_list(~h_benjamini_hochberg(1:1:inset_index)), pvalues_sorted(~h_benjamini_hochberg(1:1:inset_index)), dot_size, color_dot )
    scatter( index_list(h_benjamini_hochberg(1:1:inset_index)), pvalues_sorted(h_benjamini_hochberg(1:1:inset_index)), dot_size, color_dot_sig )
    xlim([ 0 inset_index ])
    title([ num2str(sum(h_benjamini_hochberg)) '/' num2str(num_tests) ])
    set(gca,'FontSize',fs-2)
    hold off    
end
hold off

% Save
if nargin >3
    % Save image
    if ~exist( dir_save, 'dir')
        mkdir( dir_save )
    end
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 8 6]);
    print([ dir_save '/' plot_file_name '.png' ],'-dpng')
end


end
