function compare_filter_to_poisson( obs_genes_lengths, pval_poisson, obs_mutated_genes_genenums, genenums_of_interest, ...
    plot_title, dir_save, plot_file_name )

not_genenums_of_interest = setdiff( obs_mutated_genes_genenums, genenums_of_interest );

% Make figure
figure(2); 
clf(2)
hold on
box on
% scatter: not genes of interest
scatter( obs_genes_lengths( ismember( obs_mutated_genes_genenums,not_genenums_of_interest ) ), pval_poisson( ismember( obs_mutated_genes_genenums,not_genenums_of_interest ) ) )
set(gca,'yscale','log')
% scatter: genes of interest
scatter( obs_genes_lengths( ismember( obs_mutated_genes_genenums,genenums_of_interest ) ), pval_poisson( ismember( obs_mutated_genes_genenums,genenums_of_interest ) ) )
set(gca,'yscale','log')
% labels
xlabel( 'gene length' )
ylabel( 'poisson pval' )
title( plot_title )
% legend
legend( { 'other genes', 'genes of interest' }, 'Location', 'southeast' ) 
% appearance
set(gca,'FontSize',16)
hold off

% Save
if nargin >4
    % Save image
    if ~exist( dir_save, 'dir')
        mkdir( dir_save )
    end
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 8 6]);
    print([ dir_save '/' plot_file_name ],'-dpng')
end

end