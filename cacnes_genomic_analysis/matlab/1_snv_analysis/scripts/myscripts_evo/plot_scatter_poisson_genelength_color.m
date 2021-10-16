function plot_scatter_poisson_genelength_color( obs_genes_lengths, pval_poisson, obs_feature_by_gene, color_label, ...
    dir_save, plot_file_name )


if contains( color_label, 'dN/dS' )
    
    obs_feature_by_gene = log10( obs_feature_by_gene );

    dNdS_range_max = ceil( max( obs_feature_by_gene(obs_feature_by_gene~=Inf) ) );
    dNdS_range_min = floor( min( obs_feature_by_gene(obs_feature_by_gene~=-Inf) ) );

    feature_range = max( abs(dNdS_range_min), abs(dNdS_range_max) );

    obs_feature_by_gene(obs_feature_by_gene==Inf) = feature_range;
    obs_feature_by_gene(obs_feature_by_gene==-Inf) = -feature_range;
    
    % Colors
    my_colors = parula(100);
    colors_by_gene = zeros( numel(pval_poisson), 3 );
    for g=1:numel(pval_poisson)
        color_index = 1+round( 99*( obs_feature_by_gene(g)-(-feature_range) )/(2*feature_range) );
        colors_by_gene(g,:) = my_colors( color_index,: );
    end
    
elseif contains( color_label, 'frac' )

    feature_range = 1;
    
    % Colors
    my_colors = parula(10);
    colors_by_gene = zeros( numel(pval_poisson), 3 );
    for g=1:numel(pval_poisson)
        color_index = min( 1+floor(obs_feature_by_gene(g)*10), 10);
        colors_by_gene(g,:) = my_colors( color_index,: );
    end
    
else % num
    
    feature_range = max( obs_feature_by_gene );
    
    % Colors
    my_colors = parula(feature_range+1);
    colors_by_gene = zeros( numel(pval_poisson), 3 );
    for g=1:numel(pval_poisson)
        color_index = obs_feature_by_gene(g)+1;
        colors_by_gene(g,:) = my_colors( color_index,: );
    end
    
end


% Make figure
figure(10); 
clf(10)
hold on
box on
% scatter: not genes of interest
scatter( obs_genes_lengths, pval_poisson, 100, colors_by_gene )
set(gca,'yscale','log')
% labels
xlabel( 'gene length' )
ylabel( 'poisson pval' )
% colorbar
c=colorbar;
if contains( color_label, 'dN/dS' )
    c.Ticks = [ 0 0.5 1 ];
    c.TickLabels = [ -feature_range 1 feature_range ];
    c.Label.String = color_label;
elseif contains( color_label, 'frac' )
    c.Ticks = [ 0 0.5 1 ];
    c.Label.String = color_label;
else
    c.Ticks = 0:1/feature_range:1;
    c.TickLabels = 0:1:feature_range;
    c.Label.String = color_label;
end
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