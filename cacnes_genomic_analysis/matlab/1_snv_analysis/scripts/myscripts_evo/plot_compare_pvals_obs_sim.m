function plot_compare_pvals_obs_sim( num_genes_mutated, sim_gene_genenums, obs_mutated_genes_genenums, ...
    obs_num_muts_per_gene, num_proteins_in_genome, num_muts_to_sim, num_trials, sim_num_muts_per_gene, sim_prob_gene_mutated, dir_save, plot_file_name )

% Pval obs
temp_pval_obs = zeros( num_genes_mutated,1 );
temp_genes_zero_muts = setdiff( sim_gene_genenums, obs_mutated_genes_genenums );
temp_obs_num_muts_per_gene = [ obs_num_muts_per_gene, zeros( size(temp_genes_zero_muts )) ];
temp_obs_mutated_genes_genenums = [ obs_mutated_genes_genenums, temp_genes_zero_muts ];
for n=1:num_proteins_in_genome
    temp_pval_obs(n) = poisscdf( temp_obs_num_muts_per_gene(n)-1, ...
        sim_prob_gene_mutated( temp_obs_mutated_genes_genenums(n)==sim_gene_genenums )*num_muts_to_sim, ...
        'upper' );
end
% Pval sim
temp_pval_sim = zeros( num_trials,num_proteins_in_genome );
for t=1:num_trials
    for n=1:num_proteins_in_genome
        temp_pval_sim(t,n) = poisscdf( sim_num_muts_per_gene(t,n)-1, ...
            sim_prob_gene_mutated( n )*num_muts_to_sim, ...
            'upper' );
    end
end

bin_widths = [ 0.001, 0.005, 0.01, 0.02, 0.05, 0.1 ];
for i=1:numel(bin_widths)
    % Histogram pvalues
    bin_width = bin_widths(i);
    num_bins = 1/bin_width;
    bin_max = min( 1, bin_width*50 );
    num_bins_actual = bin_max/bin_width;
    temp_pval_obs_hist = histcounts(temp_pval_obs,0:bin_width:bin_max);
    temp_pval_sim_hist = zeros( num_trials,num_bins_actual );
    for t=1:num_trials
        temp_pval_sim_hist(t,:) = histcounts(temp_pval_sim(t,:),0:bin_width:bin_max);
    end
    temp_pval_sim_hist_mean = mean(temp_pval_sim_hist);
    temp_pval_sim_hist_stddev = std( temp_pval_sim_hist );
    % Figure
    figure(100); 
    clf(100)
    bar( [ temp_pval_obs_hist; temp_pval_sim_hist_mean ]' )
    xlim([0 num_bins_actual+1])
    ylim([0 1.1*max(max(temp_pval_obs_hist),max(temp_pval_sim_hist_mean))])
    xticks(0.5:5:num_bins_actual+0.5)
    xticklabels(0:5*bin_width:bin_max)
    xtickangle(90)
    xlabel( 'poiss pval' )
    ylabel('num genes')
    title([ 'all genes (n=' num2str(num_proteins_in_genome) ')' ])
    for b=1:num_bins_actual
        if b==1
            line( [b+0.125,b+0.125], [temp_pval_sim_hist_mean(b)-temp_pval_sim_hist_stddev(b) temp_pval_sim_hist_mean(b)+temp_pval_sim_hist_stddev(b)], 'Color', 'k', 'LineWidth', 1 )
        else
            line( [b+0.125,b+0.125], [temp_pval_sim_hist_mean(b)-temp_pval_sim_hist_stddev(b) temp_pval_sim_hist_mean(b)+temp_pval_sim_hist_stddev(b)], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off' )
        end
    end
    legend( 'obs', 'sim', 'stddev', 'Location', 'northeastoutside' )
    set(gca,'FontSize',20)

    % Save
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [20 8]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 20 8]);
    print([ dir_save '/' plot_file_name '_nbins-' num2str(num_bins) '.png'],'-dpng')
    
end

end