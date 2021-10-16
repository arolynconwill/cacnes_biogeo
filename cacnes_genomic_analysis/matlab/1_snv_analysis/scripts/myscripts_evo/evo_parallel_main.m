function [ dNdS_return, CIs_return, num_genes_detected, dNdS_top_genes, CIs_top_genes ] ...
    = evo_parallel_main( data_set_name, dir_results, annotation_full_everything, ...
    dir_ref_genome, GenomeLength, positions_to_mask, mut_spec_prob, ...
    fdr_to_test, top_num_genes_to_test )

fprintf(1,['Analyzing data set: ' data_set_name '...\n'])


%% DIRECTORY SETUP
fprintf(1,['Setting up directories...' '\n'])

% Where to save results
if ~exist( dir_results, 'dir')
    mkdir( dir_results )
end

dir_geneinfo = [ dir_results '/' 'gene_info' ];
if ~exist( dir_geneinfo, 'dir')
    mkdir( dir_geneinfo )
end



%% I. OBSERVED and EXPECTED NUMBER OF GENES WITH MULTIPLE MUTATIONS


%% Get gene numbers for proteins only (exclude tRNA and rRNA)
fprintf(1,['Getting genome information...' '\n'])

cdsindex_proteins = get_coding_genenums( dir_ref_genome, positions_to_mask ); % avoid genes in regions masked in SNP calling
num_proteins_in_genome = numel( cdsindex_proteins );


%% Observed genes with multiple mutations
fprintf(1,['Finding mutated genes...' '\n'])

[ obs_num_muts_per_gene, obs_mutated_genes_genenums, obs_genes_lengths ] = ...
    compute_observed_parallelevo( annotation_full_everything, cdsindex_proteins ); % sim_gene_genenums used to exclude rRNA, tRNA
obs_num_mutsper1kb_per_gene = 1000*obs_num_muts_per_gene./obs_genes_lengths;
num_genes_mutated = numel( obs_mutated_genes_genenums );
obs_num_muts_total = sum( obs_num_muts_per_gene );


%% Expected genes with multiple mutations (simulate neutral model)
fprintf(1,['Simulating mutated genes...' '\n'])

num_muts_to_sim = obs_num_muts_total;
num_trials = 1000;
[ sim_num_muts_per_gene, sim_num_mutsper1kb_per_gene, sim_gene_lengths, sim_gene_genenums, sim_prob_gene_mutated ] = ...
    compute_expected_parallelevo( dir_ref_genome, mut_spec_prob, num_muts_to_sim, num_trials, cdsindex_proteins );


%%

% Histogram of gene lengths
% figure(10); histogram(sim_gene_lengths); xlabel('gene lengths'); ylabel('num genes'); title('all proteins'); set(gca,'FontSize',20)

% if isequal( data_set_name, 'all-muts' )
%     dir_save = dir_results; 
%     plot_file_name = [ 'parallel_' data_set_name '_doublebar-pval-obs-sim' ];
%     plot_compare_pvals_obs_sim( num_genes_mutated, sim_gene_genenums, obs_mutated_genes_genenums, ...
%         obs_num_muts_per_gene, num_proteins_in_genome, num_muts_to_sim, num_trials, sim_num_muts_per_gene, sim_prob_gene_mutated, dir_save, plot_file_name )
% end


%% COMPARE OBSERVED VS EXPECTED NUMBERS OF MULTIPLY MUTATED GENES
fprintf(1,['Computing dN/dS as a function of mutational density...' '\n'])

% Filters
min_num_muts_per_gene = 2; % only do density for multiply mutated genes
list_density_thresholds = 1:1:5; 

% Number of gene sets to investigate
num_genesets_of_interest = numel(list_density_thresholds)+1;

% % Initialize
genesets_of_interest_genenums = cell( numel(list_density_thresholds)+1,1 );

% Find genes of interest
for i=1:numel(list_density_thresholds)
    this_density_threshold = list_density_thresholds(i);
    genes_of_interest_bool = ( ...
        obs_num_muts_per_gene >= min_num_muts_per_gene ...
        & obs_num_mutsper1kb_per_gene >= this_density_threshold ...
        );
    genesets_of_interest_genenums{i} = obs_mutated_genes_genenums( genes_of_interest_bool );
end
genesets_of_interest_genenums{end} = obs_mutated_genes_genenums; % all genes

% Get dN/dS info for genes of interest
[ ~, ~, dNdS_by_geneset, dNdS_lower_by_geneset, dNdS_upper_by_geneset, ~ ] = ...
    get_dNdS_for_genesets_of_interest( num_genesets_of_interest, genesets_of_interest_genenums, annotation_full_everything, ...
    dir_ref_genome, mut_spec_prob, dir_geneinfo );

% Figure: dN/dS for genes of interest
inf_substitute_for_plotting = 100; % since we can't plot infinity
dNdS_for_plot = dNdS_by_geneset;
dNdS_for_plot( dNdS_for_plot==Inf ) = inf_substitute_for_plotting;
CI_for_plot = [ dNdS_lower_by_geneset; dNdS_upper_by_geneset ];
CI_for_plot( CI_for_plot==Inf ) = inf_substitute_for_plotting;
plot_title = 'dN/dS by mutational density';
x_axis_label = 'mutational density threshold (mut/kb)';
density_list_labels = arrayfun(@(x) { num2str(x) }, list_density_thresholds );
density_list_labels{end+1} = 'all';
x_tick_labels = density_list_labels;
dir_save = dir_results ; 
plot_file_name = [ 'parallel-sets_' data_set_name '_dNdS_density' ];
y_axis_factor = 4;
% Call plot function
plot_dnds_bar( dNdS_for_plot, CI_for_plot, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
    dir_save, plot_file_name ) 

% Save data
save( [ dir_results '/' 'data_dNdS_by_density.mat' ], ...
    'min_num_muts_per_gene', 'list_density_thresholds', ...
    'dNdS_by_geneset', 'dNdS_lower_by_geneset', 'dNdS_upper_by_geneset' )


%% IDENTIFY GENES OF INTEREST USING POISSON STATISTICS AND BENJAMINI-HOCHBERG PROCEDURE


% Poisson pvalue for the observed number of mutations on each gene
pval_poisson = zeros( num_genes_mutated,1 );
for n=1:num_genes_mutated
    pval_poisson(n) = poisscdf( obs_num_muts_per_gene(n)-1, ...
        sim_prob_gene_mutated( obs_mutated_genes_genenums(n)==sim_gene_genenums )*num_muts_to_sim, ...
        'upper' ); % poisscdf( x, labmda, 'upper' ) computes pval for >x not >=x
end


%% For supplemental figure: 
% % number of genes detected for different FDR thresholds
% % dN/dS for top genes

% Non-changing params
min_num_muts_per_gene = 1; % minimum number of mutations on a gene

% Number of genes detected at different FDR's 
num_genes_detected = zeros( size( fdr_to_test ) );

for i=1:numel(fdr_to_test)

    % Parameters
    fdr_bh = fdr_to_test(i); % false discovery rate
    % Strings to record params in fig names
    fdr_str = [ '_fdr-' num2str(100*fdr_bh) ];
    nmuts_str = [ '_minmuts-' num2str(min_num_muts_per_gene) ];
    params_str = [ fdr_str nmuts_str ];
    fprintf(1,['Params: ' params_str '\n'])

    % Benjamini-Hochberg procedure
    plot_file_name = [ 'parallel_' data_set_name '_benjhoch' fdr_str nmuts_str ];
    dir_save = dir_results; 
    [ ~, pval_max ] = benjamini_hochberg_sig_filter( pval_poisson, num_proteins_in_genome, fdr_bh, dir_save, plot_file_name );

    % Genes of interest
    genes_of_interest_bool = ( pval_poisson'<=pval_max & obs_num_muts_per_gene>=min_num_muts_per_gene );
    num_genes_of_interest = sum( genes_of_interest_bool );

    num_genes_detected(i) = num_genes_of_interest;
        
end

% dN/dS for top N genes

[ ~, pval_order ] = sort( pval_poisson, 'ascend' );
dNdS_top_genes = zeros( numel(top_num_genes_to_test),1 ); % initialize
CIs_top_genes = zeros( numel(top_num_genes_to_test),2 ); % initialize
genenums_to_save = cell( numel(top_num_genes_to_test),1 );
for n=1:numel(top_num_genes_to_test)
    top_num_genes_to_test_num = top_num_genes_to_test(n);
    top_gene_indices = pval_order( 1:top_num_genes_to_test_num );
    % Genes of interest
    genes_of_interest_genenums = obs_mutated_genes_genenums( top_gene_indices );
    genenums_to_save{n} = genes_of_interest_genenums;
    genes_of_interest_lengths = obs_genes_lengths( top_gene_indices );
    num_genes_of_interest = top_num_genes_to_test_num;
    % Get dN/dS for these genes
    [ ~, ~, dNdS_by_gene, dNdS_lower_by_gene, dNdS_upper_by_gene, ~, ~ ] = ...
        get_dNdS_for_genes_of_interest( num_genes_of_interest, genes_of_interest_genenums, annotation_full_everything, ...
        dir_ref_genome, mut_spec_prob, dir_geneinfo );
    dNdS_top_genes(n) = dNdS_by_gene(end);
    CIs_top_genes(n,:) = [ dNdS_lower_by_gene(end); dNdS_upper_by_gene(end) ];
end
% Save data for supp
save( [ dir_results '/' 'data_for_supp.mat' ], ...
    'fdr_to_test', 'min_num_muts_per_gene', 'top_num_genes_to_test', ...
    'num_genes_detected', 'dNdS_top_genes', 'CIs_top_genes', 'genenums_to_save' )


%% Broad scan through FDRs to identify genes of interest
fprintf(1,['Identifying genes of interest...' '\n'])

fdr_bh_list = [ 0.1 0.2 0.5 1 2 5 1000 ];
min_num_muts_list = [ 2 ];

for i=1:numel(fdr_bh_list)
    for j=1:numel(min_num_muts_list)
        
        % Parameters
        fdr_bh = fdr_bh_list(i); % false discovery rate
        min_num_muts_per_gene = min_num_muts_list(j); % minimum number of mutations on a gene
        % Strings to record params in fig names
        fdr_str = [ '_fdr-' num2str(100*fdr_bh) ];
        nmuts_str = [ '_minmuts-' num2str(min_num_muts_per_gene) ];
        params_str = [ fdr_str nmuts_str ];
        fprintf(1,['Params: ' params_str '\n'])

        % Benjamini-Hochberg procedure
        plot_file_name = [ 'parallel_' data_set_name '_benjhoch' fdr_str ];
        dir_save = dir_results; 
        [ ~, pval_max ] = benjamini_hochberg_sig_filter( pval_poisson, num_proteins_in_genome, fdr_bh, dir_save, plot_file_name );

        % Genes of interest
        genes_of_interest_bool = ( pval_poisson'<=pval_max & obs_num_muts_per_gene>=min_num_muts_per_gene );
        num_genes_of_interest = sum( genes_of_interest_bool );

        genes_of_interest_genenums = obs_mutated_genes_genenums( genes_of_interest_bool );
        genes_of_interest_lengths = obs_genes_lengths( genes_of_interest_bool );

        % Plots describing selection of genes of interest
        % Extra plot: scatter plot of poisson pvalue vs gene length with genes of interest highlighted
        plot_title = 'all within-lineage mutated genes';
        dir_save = dir_results; 
        plot_file_name = [ 'parallel_' data_set_name '_scatter' params_str ];
        compare_filter_to_poisson( obs_genes_lengths, pval_poisson, obs_mutated_genes_genenums, genes_of_interest_genenums, ...
            plot_title, dir_save, plot_file_name );
        
        if num_genes_of_interest>0

            % Extra plot: histogram of poisson pvalues with genes of interest highlighted
            plot_title = 'all within-lineage mutated genes';
            x_axis_label = 'poisson pvalue';
            y_axis_label = 'num of genes';
            line_label = [ 'pval_cutoff=' num2str(pval_max) ' (BH)' ];
            plot_file_name = [ 'parallel_' data_set_name '_hist' params_str ];
            plot_parallelevo_hist_exp_line_obs_poisson( ...
                log10(pval_max), log10(pval_poisson), genes_of_interest_bool, ...
                plot_title, x_axis_label, y_axis_label, line_label, dir_results, plot_file_name );

            if num_genes_of_interest<=50

                % Get dN/dS info for genes of interest
                [ Nobs_by_gene, Sobs_by_gene, dNdS_by_gene, dNdS_lower_by_gene, dNdS_upper_by_gene, name_by_gene, pval_dNdS ] = ...
                    get_dNdS_for_genes_of_interest( num_genes_of_interest, genes_of_interest_genenums, annotation_full_everything, ...
                    dir_ref_genome, mut_spec_prob, dir_geneinfo );
                length_by_gene = [ genes_of_interest_lengths, mean(genes_of_interest_lengths) ];

                % Figure: dN/dS for genes of interest
                inf_substitute_for_plotting = 100; % since we can't plot infinity
                dNdS_for_plot = dNdS_by_gene;
                dNdS_for_plot( dNdS_for_plot==Inf ) = inf_substitute_for_plotting;
                CI_for_plot = [ dNdS_lower_by_gene; dNdS_upper_by_gene ];
                CI_for_plot( CI_for_plot==Inf ) = inf_substitute_for_plotting;
                plot_title = [ 'genes of interest: ' data_set_name ' (p=' num2str(pval_dNdS) ')' ];
                x_axis_label = 'gene annotation';
                x_tick_labels = name_by_gene;
                dir_save = dir_results ; 
                plot_file_name = [ 'parallel_' data_set_name '_dNdS' params_str  ];
                text_labels_1 = arrayfun(@(i) [ num2str(Nobs_by_gene(i)) '/' num2str(Sobs_by_gene(i)) ], 1:1:num_genes_of_interest+1, 'UniformOutput', false);
                text_labels_2 = length_by_gene;
                y_axis_factor = 8;
                % Call plot function
                plot_dnds_bar( dNdS_for_plot, CI_for_plot, ...
                    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
                    dir_save, plot_file_name, ...
                    text_labels_1, text_labels_2 ) 

            end
            
        end

    end
end


%% Broad scan through FDRs to compute dN/dS
fprintf(1,['Evaluating dN/dS...' '\n'])

fdr_bh_list = [ 0.1 0.2 0.5 1 2 5 1000 ];
min_num_muts_list = [ 1 2 ];

% Save dN/dS for all mutations meeting min_num_muts_list criteria
dNdS_return = size( min_num_muts_list );
CIs_return = zeros( numel(min_num_muts_list), 2);

for i=1:numel( min_num_muts_list )
    
    % Identify genes of interest
    num_genesets_of_interest = numel( fdr_bh_list );
    genesets_of_interest_genenums = cell( num_genesets_of_interest,1 );
    genesets_of_interest_size = zeros( num_genesets_of_interest,1 );
    % Loop through BH FDR cutoffs
    for j=1:numel( fdr_bh_list )
        [ ~, pval_max ] = benjamini_hochberg_sig_filter( pval_poisson, num_proteins_in_genome, fdr_bh_list(j) );
        genes_of_interest_bool = ( pval_poisson'<=pval_max & obs_num_muts_per_gene>=min_num_muts_list(i) );
        genes_of_interest_genenums = obs_mutated_genes_genenums( genes_of_interest_bool );
        genesets_of_interest_genenums{j} = genes_of_interest_genenums;
        genesets_of_interest_size(j) = sum( genes_of_interest_bool );
    end

    % Get dN/dS info for genes of interest
    [ ~, ~, dNdS_by_geneset, dNdS_lower_by_geneset, dNdS_upper_by_geneset, dnds_by_geneset ] = ...
        get_dNdS_for_genesets_of_interest( num_genesets_of_interest, genesets_of_interest_genenums, annotation_full_everything, ...
        dir_ref_genome, mut_spec_prob, dir_geneinfo );

    % Figure: dN/dS for genes of interest
    inf_substitute_for_plotting = 100; % since we can't plot infinity
    dNdS_for_plot = dNdS_by_geneset;
    dNdS_for_plot( dNdS_for_plot==Inf ) = inf_substitute_for_plotting;
    CI_for_plot = [ dNdS_lower_by_geneset; dNdS_upper_by_geneset ];
    CI_for_plot( CI_for_plot==Inf ) = inf_substitute_for_plotting;
    plot_title = 'genes of interest: dN/dS';
    x_axis_label = 'BH FDR';
    fdr_bh_list_labels = arrayfun(@(x) { num2str(x) }, fdr_bh_list );
    fdr_bh_list_labels{end} = 'all';
    x_tick_labels = fdr_bh_list_labels;
    dir_save = dir_results ; 
    plot_file_name = [ 'parallel-sets_' data_set_name '_dNdS_nmuts-' num2str(min_num_muts_list(i)) ];
    text_labels_1 = dnds_by_geneset;
    text_labels_2 = genesets_of_interest_size;
    y_axis_factor = 4;
    % Call plot function
    plot_dnds_bar( dNdS_for_plot, CI_for_plot, ...
        plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
        dir_save, plot_file_name, ...
        text_labels_1, text_labels_2 ) 
    
    % Save
    dNdS_return(i) = dNdS_for_plot(end);
    CIs_return(i,:) = CI_for_plot(:,end);
    
end


%% COMPARE OBSERVED VS EXPECTED NUMBERS OF MULTIPLY MUTATED GENES
fprintf(1,['Comparing observed vs expected numbers of multiply mutated genes...' '\n'])

% Filter parameters
filter_min_gene_length = 0;
min_density_list = 1:1:5; % filter_min_mut_density_factor
max_num_muts_per_gene = max( max(max(sim_num_muts_per_gene)), max(obs_num_muts_per_gene) ); 

% Initialize
record_obs_num_genes_with_n_muts = zeros( numel(min_density_list), max_num_muts_per_gene );
record_sim_num_genes_with_n_muts_mean = zeros( numel(min_density_list), max_num_muts_per_gene );
record_sim_num_genes_with_n_muts_std = zeros( numel(min_density_list), max_num_muts_per_gene );

for filter_index = 1:numel(min_density_list)

    %% Compare how many genes were mutated N times (observed vs expected)

    % Filters
    filter_min_mut_density_factor = min_density_list(filter_index); % what multiple of average mutational density is required for a gene of interest
    filter_min_mut_density = filter_min_mut_density_factor*1000*obs_num_muts_total/GenomeLength; % multiple of the observed average mutation density (muts per 1kb)

    %% Number of genes with N mutations

    % Observed
    obs_num_genes_with_n_muts = bin_parallel_evo( obs_num_muts_per_gene, obs_num_mutsper1kb_per_gene, obs_genes_lengths, ...
        filter_min_mut_density, filter_min_gene_length, max_num_muts_per_gene );
    record_obs_num_genes_with_n_muts(filter_index,:) = obs_num_genes_with_n_muts;

    % Expected
    sim_num_genes_with_n_muts = bin_parallel_evo( sim_num_muts_per_gene, sim_num_mutsper1kb_per_gene, sim_gene_lengths, ...
        filter_min_mut_density, filter_min_gene_length, max_num_muts_per_gene );
    sim_num_genes_with_n_muts_mean = mean( sim_num_genes_with_n_muts );
    sim_num_genes_with_n_muts_stddev = std( sim_num_genes_with_n_muts );
    record_sim_num_genes_with_n_muts_mean(filter_index,:) = sim_num_genes_with_n_muts_mean;
    record_sim_num_genes_with_n_muts_std(filter_index,:) = sim_num_genes_with_n_muts_stddev;


    %% Figure
    
    % Figure
    plot_title = ['all within-set SNVs (den_factor=' num2str(filter_min_mut_density_factor) ')'];
    x_axis_label = '\geq n mutations';
    plot_file_name = [ 'parallel_' data_set_name '_nummuts_all-muts_doublebar_minden-' num2str(filter_min_mut_density_factor) ];
    make_inset = true; inset_start = ceil(numel(obs_num_genes_with_n_muts)/2);
    plot_parallelevo_doublebar( obs_num_genes_with_n_muts, sim_num_genes_with_n_muts_mean, sim_num_genes_with_n_muts_stddev, plot_title, x_axis_label, dir_results, plot_file_name, make_inset, inset_start ) 

end

% Save data
save( [ dir_results '/' 'data_old_method.mat' ], ...
    'filter_min_gene_length', 'min_density_list', 'max_num_muts_per_gene', ...
    'record_obs_num_genes_with_n_muts', 'record_sim_num_genes_with_n_muts_mean', 'record_sim_num_genes_with_n_muts_std' )


end
