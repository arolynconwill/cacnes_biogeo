function [ Nobs_by_gene, Sobs_by_gene, dNdS_by_gene, dNdS_lower_by_gene, dNdS_upper_by_gene, name_by_gene, pval_dNdS, locustag_by_gene ] = ...
    get_dNdS_for_genes_of_interest( num_genes_of_interest, genes_of_interest_genenums, annotation_full_everything, ...
    dir_ref_genome, mut_spec_prob, dir_save ) 


%% Initialize

Nobs_by_gene  = zeros( 1,num_genes_of_interest+1 ); % initialize
Sobs_by_gene  = zeros( 1,num_genes_of_interest+1 ); % initialize
dNdS_by_gene  = zeros( 1,num_genes_of_interest+1 ); % initialize
dNdS_lower_by_gene = zeros( 1,num_genes_of_interest+1 ); % initialize
dNdS_upper_by_gene = zeros( 1,num_genes_of_interest+1 ); % initialize
name_by_gene = cell( num_genes_of_interest+1,1 ); % initialize
locustag_by_gene = cell( num_genes_of_interest+1,1 ); % initialize

%% dN/dS for potential genes of interest

% Loop through genes of interest
for n=1:numel(genes_of_interest_genenums)
    
    % Get gene number
    my_genenum = genes_of_interest_genenums(n);
    % Find mutation entries in table
    my_indices = find( my_genenum==[annotation_full_everything.gene_num] );
    [~,my_indices_resort] = sort( [ annotation_full_everything(my_indices).pos ] );
    annotation_full_mini = annotation_full_everything(my_indices(my_indices_resort));
    save([ dir_save '/' 'gene_' num2str(my_genenum) '.mat' ], 'annotation_full_mini' )
    % Which protein
    my_protein = strtrim( [ '(' num2str(my_genenum) ') ' annotation_full_everything(my_indices(1)).protein ] );
    my_locustag = annotation_full_everything(my_indices(1)).locustag;
    fprintf(1,[ my_locustag '\n' ] )
    fprintf(1,[ my_protein '\n' ] )
    % Gene length
    my_length = annotation_full_everything(my_indices(1)).loc2 - annotation_full_everything(my_indices(1)).loc1;
    fprintf(1,[ 'gene_length: ' num2str( my_length ) '\n' ])
    
    % Compute observed N/S
    [ p_nonsyn, CI_nonsyn, num_muts_N, num_muts_S ] = ...
        compute_observed_dnds( annotation_full_everything, my_genenum );
    dnds_observed = (p_nonsyn/(1-p_nonsyn)); % N/S observed
    fprintf(1,[ 'N: ' num2str(num_muts_N) '\n' ])
    fprintf(1,[ 'S: ' num2str(num_muts_S) '\n' ])
    
    % Compute expected N/S
    probnonsyn_expected = compute_expected_dnds( dir_ref_genome, mut_spec_prob, my_genenum );
    dnds_expected_local = probnonsyn_expected/(1-probnonsyn_expected); % N/S expected from neutral model
    
    % Compute dN/dS
    dNdS = dnds_observed/dnds_expected_local; % relative to neutral model
    CI_lower = ( CI_nonsyn(1)/(1-CI_nonsyn(1)) ) / dnds_expected_local;
    CI_upper = ( CI_nonsyn(2)/(1-CI_nonsyn(2)) ) / dnds_expected_local;
    
    % Record
    Nobs_by_gene(n) = num_muts_N;
    Sobs_by_gene(n) = num_muts_S;
    dNdS_by_gene(n)  = dNdS;
    dNdS_lower_by_gene(n) = CI_lower;
    dNdS_upper_by_gene(n) = CI_upper;
    name_by_gene{n} = my_protein; % num2str(my_genenum); 
    locustag_by_gene{n} = my_locustag; 
    dnds_by_gene{n} = [ num2str(num_muts_N) '/' num2str(num_muts_S) ]; 
    
    fprintf(1,'\n')
    
end

%% dN/dS across all potential genes of interest

% Compute observed N/S
[ p_nonsyn_all, CI_nonsyn_all, num_muts_N_all, num_muts_S_all ] = ...
    compute_observed_dnds( annotation_full_everything, genes_of_interest_genenums );
dnds_observed_all = (p_nonsyn_all/(1-p_nonsyn_all)); % N/S observed; equivalent to num_muts_N/num_muts_S

% Compute expected N/S
probnonsyn_expected = compute_expected_dnds( dir_ref_genome, mut_spec_prob, genes_of_interest_genenums );
dnds_expected_local = probnonsyn_expected/(1-probnonsyn_expected); % N/S expected from neutral model

% Compute dN/dS with 95% CIs
dNdS_all = dnds_observed_all/dnds_expected_local; % relative to neutral model
CI_lower_all = ( CI_nonsyn_all(1)/(1-CI_nonsyn_all(1)) ) / dnds_expected_local;
CI_upper_all = ( CI_nonsyn_all(2)/(1-CI_nonsyn_all(2)) ) / dnds_expected_local;

% Binomial test: pval for Nobs and Sobs given probNexp
pval_dNdS = binocdf( num_muts_N_all-1, num_muts_N_all+num_muts_S_all, probnonsyn_expected, 'upper');

% Record info across all genes
Nobs_by_gene(num_genes_of_interest+1) = num_muts_N_all;
Sobs_by_gene(num_genes_of_interest+1) = num_muts_S_all;
dNdS_by_gene(num_genes_of_interest+1) = dNdS_all;
dNdS_lower_by_gene(num_genes_of_interest+1) = CI_lower_all;
dNdS_upper_by_gene(num_genes_of_interest+1) = CI_upper_all;
name_by_gene{num_genes_of_interest+1} = 'all genes of interest'; 
locustag_by_gene{num_genes_of_interest+1} = 'all'; 


end