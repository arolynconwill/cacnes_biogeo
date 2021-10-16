function [ Nobs_by_geneset, Sobs_by_geneset, dNdS_by_geneset, dNdS_lower_by_geneset, dNdS_upper_by_geneset, dnds_by_geneset, pval_dNdS ] = ...
    get_dNdS_for_genesets_of_interest( num_genesets_of_interest, genesets_of_interest_genenums, annotation_full_everything, ...
    dir_ref_genome, mut_spec_prob, dir_save ) 


%% Initialize

Nobs_by_geneset  = zeros( 1,num_genesets_of_interest ); % initialize
Sobs_by_geneset  = zeros( 1,num_genesets_of_interest ); % initialize
dnds_by_geneset = cell( num_genesets_of_interest,1 ); % initialize
dNdS_by_geneset  = zeros( 1,num_genesets_of_interest ); % initialize
dNdS_lower_by_geneset = zeros( 1,num_genesets_of_interest ); % initialize
dNdS_upper_by_geneset = zeros( 1,num_genesets_of_interest ); % initialize


%% dN/dS for potential genes of interest

% Loop through genes of interest
for n=1:numel(genesets_of_interest_genenums)
    
    % Get gene number
    my_genenums = genesets_of_interest_genenums{n};
    
    % Compute observed N/S
    [ p_nonsyn, CI_nonsyn, num_muts_N, num_muts_S ] = ...
        compute_observed_dnds( annotation_full_everything, my_genenums );
    dnds_observed = (p_nonsyn/(1-p_nonsyn)); % N/S observed
    fprintf(1,[ 'N: ' num2str(num_muts_N) '\n' ])
    fprintf(1,[ 'S: ' num2str(num_muts_S) '\n' ])
    
    % Compute expected N/S
    probnonsyn_expected = compute_expected_dnds( dir_ref_genome, mut_spec_prob, my_genenums );
    dnds_expected_local = probnonsyn_expected/(1-probnonsyn_expected); % N/S expected from neutral model
    
    % Compute dN/dS
    dNdS = dnds_observed/dnds_expected_local; % relative to neutral model
    CI_lower = ( CI_nonsyn(1)/(1-CI_nonsyn(1)) ) / dnds_expected_local;
    CI_upper = ( CI_nonsyn(2)/(1-CI_nonsyn(2)) ) / dnds_expected_local;
    
    % Record
    Nobs_by_geneset(n) = num_muts_N;
    Sobs_by_geneset(n) = num_muts_S;
    dNdS_by_geneset(n)  = dNdS;
    dNdS_lower_by_geneset(n) = CI_lower;
    dNdS_upper_by_geneset(n) = CI_upper;
    dnds_by_geneset{n} = [ num2str(num_muts_N) '/' num2str(num_muts_S) ]; 
    
    fprintf(1,'\n')
    
end


end