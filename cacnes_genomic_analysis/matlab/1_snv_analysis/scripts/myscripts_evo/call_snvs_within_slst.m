function my_set_names = call_snvs_within_slst( cluster_names, clusters_all, dir_clusters, dir_ref_genome, dir_save_strain )

%% Summary

% This function uses lineage MRCAs and lineage SLST assignments to look at
% within-superSLST evolution by calling SNPs that distinguish lineage MRCAs
% from each other within the same superSLST.


%% Directory setup

if ~exist( dir_save_strain, 'dir')
    mkdir( dir_save_strain )
end

%% Load additional data

% Load SLST for each cluster
load( [ dir_clusters '/' 'sample_names.mat' ], 'SampleNamesLong_all' )
slst_all = cellfun(@(x) x(end-6:end-5), SampleNamesLong_all, 'UniformOutput', false);
clusters_all_slsts = {};
for i=1:numel(clusters_all)
    this_cluster_slsts = unique(slst_all(clusters_all{i}));
    if numel(this_cluster_slsts)>1
        clusters_all_slsts{end+1} = this_cluster_slsts{1};
    else
        clusters_all_slsts{end+1} = this_cluster_slsts{1};
    end
end
clusters_all_slsts_coarse = cellfun(@(x) x(1), clusters_all_slsts);

% Load MRCA for each cluster 
load( [ '../data_snps/' 'cluster_step_variables.mat' ], 'p_all' )
Calls_clades_all = zeros( numel(clusters_all), numel(p_all) );
for c=1:numel(cluster_names)
    load( [ dir_clusters '/' cluster_names{c} '/' 'data_' cluster_names{c} '.mat' ], 'Calls_lineage_ancestor' )
    Calls_clades_all(c,:) = Calls_lineage_ancestor';
end
% Save
save([ dir_save_strain '/' 'lineage_MRCAs' '.mat'], ...
    'Calls_clades_all', 'cluster_names' ) % save this for comparison later


% Load reference genome info
NTs = 'ATCG'; % Nucleotides: 1=A, 2=T, 3=C, 4=G
[ChrStarts, GenomeLength, ~, ScafNames] = genomestats(dir_ref_genome);
refnt_all = extract_outgroup_mutation_positions(dir_ref_genome, p2chrpos(p_all,ChrStarts));
[~,refnti_all]=ismember(refnt_all,NTs); 


%% Find and annotate mutations within each set

% Group lineages by superSLST
slst_list = unique( clusters_all_slsts_coarse );
slst_sets = arrayfun(@(x) find(x==clusters_all_slsts_coarse), slst_list, 'UniformOutput', false);
% Only keep SLSTs with more than one clade
sets_keep = ( cellfun(@(x) numel(x), slst_sets ) > 1 );
slsts_keep = slst_list( sets_keep );
slst_sets_keep = slst_sets( sets_keep );

% Call SNPs within superSLST
my_set_names = {};
for this_set = 1:numel(slst_sets_keep) 

    % Get info for this set
    this_set_name = [ 'SLST-' slsts_keep(this_set) ];
    my_set_names{end+1} = this_set_name;
    fprintf(1,[ 'Analyzing ' this_set_name ' clades...\n' ])
    this_set_indices = slst_sets_keep{this_set};
    Calls_this_set = Calls_clades_all( this_set_indices,: ); %Nx239547; includes all positions in p_all

    % Define "ancestor" based on outgroup (using all lineage MRCAs from other superSLSTs; not super important because we only want N/S and don't care about mutational direction for this analysis)
    outgroup_isolates = setdiff(1:1:numel(cluster_names),this_set_indices); 
    outgroup_calls = Calls_clades_all( outgroup_isolates,: );
    outgroup_calls( outgroup_calls==0 ) = nan; % mode takes the lowest value when there is a tie, setting to nan avoids this
    anc_nti = mode(outgroup_calls,1); % set ancestor nucleotides to mode of outgroup sample calls at that position
    % Replace with reference if no value exists from outgroup
    refnti_clades = refnti_all;
    anc_nti(isnan(anc_nti)) = refnti_clades(isnan(anc_nti)); 

    % How many positions actually vary WITHIN the set??
    all_A = sum( Calls_this_set == 1, 1) == sum( Calls_this_set ~=0, 1 );
    all_T = sum( Calls_this_set == 2, 1) == sum( Calls_this_set ~=0, 1 );
    all_C = sum( Calls_this_set == 3, 1) == sum( Calls_this_set ~=0, 1 );
    all_G = sum( Calls_this_set == 4, 1) == sum( Calls_this_set ~=0, 1 );
    all_same = all_A + all_T + all_C + all_G;
    
    % Preliminary SNP positions
    goodpos_preliminary = find(~all_same); % only want variation WITHIN group of clades with the same SLST assignment
    fprintf(1,['Preliminary number of SNP positions: ' num2str(numel(goodpos_preliminary)) '\n']);

    % Filter out positions with suspected recombination
    fixedmutation_preliminary = ( (Calls_this_set(:,goodpos_preliminary)~=repmat(anc_nti(goodpos_preliminary),numel(this_set_indices),1)) & Calls_this_set(:,goodpos_preliminary)>0 ); 
    save_plots = true;
    dir_figs = dir_save_strain;
    % Parameters for filtering recombination regions
    recombination_block_size=150;
    correlation_cutoff=0.75;
    p_involved_in_non_snp_event = identify_non_snp_events_strain( ...
        recombination_block_size, correlation_cutoff, ...
        goodpos_preliminary, p_all, fixedmutation_preliminary', GenomeLength, ...
        save_plots, dir_figs, this_set_name );
    % Update goodpos by removing positions with suspected recombination
    goodpos = setdiff( goodpos_preliminary, find(ismember(p_all,p_involved_in_non_snp_event)) );
    p_all_set = p_all( goodpos ); 
    fprintf(1,['Number of SNP positions after filtering suspected recombination: ' num2str(numel(goodpos)) '\n']);
    
    % Lineage-specific info
    p4annotations = p_all;
    NTs4annotations = anc_nti;
    % Clade data
    anc_nti_goodpos = anc_nti(goodpos); 
    Calls_this_set_goodpos = Calls_this_set(:,goodpos); % Calls
    counts_this_set = ones( [8 numel(goodpos) numel(this_set_indices) ] ); % counts placeholder NOT REAL COUNTS, JUST NEED SOMETHING FOR ANNOTATION FUNCTIONS TO RUN!!!
    fixedmutation = ( ( Calls_this_set_goodpos~=repmat(anc_nti_goodpos,numel(this_set_indices),1) ) & Calls_this_set_goodpos>0 )'; % fixedmutation placeholder
    % Remove recombination regions from fake fixedmutation
    recombinationregion_L11 = find((p_all_set>=1793873 & p_all_set<=1793881) & sum(fixedmutation,2)>0); 
    recombinationregion_prophage = find((p_all_set>=1398569 & p_all_set<=1418224) & sum(fixedmutation,2)>0); 
    if ~isempty(recombinationregion_L11)
        fixedmutation(recombinationregion_L11,:)=0;
        fprintf(1,['Removed L11 positions. \n'])
    end
    if ~isempty(recombinationregion_prophage)
        fixedmutation(recombinationregion_prophage,:)=0;
        fprintf(1,['Removed prophage positions. \n'])
    end

    % Make annotation table
    promotersize = 150; % bp
    annotations = annotate_mutations_gb_lineage(p2chrpos(p_all(goodpos),ChrStarts),dir_ref_genome,p4annotations, NTs4annotations) ;
    annotation_full = append_annotations(annotations, anc_nti_goodpos, Calls_this_set_goodpos', counts_this_set, fixedmutation, promotersize) ; 

    % Get mutation info
    directory_save = 'none';
    save_data = false;
    [mutationmatrix, mut_observed, typecounts, prob_nonsyn ] = ...
        mutation_spectrum_module( goodpos, anc_nti_goodpos, ...
        Calls_this_set_goodpos', annotation_full, ...
        p_all, NTs, ...
        directory_save, save_data );

    % Save
    save([ dir_save_strain '/' this_set_name '_mutations' '.mat'], ...
        'mutationmatrix', 'mut_observed', 'typecounts', 'prob_nonsyn', ...
        'annotation_full', 'Calls_this_set', ...
        'this_set_name', 'this_set_indices' ) % save this for comparison later

end


end