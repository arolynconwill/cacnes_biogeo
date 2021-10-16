function [ dist_to_lineage_mrca, total_N_muts_by_dist, total_S_muts_by_dist, ...
    mutgenes_genenums_bydist, mutgenes_nummuts_bydist ] = ...
    filter_snps_by_dist( this_lineage_name, dir_clusters )

%% Summary

% This function considered one lineage and estimates the age of each
% within-lineage SNP. This function returns the number of nonsynonymous (N)
% and synonymous mutations based on their age relative to the inferred
% lineage ancestor. It also returns which genes were mutated and how many
% times they were mutated.

%% Version history

% 2020.01.25, Arolyn


%% Load data

load([ dir_clusters '/' this_lineage_name '/' 'data_' this_lineage_name '.mat' ] )


%% Get dMRCAs within the lineage

if isequal( this_lineage_name,'Cluster-1-A-212' )
    hypermutators = specimen_numbers == 31 | specimen_numbers == 41;
    samples_to_include = ~hypermutators' & ~outgroup_isolates;
    Calls_lineage = Calls_for_analysis(:,samples_to_include); % exclude hypermutators and outgroup isolates
    num_samples_in_lineage = sum(samples_to_include);
    fprintf(1,['Note: excluding hypermutator clade in lineage 1...\n'])
else
    Calls_lineage = Calls_for_analysis(:,~outgroup_isolates); % exclude outgroup isolates
    num_samples_in_lineage = sum(~outgroup_isolates);
end

% Boolean showing which positions differ from the inferred lineage MRCA for each colony in that lineage
diff_mrca = ( ...
    ( Calls_lineage ~= repmat(anc_nti_goodpos,1,num_samples_in_lineage) ) ...
    & Calls_lineage > 0 ...
    );


%% Modify diff_mrca to include inferred internal nodes

% Get diff_mrca for all tree leaves
diff_mrca_leaves = diff_mrca'; % initialize and transpose
% Reduce diff_mrca to unique columns
diff_mrca_leaves = unique( diff_mrca_leaves, 'rows' );

% Create diff_mrca for shared mutations for all pairs of colonies
num_leaves = size(diff_mrca_leaves,1);
num_pos = size(diff_mrca_leaves,2);
diff_mrca_nodes = zeros( num_leaves^2, num_pos, 'logical' );
for i=1:num_leaves
    for j=1:num_leaves
        diff_mrca_nodes( i*(num_leaves-1)+j,: ) = diff_mrca_leaves(i,:) & diff_mrca_leaves(j,:);
    end
end
% Keep unique rows only
diff_mrca_nodes = unique( diff_mrca_nodes, 'rows' );

% Combine diff_mrca_leaves and diff_mrca_nodes while removing duplicate rows
diff_mrca_all = [ diff_mrca_leaves; diff_mrca_nodes ]; % combine
diff_mrca_all = unique( diff_mrca_all, 'rows' ); % remove duplicates

% Get "age" of each row
dist_to_lineage_mrca = sum( diff_mrca_all,2 ); 
% Get maximum "age"
lineage_max_dmrca = max(dist_to_lineage_mrca); 


%% Date each mutation

mut_age = zeros( num_pos,1 );
for m=1:num_pos
    has_mut = diff_mrca_all(:,m); % which rows have the mutation
    if sum(has_mut)>0
        mut_age(m) = min( dist_to_lineage_mrca( has_mut ) );
    else
        mut_age(m) = nan; % mutation might not be present if samples removed (ie was only in hypermutators)
    end
end

% Make a historgram for lineage 2
if contains( this_lineage_name, 'Cluster-2' )
    figure(10); clf(10); histogram(mut_age,0:1:max(mut_age)); xlabel('mut age'); title(this_lineage_name)
end

%% Find mutations within nnn distance of lineage MRCA and count how many of them are N or S

% Initialize
dist_to_lineage_mrca = zeros( lineage_max_dmrca,1 );
total_N_muts_by_dist = zeros( lineage_max_dmrca,1 );
total_S_muts_by_dist = zeros( lineage_max_dmrca,1 );
mutgenes_genenums_bydist = cell(lineage_max_dmrca,1);
mutgenes_nummuts_bydist = cell(lineage_max_dmrca,1);

% Loop through distances and count mutations
for next_dist = 1:lineage_max_dmrca

    % Which mutations are next_dist away from the lineage ancestor
    this_set_muts = ( mut_age == next_dist );
    num_muts = sum( this_set_muts );

    % Get mutation info
    if num_muts > 0
        these_types = [ annotation_full( this_set_muts ).type ];
        num_N_muts = sum( these_types=='N' );
        num_S_muts = sum( these_types=='S' );
        these_genes = [ annotation_full( this_set_muts ).gene_num ]; 
        these_genes = these_genes( these_genes==ceil(these_genes) ); % only coding regions
    else
        num_N_muts = 0;
        num_S_muts = 0;
        these_genes = [];
    end
    
    % Get list of unique genes and how many times each was mutated
    these_genes_list = unique( these_genes );
    these_genes_nummuts = arrayfun(@(x) sum(x==these_genes), these_genes_list );

    % Record
    dist_to_lineage_mrca(next_dist) = next_dist;
    total_N_muts_by_dist(next_dist) = num_N_muts;
    total_S_muts_by_dist(next_dist) = num_S_muts;
    mutgenes_genenums_bydist{next_dist} = these_genes_list;
    mutgenes_nummuts_bydist{next_dist} = these_genes_nummuts;

end


end