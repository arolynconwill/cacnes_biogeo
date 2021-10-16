function [ mutated_genes_nummuts, mutated_genes_genenums, mutated_genes_lengths ] = compute_observed_parallelevo( annotation_full, gene_nums_of_interest )


%% Summary

% Finds genes mutated multiple times (either by absolute number of
% mutations or by mutational density) based onn all the mutations provided
% in annotation_full

% Optional third input focuses specifically on genes of interest


%% Evaluate mutations in annotation_full

% Number of mutations
num_muts_total = length( annotation_full );

% Which mutations to consider
muts_genenums = [annotation_full.gene_num];
if nargin == 2 % consider specific genes if provided
    muts_filter = ismember( muts_genenums, gene_nums_of_interest );
else % otherwise consider all genes
    muts_filter = ceil(muts_genenums)==muts_genenums; % only keep coding regions
end

% Filter mutations
muts_genenums_filtered = muts_genenums( muts_filter );

% Get gene lengths
muts_filter_indices = find( muts_filter );
muts_genelength_filtered = zeros( sum( muts_filter ), 1 );
for g=1:numel(muts_filter_indices)
    next_index = muts_filter_indices(g);
    muts_genelength_filtered(g) = annotation_full(next_index).loc2 - annotation_full(next_index).loc1 + 1;
end


%% Tally how many times each gene is mutated

% Gene numbers and lengths of any gene mutated
[ mutated_genes_genenums, indices ] = unique( muts_genenums_filtered );
mutated_genes_lengths = muts_genelength_filtered( indices )';

% How many times each gene was mutated
mutated_genes_nummuts = arrayfun(@(x) sum(x==muts_genenums_filtered), mutated_genes_genenums );


end