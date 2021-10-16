function num_genes_with_n_muts = bin_parallel_evo( num_muts_per_gene, num_muts_per1kb, gene_lengths, ...
    filter_min_mut_density, filter_min_gene_length, max_num_muts_per_gene )

%% Summary

% This function determines how many genes have >=N mutations (up to N = 
% max_num_muts_per_gene ), with the option to filter out genes that are too 
% short (length < filter_min_gene_length) or genes that have too low a 
% mutational density (den < filter_min_mut_density ) in muts per 1kb. 


%%

% Input data characteristics
num_sets = size( num_muts_per_gene,1 );

% Initialize
num_genes_with_n_muts = zeros( num_sets, max_num_muts_per_gene );

% Count number of genes with >= n mutations, subject to filters
for i=1:num_sets
    this_set_num_muts_per_gene = num_muts_per_gene(i,:);
    this_set_num_muts_per1kb = num_muts_per1kb(i,:);
    for n=1:max_num_muts_per_gene
        num_genes_with_n_muts(i,n) = ...
            sum( ...
            this_set_num_muts_per_gene >= n ...
            & this_set_num_muts_per1kb >= filter_min_mut_density ...
            & gene_lengths >= filter_min_gene_length ...
            );
    end
end



end