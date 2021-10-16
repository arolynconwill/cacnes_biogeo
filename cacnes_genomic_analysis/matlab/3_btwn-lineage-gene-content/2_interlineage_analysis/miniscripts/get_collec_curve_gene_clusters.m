function [ num_lineages, num_of_gene_clusters ] = get_collec_curve_gene_clusters( gene_cluster_matrix )

num_lineages = size( gene_cluster_matrix, 2 );
num_permutations = 500; % could explicitly find all combinations but lazy
num_of_gene_clusters_trials = zeros( num_lineages, num_permutations );

for i=1:num_lineages
    for j=1:num_permutations
        num_of_gene_clusters_trials( i,j ) = sum( sum( gene_cluster_matrix(:,randperm(num_lineages,i)),2 )>0 );
    end
end

num_of_gene_clusters = mean( num_of_gene_clusters_trials,2 );

end