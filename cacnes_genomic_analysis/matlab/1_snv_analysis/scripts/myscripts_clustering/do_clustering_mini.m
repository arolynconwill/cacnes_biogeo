function trial_clusters = do_clustering_mini( trial_distance_matrix, dist_cutoff, clustering_minpts, allow_cluster_merging, threshold_merging )

dm=trial_distance_matrix; % distance matrix
dm_i=max(dm(:))-dm; % similarity matrix
clusters=dbscan(dm_i,max(dm_i(:))-dist_cutoff,clustering_minpts);
maxclusters=max(clusters);
% Merge clusters that are close
if allow_cluster_merging % only merge if merging is on
    for i=1:maxclusters          
        for j=(i+1):maxclusters
            distance_ij=dm(clusters==i,clusters==j);
            if mean(distance_ij(:)) < threshold_merging
                clusters(clusters==j)=i; 
            end
        end
    end
end
% Save cluster information
trial_clusters = {};
for i=1:maxclusters
    if sum(clusters==i)>0 %in cases cluster was removed by above step
        trial_clusters(end+1)={find(clusters==i)};
    end
end