function [ clusters_all_sorted, clusters_all_subjects ] = sort_cluster_order( clusters_all, SampleNames_all )

% This script reorders clusters from biggest to smallest and returns the
% subject(s) that are part of each cluster

% Input:
% % clusters_all
% % SampleNames_all

% Output: 
% % clusters_all_sorted: clusters in order from biggest to smallest
% % clusters_all_subjects: string subject tags for all subjects included in
% each cluster (same indexing as sorted clusters)


%% Version History
% Arolyn, 2018.11.13: Initial version


%% Get sizes of all clusters

% Size of clusters
cluster_sizes = cellfun(@(x) numel(x), clusters_all);
% Sort by size
[cluster_size_order, cluster_index_order] = sort(cluster_sizes,'descend');

% Reorder clusters with largest one first
clusters_all_sorted = {}; % new cell for clusters
for c=1:numel(cluster_index_order);
    clusters_all_sorted{end+1} = clusters_all{cluster_index_order(c)};
end


%% Now get a representative subject from each cluster

clusters_all_subjects = {}; % cell for representative subjects
for c=1:numel(clusters_all_sorted)
    this_cluster = clusters_all_sorted{c};
    this_cluster_samplenames = SampleNames_all(this_cluster);
    this_cluster_subjects = cellfun(@(x) x(1), this_cluster_samplenames);
    %this_subject = mode(this_cluster_subjects); % most common subject in this cluster
    %clusters_all_subjects{end+1} = this_subject;
    these_subjects = unique(this_cluster_subjects);
    clusters_all_subjects{end+1} = these_subjects;
end

