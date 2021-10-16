function test_cluster_params( allow_clustering_btwn_subjects, allow_cluster_merging, ...
    temp_dist_cutoffs, temp_dist_cutoffs_high, threshold_merging, ...
    subjects_all_original, subject_fake, distance_matrix, SampleNames_all, ...
    dirFigName, figFileName)


%% Purpose

% Looks at the number of clusters and fraction of samples clustered as a
% function of the dbscan distance cutoff...


%% Version History

% Arolyn, 2019.11.21: Original version


%% Do clustering at test distance cutoffs

% Low cutoffs (actually useful)
temp_dist_cutoffs = 1:1:100;
temp_num_clusters = temp_dist_cutoffs; % initialize same size
temp_frac_clustered = temp_dist_cutoffs; % initialize same size
for i=1:numel(temp_dist_cutoffs)
    % Clustering
    temp_dist_cutoff = temp_dist_cutoffs(i);
    clusters_all_temp = do_clustering ( ...
        allow_clustering_btwn_subjects, temp_dist_cutoff, allow_cluster_merging, threshold_merging, ...
        subjects_all_original, subject_fake, distance_matrix, SampleNames_all );
    % Record
    temp_num_clusters(i) = numel(clusters_all_temp);
    temp_frac_clustered(i) = sum( cellfun(@(x) numel(x), clusters_all_temp ) )/numel(SampleNames_all(subjects_all_original~=subject_fake));
end

% High cutoffs (reality check of limits)
temp_num_clusters_high = temp_dist_cutoffs_high; % initialize same size
temp_frac_clustered_high = temp_dist_cutoffs_high; % initialize same size
for i=1:numel(temp_dist_cutoffs_high)
    % Clustering
    temp_dist_cutoff = temp_dist_cutoffs_high(i);
    clusters_all_temp = do_clustering ( ...
        allow_clustering_btwn_subjects, temp_dist_cutoff, allow_cluster_merging, threshold_merging, ...
        subjects_all_original, subject_fake, distance_matrix, SampleNames_all );
    % Record
    temp_num_clusters_high(i) = numel(clusters_all_temp);
    temp_frac_clustered_high(i) = sum( cellfun(@(x) numel(x), clusters_all_temp ) )/numel(SampleNames_all(subjects_all_original~=subject_fake));
end


%% Make a plot

% Plot
figure(10)
clf(10)
subplot(2,2,1)
scatter( temp_dist_cutoffs,temp_num_clusters )
ylim([0 max(temp_num_clusters)+5])
xlabel('dbscan dist cutoff')
ylabel('number of clusters')
subplot(2,2,2)
scatter( temp_dist_cutoffs_high,temp_num_clusters_high )
ylim([0 max(temp_num_clusters)+5])
xlabel('dbscan dist cutoff')
ylabel('number of clusters')
subplot(2,2,3)
scatter( temp_dist_cutoffs,temp_frac_clustered )
ylim([0 1])
xlabel('dbscan dist cutoff')
ylabel('fraction clustered')
subplot(2,2,4)
scatter( temp_dist_cutoffs_high,temp_frac_clustered_high )
ylim([0 1])
xlabel('dbscan dist cutoff')
ylabel('fraction clustered')

%% Save plot

% Check directory existance
if ~exist(dirFigName,'dir')
    mkdir(dirFigName)
end
% Save
print([ dirFigName '/' figFileName ],'-dpng')


%%

end