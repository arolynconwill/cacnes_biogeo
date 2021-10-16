function subjects_all = do_within_cluster_filtering ( subjects_all, ...
    clusters_all_subject, clusters_all_initial_sorted, clusters_all_for_filtering, ...
    SampleNames_all, ...
    Calls_all, coverage_all, maf_all, Quals_all, ...
    distance_matrix, p_all, ...
    figures_boolean, pause_boolean, subject_fake, ...
    min_size_to_examine_cluster, min_pos_to_examine_cluster, ...
    filtering_mean_maf_per_sample, filtering_min_cov_to_consider_pos, ...
    dirFigName)

%% Summary

% This function is a wrapper function for within_cluster filtering.

% Inputs:
% % subjects_all: Subjects for each sample
% % clusters_all_subject: Subject tag for each initial cluster
% % clusters_all_initial_sorted: Initial clusters
% % SampleNames_all: Names of samples
% % Calls_all: Used in filtering
% % coverage_all: Used in filtering
% % maf_all: Used in filtering
% % distance_matrix: Used in filtering
% % p_all: Used in filtering
% % figures_boolean: Whether or not to save figures
% % pause_boolean: Whether or not to pause at each cluster for manual
% inspection
% % subject_fake: Numeric tag to indicate a bad sample
% % min_size_to_examine_cluster: Only filter cluster if above this size
% % min_pos_to_examine_cluster: Only filter if enough positions identified
% % filtering_mean_maf_per_sample; % minimum allowed mean major allele frequency over candidate SNPs for a sample to pass filtering  % passed through
% % filtering_min_cov_to_consider_pos; % minimum coverage for a position to be included in sample filtering  % passed through
% % dirFigName: Name of directory to put figures in

% Outputs:
% % subjects_all: Subjects for each sample, where subject ID for bad
% samples has been changed to subject_fake


%% Filtering

% Filter each cluster
for this_cluster = 1:numel(clusters_all_initial_sorted) % loop through clusters...
    this_subject = [clusters_all_subject{this_cluster} '-' num2str(numel(clusters_all_initial_sorted{this_cluster}))]; % Aro_Change: 0 to indicate that initial clustering was not done by subject
    %this_subject = SampleNames_all{clusters_all_original{this_cluster}(1)}(1); % grab letter indicating subject name off of first sample name
    if clusters_all_for_filtering(this_cluster) && ...
            (numel(clusters_all_initial_sorted{this_cluster}) >= min_size_to_examine_cluster) % examine clusters that are big enough
        % Aro_Note: The following module filters samples within each
        % cluster. VERY IMPORTANT: There are filtering parameters inside
        % these scripts, which you may want to change depending on your
        % sample set!!!
        [ mySubjectSamplesToRemoveIndexedAll, mySetPositionsImportantNum ] = ...
            examine_samples_within_set( this_cluster, clusters_all_initial_sorted, this_subject, ...
            SampleNames_all, Calls_all, coverage_all, maf_all, Quals_all, distance_matrix, p_all, ...
            figures_boolean, min_pos_to_examine_cluster, filtering_mean_maf_per_sample, filtering_min_cov_to_consider_pos, ...
            dirFigName);
        % Pause after analyzing each cluster, if requested by user
        if pause_boolean
            pause % pause after each filtering step to inspect figures % comment this out if you do not want to inspect each cluster
        end
        % Tag bad samples by changing their subject number, but only if
        % sufficient positions were considered when examining samples
        % within the set
        if numel(mySubjectSamplesToRemoveIndexedAll) > 0 && mySetPositionsImportantNum >= min_pos_to_examine_cluster % only remove if samples to remove
            subjects_all(mySubjectSamplesToRemoveIndexedAll) = subject_fake;
        end
    else % otherwise skip examining this cluster
        if ~clusters_all_for_filtering(this_cluster)
            fprintf(1,['(Filtering) Skipping this cluster (no change from previous round): Cluster-' num2str(this_cluster) '-' this_subject '\n'])
        else
            fprintf(1,['(Filtering) Skipping this cluster (too small): Cluster-' num2str(this_cluster) '-' this_subject '\n'])
        end
    end
end

% Print total number of samples across all clusters that were filtered in this step
fprintf(1,['(Filtering) Total number of samples filtered: ' num2str(sum(subjects_all == subject_fake)) '\n'])


end