function subjects_all = do_within_cluster_examination_woutgroup( ...
    subjects_all, clusters_all_subject, clusters_all, outgroup, SampleNames_all, ...
    Calls_all, coverage_all, maf_all, distance_matrix, p_all, ...
    figures_boolean, pause_boolean, subject_fake, ...
    min_size_to_examine_cluster, min_pos_to_examine_cluster, ...
    dirFigName)

%% Summary

% This function makes diagnostic plots of final clusters to make sure
% filtering went well.

% Inputs:
% % subjects_all: Subjects for each sample
% % clusters_all_subject: Subject tag for each initial cluster
% % clusters_all: Final clusters
% % outgroups: Final outgroups
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
% % dirFigName: Name of directory to put figures in

% Outputs:
% % subjects_all: Subjects for each sample, where subject ID for bad
% samples has been changed to subject_fake



%% Check filtering

% Examine each cluster
for this_cluster = 1:numel(clusters_all) % loop through clusters...
    this_subject = [clusters_all_subject{this_cluster} '-' num2str(numel(clusters_all{this_cluster}))]; 
    if numel(clusters_all{this_cluster}) >= min_size_to_examine_cluster % examine clusters that are big enough
        mySetPositionsImportantNum = ...
            examine_samples_within_set_woutgroup( this_cluster, clusters_all, outgroup, this_subject, ...
            SampleNames_all, Calls_all, coverage_all, maf_all, distance_matrix, p_all, ...
            figures_boolean, min_pos_to_examine_cluster, ...
            dirFigName);
        % Pause after analyzing each cluster, if requested by user
        if pause_boolean
            pause % pause after each filtering step to inspect figures % comment this out if you do not want to inspect each cluster
        end
    else % otherwise skip examining this cluster
        fprintf(1,['(Filtering) Skipping this cluster (too small): ' 'Subject-' this_subject '_Cluster-' num2str(this_cluster) '\n'])
    end
end


end
