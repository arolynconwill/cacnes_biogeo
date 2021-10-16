function clusters_all = do_clustering ( allow_clustering_btwn_subjects, dist_cutoff, clustering_minpts, allow_cluster_merging, threshold_merging, subjects_all, subject_fake, distance_matrix, SampleNames_all )

% Inputs
% % allow_clustering_btwn_subjects: whether or not to allow clustering
% between subjects; 0 = only cluster within each subject; 1 = allow samples
% from different subjects to cluster together
% % dist_cutoff: for dbscan
% % allow_cluster_merging: whether or not to allow clusters to merge
% together
% % threshold_merging: threshold for merging
% % subjects_all: subjects corresponding to each sample
% % subject_fake: subject tag for a bad sample (not a real subject)
% % distance_matrix: distance matrix
% % SampleNames_all: names of all samples; indexed like distance_matrix and
% subjects_all


% Outputs
% % clusters_all: clusters generated here; cell with lists of sample
% indices


%% Version history
% 2018.11.13 Arolyn: Initial version
% 2018.11.14 Arolyn: Allow for multiple types of clustering


%% Initialize variables

% Initialize variable for storing clusters
clusters_all={};


%% Case 1: Cluster each subject separately

if ~allow_clustering_btwn_subjects

    % Find max real subject number
    if max(subjects_all) ~= subject_fake % no s=100's
        max_subject_index = max(subjects_all); % take max subject number
    else
        all_subject_indices = sort(unique(subjects_all)); % s=100 exists
        max_subject_index = all_subject_indices(end-1); % take second to max subject number
    end

    % Cluster each subject separately
    for s=1:max_subject_index % Loop through all subjects
        indices_this_subject=find(subjects_all==s);
        n=SampleNames_all(indices_this_subject);

        % Only perform clustering if there are at least four samples for this
        % subject
        if sum(subjects_all==s) >= clustering_minpts % >4

            %find clusters using a strict cutoff
            dm=distance_matrix(subjects_all==s, subjects_all==s);
            dm_i=max(dm(:))-dm;
                % distance matrix is max of distance matrix minus distance
                % matrix, so actually similarity
            clusters=dbscan(dm_i,max(dm_i(:))-dist_cutoff,clustering_minpts);
                % "Neighborhood" is max similarity - cutoff, so things need to  
                % be more similar than this to be clustered together;
                % i.e. a smaller cutoff means a higher similarity threshold!
                % Aro_Change: max(dm(:))-dist_cutoff changed to max(dm_i(:))-dist_cutoff
                % (This actually gives the same value but easier is to read if 
                % it is relative to max similarity rather than dissimilarity.)
            maxclusters=max(clusters);
            biggestcluster=mode(clusters(clusters>0));
            numbiggestcluster=sum(clusters==biggestcluster);

            % Merge clusters from this subject that are close
            if allow_cluster_merging % only merge if merging is on
                for i=1:maxclusters          
                    for j=(i+1):maxclusters
                        distance_ij=dm(clusters==i,clusters==j);
                        if mean(distance_ij(:)) < threshold_merging
                            fprintf(1,['Merging clusters ' num2str(i) ' and ' num2str(j) '\n'])
                            disp([n(find(clusters==i,1)) n(find(clusters==j,1))])
                            disp([min(distance_ij(:)) max(distance_ij(:)) median(distance_ij(:)) mean(distance_ij(:))])
                            clusters(clusters==j)=i; 
                        end
                    end
                end
            end

            % Save cluster information
            for i=1:maxclusters
                if sum(clusters==i)>0 %in cases cluster was removed by above step
                    clusters_all(end+1)={indices_this_subject(clusters==i)};
                end
            end

        end % end for if statement on having enough samples for this subject
    end % end of for loop through subjects

end % end of if statement


%% Case 2: Cluster all samples together

if allow_clustering_btwn_subjects

    % Cluster EVERYTHING
    indices_this_subject = find( subjects_all < subject_fake ); % Cluster ALL samples that passed filters, excluding s=100
    n = SampleNames_all(indices_this_subject);

    %find clusters using a strict cutoff
    dm=distance_matrix(indices_this_subject,indices_this_subject);
    dm_i=max(dm(:))-dm;
        % distance matrix is max of distance matrix minus distance
        % matrix, so actually similarity
    clusters=dbscan(dm_i,max(dm_i(:))-dist_cutoff,3);
        % "Neighborhood" is max similarity - cutoff, so things need to  
        % be more similar than this to be clustered together;
        % i.e. a smaller cutoff means a higher similarity threshold!
        % Aro_Change: max(dm(:))-dist_cutoff changed to max(dm_i(:))-dist_cutoff
        % (This actually gives the same value but easier is to read if 
        % it is relative to max similarity rather than dissimilarity.)
    maxclusters=max(clusters);
    biggestcluster=mode(clusters(clusters>0));
    numbiggestcluster=sum(clusters==biggestcluster);

    % Merge clusters that are close
    if allow_cluster_merging % only merge if merging is on
        for i=1:maxclusters          
            for j=(i+1):maxclusters
                distance_ij=dm(clusters==i,clusters==j);
                %fprintf(1,['Distance ' num2str(i) '-' num2str(j) ': ' num2str(mean(distance_ij(:))) '\n'])
                if mean(distance_ij(:)) < threshold_merging
                    fprintf(1,['Merging clusters ' num2str(i) ' and ' num2str(j) '\n'])
                    %disp([n(clusters==i)]); disp([ n(clusters==j)]);
                    disp([n(find(clusters==i,1)) n(find(clusters==j,1))])
                    disp([min(distance_ij(:)) max(distance_ij(:)) median(distance_ij(:)) mean(distance_ij(:))])
                    clusters(clusters==j)=i; 
                end
            end
        end
    end

    % Save cluster information
    for i=1:maxclusters
        if sum(clusters==i)>0 %in cases cluster was removed by above step
            clusters_all(end+1)={indices_this_subject(clusters==i)};
        end
    end
    
end