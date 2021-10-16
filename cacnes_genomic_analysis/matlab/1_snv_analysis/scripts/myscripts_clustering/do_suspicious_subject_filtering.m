function [ subjects_all_filtered_all, removed_list ] = do_suspicious_subject_filtering( ...
    min_num_colonies_from_minority_subject, clusters_all, ...
    subjects_all_filtered_2, subject_fake, ...
    SampleNames_all, ...
    dir_clustering )

% Rename original clusters and subjects objects
clusters_all_original = clusters_all;
subjects_all = subjects_all_filtered_2;

% Log file (open)
fid_log = fopen( [dir_clustering '/Log_ClusterFilter_SuspiciousSubjects.txt'], 'wt' ); 
fprintf(fid_log,[ 'Suspicious Sample Removal' '\n' 'Starting log file...' '\n']);
fprintf(fid_log, datestr(now, 'yyyy-mm-dd HH:MM:SS') );
fprintf(fid_log,'\n \n');

% Track samples removed in this step for filter tracking
removed_list = [];

fprintf(1,['First checking for suspicious subjects...' '\n']);
fprintf(fid_log,['\n \nFirst checking for suspicious subjects...' '\n']);
for k=1:numel(clusters_all_original) % Loop through a bunch of clusters
    % Identify samples in this cluster
    clustersamples = clusters_all{k}; % in cluster and in outrgoup
    fprintf(1,['Cluster-' num2str(k) ':\n']);
    fprintf(fid_log,['\nCluster-' num2str(k) ':\n']);
    fprintf(1,['Number of samples: ' num2str(sum(clustersamples>0)) '\n']);
    fprintf(fid_log,['Number of samples: ' num2str(sum(clustersamples>0)) '\n']);
    clusterSampleNames=SampleNames_all(clustersamples);
    subjects=subjects_all(clustersamples);
    remove = zeros( numel(clustersamples),1 ); % boolean to keep track of samples to remove
    % Get subjects in this clade
    clade_subjects = unique(subjects);
    clade_subject_counts = arrayfun(@(x) sum(x==subjects), clade_subjects);
    % Find the identify of the most common subject
    [max_counts,max_index] = max(clade_subject_counts);
    fprintf(1,[ 'Majority subject: Subject ' char(64+clade_subjects(max_index)) ', ' ...
        num2str(max_counts) '/' num2str(sum(clade_subject_counts)) ' samples in clade' '\n' ])
    fprintf(fid_log,[ 'Majority subject: Subject ' char(64+clade_subjects(max_index)) ', ' ...
        num2str(max_counts) '/' num2str(sum(clade_subject_counts)) ' samples in clade' '\n' ]);
    if numel(clade_subjects)>1 % if clade has multiple subjects
        for i=1:numel(clade_subjects)
            if clade_subjects(i) ~= clade_subjects(max_index) % if not majority subject
                minority_counts_clade = clade_subject_counts(i);
                minority_counts_all = sum( subjects_all == clade_subjects(i) );
                fprintf(1,[ 'Minority subject: Subject ' char(64+clade_subjects(i)) ', ' ...
                    num2str(minority_counts_clade) '/' num2str(sum(clade_subject_counts)) ' samples in clade, ' ...
                    num2str(minority_counts_clade) '/' num2str(minority_counts_all) ' samples for subject' '\n' ])
                fprintf(fid_log,[ 'Minority subject: Subject ' char(64+clade_subjects(i)) ', ' ...
                    num2str(minority_counts_clade) '/' num2str(sum(clade_subject_counts)) ' samples in clade, ' ...
                    num2str(minority_counts_clade) '/' num2str(minority_counts_all) ' samples for subject' '\n' ]);
                % Print out a warning if some of the subjects are suspicious...
                if minority_counts_clade < min_num_colonies_from_minority_subject
                % Remove all minority subjects from clades % Aro_Change 2020.01.04
                    fprintf(1,'Warning! Suspicious subject assignment! Removing sample(s). \n')
                    fprintf(fid_log,'Warning! Suspicious subject assignment! Removing sample(s). \n');
                    % Flag suspicious samples for removal
                    to_remove = find( clade_subjects(i)==subjects );
                    remove(to_remove) = 1;
                    for s=1:numel(to_remove)
                        fprintf(1,['     ' SampleNames_all{clustersamples(to_remove(s))} '\n'])
                        fprintf(fid_log,['     ' SampleNames_all{clustersamples(to_remove(s))} '\n']);
                    end
                end % if minority subject is suspicious
            end % if not the majority subject
        end % loop through subjects in clade
    end % if clade has multiple subjects
    % Remove suspicous samples here
    clusters_all{k} = clustersamples(~remove);
    % Record samples removed in filter tracking
    if sum(remove)>0
        removed_list = [ removed_list, clustersamples(remove>0) ];
    end
end % loop through clusters

% Sort clusters by size and find subject
[ clusters_all, clusters_all_subjects ] = ...
    sort_cluster_order( clusters_all, SampleNames_all );

% Log file (close)
fclose( fid_log ); 

% Update subjects_all
subjects_all_filtered_all = subjects_all_filtered_2;
subjects_all_filtered_all( removed_list ) = subject_fake;


end