function make_clustering_collectors_curves_downsample_specimens( ...
    allow_clustering_btwn_subjects, clustering_dist_cutoff, clustering_minpts, allow_cluster_merging, clustering_merging_threshold, ...
    subjects_final, specimen_number_final, distance_matrix_mini, clusters_final, unclustered_final, dir_collcurve )

% Parameters
test_num_trials = 100;
min_num_colonies_to_eval_subject = 20;

% list of subjects to downsample
my_subject_list = unique( subjects_final );
my_subject_list_num_samples = arrayfun(@(x) sum( x==subjects_final ), my_subject_list );
my_subject_list_keep = ( my_subject_list_num_samples >= min_num_colonies_to_eval_subject );
my_subject_list = my_subject_list( my_subject_list_keep );
my_subject_list_num_samples = my_subject_list_num_samples( my_subject_list_keep );

% For saving simulation results to compare later
my_subject_list_test_num_samples = {};
my_subject_list_mean_frac_unclustered = {};
my_subjects_list_mean_num_clusters = {}; % reclustered
my_subjects_list_mean_num_sets = {}; % just count if we have a colony from each cluster

% General plot appearance
fs = 20; % font size

% Loop through subjects
for s=1:numel(my_subject_list)

    % This subject
    my_subject = my_subject_list(s);
    my_samples = ( subjects_final==my_subject );
    % Get cluster for each sample
    clusters_final_plus = clusters_final; clusters_final_plus{end+1} = unclustered_final;
    my_samples_set = cell2mat( arrayfun(@(y) find(cellfun(@(x) ismember(y,x), clusters_final_plus )), find(my_samples), 'UniformOutput', false ) );
    my_samples_set( my_samples_set==numel(clusters_final)+1 ) = 0;
    % Specimens for this subject
    my_specimens = specimen_number_final( my_samples );
    my_specimens_list = unique( my_specimens );
    
    my_num_samples = sum(my_samples);
    my_num_specimens = numel(my_specimens_list);
    
    % Downsize distance matrix
    my_dist_matrix = distance_matrix_mini( my_samples, my_samples );

    % Things to test
    test_num_specimens = 1:1:my_num_specimens;
    record_frac_unclustered = zeros( my_num_specimens, test_num_trials );
    record_num_clusters = zeros( my_num_specimens, test_num_trials );
    record_num_sets = zeros( my_num_specimens, test_num_trials );

    for trial_num_specimens=1:my_num_specimens

        for j=1:test_num_trials

            trial_specimens = randperm( my_num_specimens, trial_num_specimens );
            trial_samples = find( ismember( my_specimens,my_specimens_list(trial_specimens)) );
            record_num_sets(trial_num_specimens,j) = numel( setdiff( unique( my_samples_set( trial_samples ) ), [0] ) );
            
            trial_distance_matrix = my_dist_matrix( trial_samples, trial_samples );

            trial_clusters = do_clustering_mini( trial_distance_matrix, clustering_dist_cutoff, clustering_minpts, allow_cluster_merging, clustering_merging_threshold );
            trial_num_clusters = numel( trial_clusters );
            trial_num_clustered = sum( cellfun(@(x) numel(x), trial_clusters) );
            trial_num_unclustered = numel(trial_samples) - trial_num_clustered;

            record_frac_unclustered(trial_num_specimens,j) = trial_num_unclustered/numel(trial_samples);
            record_num_clusters(trial_num_specimens,j) = trial_num_clusters;

        end

    end
    
    % Make a figure
    figure(20)
    clf(20)
    hold on
    box on
    scatter( test_num_specimens, 100*mean(record_frac_unclustered,2)')
    xlabel('number of samples (downsampled)')
    ylabel('percent of colonies that do not cluster')
    ylim([0 100])
    title(['subject ' char(my_subject+64) ' (n=' num2str(my_num_samples) ')'])
    set(gca,'FontSize',fs)
    hold off

    % Save figure
    print([ dir_collcurve '/CollectorCurve_DownSpec_PercentUnclustered_Subject-' char(my_subject+64)],'-dpng')

    
    % Make a figure
    figure(21)
    clf(21)
    hold on
    box on
    scatter( test_num_specimens, mean(record_num_clusters,2)')
    xlabel('number of samples (downsampled)')
    ylabel('percent of colonies that do not cluster')
    ylim([0 max(mean(record_num_clusters,2))+5])
    title(['subject ' char(my_subject+64) ' (n=' num2str(my_num_samples) ')'])
    set(gca,'FontSize',fs)
    hold off

    % Save figure
    print([ dir_collcurve '/CollectorCurve_DownSpec_NumClusters_Subject-' char(my_subject+64)],'-dpng')

    % Save data for comparison later
    my_subject_list_test_num_samples{end+1} = test_num_specimens;
    my_subject_list_mean_frac_unclustered{end+1} = 100*mean(record_frac_unclustered,2)';
    my_subjects_list_mean_num_clusters{end+1} = mean(record_num_clusters,2)';
    my_subjects_list_mean_num_sets{end+1} = mean(record_num_sets,2)';
    
end

%% Fraction unclustered: all subjects

my_colors = flipud(jet(numel(my_subject_list)));
[~,my_order] = sort(my_subject_list_num_samples, 'descend');
figure(22)
clf(22)
hold on
box on
for s=1:numel(my_subject_list)
    s_reorder = my_order(s);
    scatter( my_subject_list_test_num_samples{s_reorder}, my_subject_list_mean_frac_unclustered{s_reorder}, ...
        20, my_colors(s,:), 'MarkerEdgeAlpha', 0.5 )
end
lgd = legend( arrayfun(@(x) char(x+64), my_subject_list(my_order)), 'Location', 'northeastoutside' );
legend('boxoff')
title(lgd,'Subject')
xlabel('number of samples (downsampled)')
ylabel({'percent of colonies that are', 'not members of a lineage'})
ylim([0 100])
%title(['all subjects (num trials = ' num2str(test_num_trials) ')'])
set(gca,'FontSize',fs)
hold off

% Save figure
print([ dir_collcurve '/CollectorCurve_DownSpec_PercentUnclustered_Subject-all'],'-dpng')


%% Number of clusters: all subjects

my_colors = flipud(jet(numel(my_subject_list)));
[~,my_order] = sort(my_subject_list_num_samples, 'descend');
figure(23)
clf(23)
hold on
box on
for s=1:numel(my_subject_list)
    s_reorder = my_order(s);
    scatter( my_subject_list_test_num_samples{s_reorder}, my_subjects_list_mean_num_clusters{s_reorder}, ...
        20, my_colors(s,:), 'MarkerEdgeAlpha', 0.5)
end
lgd = legend( arrayfun(@(x) char(x+64), my_subject_list(my_order)), 'Location', 'northeastoutside' );
legend('boxoff')
title(lgd,'Subject')
xlabel('number of samples (downsampled)')
ylabel('number of lineages detected')
ylim([0 10])
%title(['all subjects (num trials = ' num2str(test_num_trials) ')'])
set(gca,'FontSize',fs)
hold off

% Save figure
print([ dir_collcurve '/CollectorCurve_DownSpec_NumClusters_Subject-all'],'-dpng')


%% Number of original clusters: all subjects

my_colors = flipud(jet(numel(my_subject_list)));
[~,my_order] = sort(my_subject_list_num_samples, 'descend');
figure(24)
clf(24)
hold on
box on
for s=1:numel(my_subject_list)
    s_reorder = my_order(s);
    scatter( my_subject_list_test_num_samples{s_reorder}, my_subjects_list_mean_num_sets{s_reorder}, ...
        20, my_colors(s,:), 'MarkerEdgeAlpha', 0.5)
end
lgd = legend( arrayfun(@(x) char(x+64), my_subject_list(my_order)), 'Location', 'northeastoutside' );
legend('boxoff')
title(lgd,'Subject')
xlabel('number of samples (downsampled)')
ylabel('number of original lineages hit')
ylim([0 10])
%title(['all subjects (num trials = ' num2str(test_num_trials) ')'])
set(gca,'FontSize',fs)
hold off

% Save figure
print([ dir_collcurve '/CollectorCurve_DownSpec_NumSets_Subject-all'],'-dpng')


end