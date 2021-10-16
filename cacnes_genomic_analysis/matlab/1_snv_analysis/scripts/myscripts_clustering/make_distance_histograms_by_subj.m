function make_distance_histograms_by_subj( clusters_all, unclustered_all, distance_matrix, clusters_all_subjects, subjects_all, subject_fake, dir_clustering )


%% Directory setup

% Where to save figures
dir_save = 'SummaryHistograms';
if ~exist( [ dir_clustering '/' dir_save ], 'dir' )
    mkdir( [ dir_clustering '/' dir_save ] )
end


%% Binning instructions for distance histograms

% Bins for histograms (log10 distance)
bin_min = 0;
bin_max = 5;
bin_width = 0.25;
bin_edges = [ -Inf, bin_min:bin_width:bin_max ];
num_bins = numel(bin_edges)-1;


%% Get subject info for clusters

% Cluster membership for each sample
cluster_membership_all = zeros( size( subjects_all ) );
for c=1:numel(clusters_all)
    cluster_membership_all( clusters_all{c} ) = c;
end

% List all subjects
subjects_list = setdiff( unique(subjects_all), subject_fake );


%% Initialize counters for across all subjects

% Initialize counters for binning distances over all subjects
counts_all_sameclus = zeros( 1,num_bins ); % initialize
counts_all_samesubj = zeros( 1,num_bins ); % initialize
counts_all_diffsubj = zeros( 1,num_bins ); % initialize
counts_all_samesubj_unclus = zeros( 1,num_bins ); % initialize
counts_all_diffsubj_unclus = zeros( 1,num_bins ); % initialize


%% Loop through all subjects

for s=1:numel(subjects_list) 

    % Subject info
    next_subject = subjects_list(s);
    next_subject_char = char( next_subject+64 );
    next_subject_clusters = find( cellfun(@(x) ismember(next_subject_char,x), clusters_all_subjects ) );
%     if numel(next_subject_clusters) <= 1
%         continue
%     end

    % Distance matrix for this subject
    %dm = distance_matrix( subjects_all==next_subject, subjects_all==next_subject );

    % Minimum distance from each colony to another colony in the same cluster
    within_cluster_min_dists = cell( numel(next_subject_clusters),1 );
    for k=1:numel(next_subject_clusters)
        this_cluster_indices = clusters_all{next_subject_clusters(k)};
        within_cluster_min_dists{k} = arrayfun(@(x) ...
            min( distance_matrix( x, setdiff(this_cluster_indices,x) ) ), ...
            this_cluster_indices ); 
    end

    % Mininum distance from each colony to another colony in a different cluster
    btwn_cluster_min_dists = cell( numel(next_subject_clusters),1 );
    if numel(next_subject_clusters) > 1
        for k=1:numel(next_subject_clusters)
            this_cluster_indices = clusters_all{next_subject_clusters(k)};
            other_cluster_indices = setdiff( find( subjects_all==next_subject & cluster_membership_all>0 ), this_cluster_indices );
            btwn_cluster_min_dists{k} = arrayfun(@(x) ...
                min( distance_matrix( x, other_cluster_indices ) ), ...
                this_cluster_indices ); 
        end
    else
        btwn_cluster_min_dists{1} = [];
    end
    
    % Minimum distance from each colony to another colony from a different subject
    btwn_subj_min_dists = cell( numel(next_subject_clusters),1 );
    for k=1:numel(next_subject_clusters)
        this_cluster_indices = clusters_all{next_subject_clusters(k)};
        other_subj_indices = setdiff( find( subjects_all~=next_subject & cluster_membership_all>0 ), this_cluster_indices );
        btwn_subj_min_dists{k} = arrayfun(@(x) ...
            min( distance_matrix( x, other_subj_indices ) ), ...
            this_cluster_indices ); 
    end
    
    % Minimum distance from each unclustered colony to another colony in any cluster (same subj)
    unclustered_indices = intersect( unclustered_all, find(subjects_all==next_subject) );
    clustered_indices = find( subjects_all==next_subject & cluster_membership_all>0 );
    unclustered_to_clustered_samesubj_min_dists = arrayfun(@(x) ...
        min( distance_matrix( x, clustered_indices ) ), ...
        unclustered_indices );

    % Minimum distance from each unclustered colony to another colony in any cluster (same subj)
    unclustered_indices = intersect( unclustered_all, find(subjects_all==next_subject) );
    clustered_indices = find( subjects_all~=next_subject & cluster_membership_all>0 );
    unclustered_to_clustered_diffsubj_min_dists = arrayfun(@(x) ...
        min( distance_matrix( x, clustered_indices ) ), ...
        unclustered_indices );

    
    %% Plotting

    % Make a plot
    figure(1)
    clf(1)
    % appearance
    color_gray = 0.75*[ 1 1 1 ];
    if numel(next_subject_clusters) > 1
        colors_clusters = parula( numel(next_subject_clusters) );
    else
        colors_clusters = parula( 5 );
        colors_clusters = colors_clusters( 3,: );
    end
    fs = 12;
    legend_labels = arrayfun(@(x) [ 'cluster ' num2str(x) ], next_subject_clusters, 'UniformOutput', false );
    %
    % same cluster
    subplot(5,1,1)
    hold on
    box on
    counts = cell2mat( cellfun(@(x) histcounts(log10(x),bin_edges), within_cluster_min_dists, 'UniformOutput', false) );
    if numel(next_subject_clusters) > 1
        counts_all_sameclus = counts_all_sameclus + sum(counts);
    else
        counts_all_sameclus = counts_all_sameclus + counts;
    end
    b=bar( counts', 'stacked', 'FaceColor', 'flat' );
    for k=1:numel(next_subject_clusters)
        b(k).FaceColor = colors_clusters( k,: );
    end
    xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
    xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
    xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
    ylim( [ 0 1.1*max(sum(counts)) ] )
    set(gca, 'FontSize', fs )
    title('distance to the closest colony in the same cluster')
    ylabel('num of colonies')
    l=legend( legend_labels, 'Location', 'northeastoutside' );
    l.Title.String = [ 'subject ' next_subject_char ' clusters' ];
    l.Title.FontSize = fs;
    hold off
    %
    % diff cluster, same subj
    subplot(5,1,2)
    hold on
    box on
    counts = cell2mat( cellfun(@(x) histcounts(log10(x),bin_edges), btwn_cluster_min_dists, 'UniformOutput', false) );
    if numel(next_subject_clusters) > 1
        counts_all_samesubj = counts_all_samesubj + sum(counts);
    else
        counts_all_samesubj = counts_all_samesubj + counts;
    end
    b=bar( counts', 'stacked', 'FaceColor', 'flat' );
    for k=1:numel(next_subject_clusters)
        b(k).FaceColor = colors_clusters( k,: );
    end
    xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
    xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
    xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
    ylim( [ 0 max(1,1.1*max(sum(counts))) ] )
    set(gca, 'FontSize', fs )
    title('distance to the closest colony in a different cluster (same subject)')
    ylabel('num of colonies')
    l=legend( legend_labels, 'Location', 'northeastoutside' );
    l.Title.String = [ 'subject ' next_subject_char ' clusters' ];
    l.Title.FontSize = fs;
    hold off
    %
    % diff cluster, diff subj
    subplot(5,1,3)
    hold on
    box on
    counts = cell2mat( cellfun(@(x) histcounts(log10(x),bin_edges), btwn_subj_min_dists, 'UniformOutput', false) );
    if numel(next_subject_clusters) > 1
        counts_all_diffsubj = counts_all_diffsubj + sum(counts);
    else
        counts_all_diffsubj = counts_all_diffsubj + counts;
    end
    b=bar( counts', 'stacked', 'FaceColor', 'flat' );
    for k=1:numel(next_subject_clusters)
        b(k).FaceColor = colors_clusters( k,: );
    end
    xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
    xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
    xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
    ylim( [ 0 1.1*max(sum(counts)) ] )
    set(gca, 'FontSize', fs )
    title('distance to closest colony in a different cluster (different subject)')
    ylabel('num of colonies')
    l=legend( legend_labels, 'Location', 'northeastoutside' );
    l.Title.String = [ 'subject ' next_subject_char ' clusters' ];
    l.Title.FontSize = fs;
    hold off
    %
    % unclustered (same subj)
    subplot(5,1,4)
    hold on
    box on
    counts = histcounts(log10(unclustered_to_clustered_samesubj_min_dists),bin_edges);
    counts_all_samesubj_unclus = counts_all_samesubj_unclus + counts;
    bar( counts, 'FaceColor', color_gray )
    xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
    xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
    xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
    ylim( [ 0 max(1,1.1*max(counts)) ] )
    yticks( 0:1:max(1,1.1*max(counts)) )
    set(gca, 'FontSize', fs )
    title('distance to closest colony in any cluster (same subject)')
    ylabel('num of colonies')
    l=legend( 'unclustered', 'Location', 'northeastoutside' );
    l.Title.String = [ 'subject ' next_subject_char ' clusters' ];
    l.Title.FontSize = fs;
    hold off
    %
    % unclustered (diff subj)
    subplot(5,1,5)
    hold on
    box on
    counts = histcounts(log10(unclustered_to_clustered_diffsubj_min_dists),bin_edges);
    counts_all_diffsubj_unclus = counts_all_diffsubj_unclus + counts;
    bar( counts, 'FaceColor', color_gray )
    xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
    xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
    xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
    ylim( [ 0 max(1,1.1*max(counts)) ] )
    yticks( 0:1:max(1,1.1*max(counts)) )
    set(gca, 'FontSize', fs )
    title('distance to closest colony in any cluster (different subject)')
    xlabel('distance (~#SNVs, log10)')
    ylabel('num of colonies')
    l=legend( 'unclustered', 'Location', 'northeastoutside' );
    l.Title.String = [ 'subject ' next_subject_char ' clusters' ];
    l.Title.FontSize = fs;
    hold off
    %%

    % Save
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [10 12]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 10 12]);
    print([ dir_clustering '/' dir_save '/' 'cluster-hists-min_subj-' next_subject_char '.png' ],'-dpng')

end


%% Make figure for ALL subjects combined

% Make a plot
figure(2)
clf(2)
% appearance
color_gray = 0.75*[ 1 1 1 ];
% colors_clusters = parula( 5 );
% color_cluster = colors_clusters( 3,: );
color_cluster = 0.25*[ 1 1 1 ];
fs = 16;
yl_pos = -0.0625;
ymax_1 = 420;
ymax_2 = 27.5;
%
% same cluster
subplot(5,1,1)
hold on
box on
b=bar( counts_all_sameclus', 'FaceColor', color_cluster );
xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
ylim( [ 0 ymax_1 ] )
yticks( 0:100:ymax_1 )
set(gca, 'FontSize', fs )
title('minimum distance to closest colony in the same cluster')
yl=ylabel('# colonies');
set(yl, 'Units', 'Normalized', 'Position', [yl_pos, 0.5, 0]);
hold off
%
% diff cluster, same subj
subplot(5,1,2)
hold on
box on
b=bar( counts_all_samesubj', 'FaceColor', color_cluster );
xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
ylim( [ 0 ymax_1 ] )
yticks( 0:100:ymax_1 )
set(gca, 'FontSize', fs )
title('minimum distance to closest colony in another cluster (same subject)')
yl=ylabel('# colonies');
set(yl, 'Units', 'Normalized', 'Position', [yl_pos, 0.5, 0]);
hold off
%
% diff cluster, diff subj
subplot(5,1,3)
hold on
box on
b=bar( counts_all_diffsubj', 'FaceColor', color_cluster );
xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
ylim( [ 0 ymax_1 ] )
yticks( 0:100:ymax_1 )
set(gca, 'FontSize', fs )
title('minimum distance to closest colony in another cluster (different subject)')
yl=ylabel('# colonies');
set(yl, 'Units', 'Normalized', 'Position', [yl_pos, 0.5, 0]);
hold off
%
% unclustered (same subj)
subplot(5,1,4)
hold on
box on
bar( counts_all_samesubj_unclus, 'FaceColor', color_gray )
xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
ylim( [ 0 ymax_2 ] )
yticks( 0:10:ymax_2 )
set(gca, 'FontSize', fs )
title('minimum distance to closest colony in any cluster (same subject)')
yl=ylabel('# colonies');
set(yl, 'Units', 'Normalized', 'Position', [yl_pos, 0.5, 0]);
hold off
%
% unclustered (diff subj)
subplot(5,1,5)
hold on
box on
bar( counts_all_diffsubj_unclus, 'FaceColor', color_gray )
xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
ylim( [ 0 ymax_2 ] )
yticks( 0:10:ymax_2 )
set(gca, 'FontSize', fs )
title('minimum distance to closest colony in any cluster (different subject)')
xlabel('distance (# SNVs)')
yl=ylabel('# colonies');
set(yl, 'Units', 'Normalized', 'Position', [yl_pos, 0.5, 0]);
hold off
%%

% Save
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [12 12]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 12 12]);
print([ dir_clustering '/' dir_save '/' 'cluster-hists-min_subj-ALL.png' ],'-dpng')

%% Legend

figure(3)
clf(3)
% params
sp = 0.1;
sp_text = 5000;
fs_leg = 100;
dx_leg = 33000;
xlim_leg = 200000;
% figure
hold on
axis off
r1=rectangle( 'Position', [ 0+sp 1+sp dx_leg-2*sp 1-2*sp ], 'FaceColor', color_cluster ); 
t1=text( dx_leg+sp_text, 1.5, 'clustered colonies', 'FontSize', fs_leg, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );
r2=rectangle( 'Position', [ 0+sp 0+sp dx_leg-2*sp 1-2*sp ], 'FaceColor', color_gray ); 
t2=text( dx_leg+sp_text, 0.5, 'unclustered colonies', 'FontSize', fs_leg, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );
xlim([0 xlim_leg])
hold off

% Save
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [20 4]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 20 4]);
print([ dir_clustering '/' dir_save '/' 'cluster-hists-min_subj-ALL_legend.png' ],'-dpng')

end