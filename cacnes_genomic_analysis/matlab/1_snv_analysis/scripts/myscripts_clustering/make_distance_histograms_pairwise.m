function make_distance_histograms_pairwise( distance_matrix_mini, subjects_final, clustered_final, unclustered_final, SampleNamesLong_final, dir_save )


%% Directory setup

% Where to save figures
if ~exist( dir_save, 'dir' )
    mkdir( dir_save )
end


%% Pairwise distance

pairwise_distance_list =  distance_matrix_mini( ...
    1:1:length(distance_matrix_mini) > [ 1:1:length(distance_matrix_mini) ]' ... % only count each pair once
    & ~eye(length(distance_matrix_mini)) ... % avoid self-pairs
    & subjects_final == subjects_final' ... % require pairs to be from the same subject
    );

% Split off pairs that include an unclustered sample
is_unclustered = zeros( length(distance_matrix_mini),1, 'logical' );
is_unclustered( unclustered_final ) = 1;
pairwise_distance_list_unclustered =  distance_matrix_mini( ...
    1:1:length(distance_matrix_mini) > [ 1:1:length(distance_matrix_mini) ]' ... % only count each pair once
    & ~eye(length(distance_matrix_mini)) ... % avoid self-pairs
    & subjects_final == subjects_final' ... % require pairs to be from the same subject
    & ( is_unclustered | is_unclustered ) ...
    ); % pair contains at least one unclustered colony
pairwise_distance_list_clustered =  distance_matrix_mini( ...
    1:1:length(distance_matrix_mini) > [ 1:1:length(distance_matrix_mini) ]' ... % only count each pair once
    & ~eye(length(distance_matrix_mini)) ... % avoid self-pairs
    & subjects_final == subjects_final' ... % require pairs to be from the same subject
    & ~( is_unclustered | is_unclustered ) ...
    ); % pair contains two clustered colonies

% Split off pairs for same subject and distance subject
pairwise_distance_list_samesubj =  distance_matrix_mini( ...
    1:1:length(distance_matrix_mini) > [ 1:1:length(distance_matrix_mini) ]' ... % only count each pair once
    & ~eye(length(distance_matrix_mini)) ... % avoid self-pairs
    & subjects_final == subjects_final' ... % require pairs to be from the same subject
    ); % pair contains at least one unclustered colony
pairwise_distance_list_diffsubj =  distance_matrix_mini( ...
    1:1:length(distance_matrix_mini) > [ 1:1:length(distance_matrix_mini) ]' ... % only count each pair once
    & ~eye(length(distance_matrix_mini)) ... % avoid self-pairs
    & subjects_final ~= subjects_final' ... % require pairs to be from the different subjects
    ); % pair contains two clustered colonies
% Grab some numbers
sum(pairwise_distance_list_diffsubj>35)/numel(pairwise_distance_list_diffsubj) % 1
sum(pairwise_distance_list_samesubj<=35)/numel(pairwise_distance_list_samesubj) %  0.2788

%% Identify colony pairs that are suspiciously close

cluster_membership = cell2mat( cellfun(@(x) str2num(x(end-1:end)), SampleNamesLong_final, 'UniformOutput', false ) );

% Pairs from same subjects <100 SNVs apart
[ p1, p2 ] = find( ...
    1:1:length(distance_matrix_mini) > [ 1:1:length(distance_matrix_mini) ]' ... % only count each pair once
    & ~eye(length(distance_matrix_mini)) ... % avoid self-pairs
    & subjects_final == subjects_final' ... % require pairs to be from the different subjects
    & cluster_membership ~= cluster_membership' ... % not already in the same cluster
    & distance_matrix_mini <= 100 ... % closer than 100 SNVs
    ); % pair contains two clustered colonies
for i=1:numel(p1)
    fprintf(1, [ SampleNamesLong_final{p1(i)} ', ' SampleNamesLong_final{p2(i)} ', ' num2str(distance_matrix_mini(p1(i),p2(i))) '\n' ])
end
fprintf(1,['\n\n'])

% Pairs from different subjects <100 SNVs apart
[ p1, p2 ] = find( ...
    1:1:length(distance_matrix_mini) > [ 1:1:length(distance_matrix_mini) ]' ... % only count each pair once
    & ~eye(length(distance_matrix_mini)) ... % avoid self-pairs
    & subjects_final ~= subjects_final' ... % require pairs to be from the different subjects
    & distance_matrix_mini <= 100 ... % closer than 100 SNVs
    ); % pair contains two clustered colonies
for i=1:numel(p1)
    fprintf(1, [ SampleNamesLong_final{p1(i)} ', ' SampleNamesLong_final{p2(i)} ', ' num2str(distance_matrix_mini(p1(i),p2(i))) '\n' ])
end

% % Checking specific clusters
% thing = distance_matrix_mini( cluster_membership==37, cluster_membership==49 );
% mean( thing(:) )
% 
% thing = distance_matrix_mini( cluster_membership==33, cluster_membership==24 )
% mean( thing(:) )
% 
% thing = distance_matrix_mini( cluster_membership==21, cluster_membership==5 )
% mean( thing(:) )
% 
% thing = distance_matrix_mini( cluster_membership==16, cluster_membership==50 )
% mean( thing(:) )

thing = distance_matrix_mini( cluster_membership==38, cluster_membership==1 )
mean( thing(:) )


%% Binning instructions for distance histograms

% Bins for histograms (log10 distance)
bin_min = 0;
bin_max = 5;
bin_width = 0.25;
bin_edges = [ -Inf, bin_min:bin_width:bin_max ];
num_bins = numel(bin_edges)-1;

% Bin pairwise distances
counts = histcounts(log10(pairwise_distance_list),bin_edges);
counts_unclustered = histcounts(log10(pairwise_distance_list_unclustered),bin_edges);
counts_clustered = histcounts(log10(pairwise_distance_list_clustered),bin_edges);
counts_samesubj = histcounts(log10(pairwise_distance_list_samesubj),bin_edges);
counts_diffsubj = histcounts(log10(pairwise_distance_list_diffsubj),bin_edges);

% Print some bins
for i=1:numel(counts_samesubj)
    fprintf(1, [ num2str(counts_samesubj(i)) ', ' ])
end
fprintf(1, [ '\n' ])
for i=1:numel(counts_diffsubj)
    fprintf(1, [ num2str(counts_diffsubj(i)) ', ' ])
end
fprintf(1, [ '\n' ])

%% Make figure for ALL subjects combined

% Make a plot
figure(1)
clf(1)
% appearance
color_gray = 0.5*[ 1 1 1 ];
fs = 18;
%
% same cluster
hold on
box on
b=bar( counts, 'FaceColor', color_gray );
xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
%ylim( [ 0 10^5 ] )
set(gca, 'FontSize', fs )
ylabel('# pairs');
xlabel('distance (# SNVs)')
title('pairwise distance between colonies from the same subject')
hold off


%%

% Save
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [12 4]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 12 4]);
print([ dir_save '/' 'cluster-hists-pairwise_subj-ALL.png' ],'-dpng')


%% Make figure for ALL subjects combined: split unclustered

% Make a plot
figure(2)
clf(2)
% appearance
color_darkgray = 0.25*[ 1 1 1 ];
color_lightgray = 0.75*[ 1 1 1 ];
fs = 18;
%
% same cluster
hold on
box on
b=bar( [ counts_clustered; counts_unclustered ]', 'stacked', 'FaceColor', 'flat' );
b(1).FaceColor = color_darkgray;
b(2).FaceColor = color_lightgray;
xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
%ylim( [ 0 10^5 ] )
set(gca, 'FontSize', fs )
ylabel('# pairs');
xlabel('distance (# SNVs)')
title('pairwise distance between colonies from the same subject')
% Save
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [12 4]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 12 4]);
print([ dir_save '/' 'cluster-hists-pairwise_subj-ALL-split-unclus.png' ],'-dpng')
% Add legend
legend({'clustered+clustered','unclustered+any'}, 'Location', 'northeastoutside')
% end fig
hold off
% Save again
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [16 4]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 16 4]);
print([ dir_save '/' 'cluster-hists-pairwise_subj-ALL-split-unclus-wleg.png' ],'-dpng')


%% Make figure for ALL subjects combined: split same vs diff subject

% Make a plot
figure(3)
clf(3)
% appearance
color_darkgray = 0.25*[ 1 1 1 ];
color_lightgray = 0.75*[ 1 1 1 ];
fs = 18;
%
% same cluster
hold on
box on
b=bar( [ counts_samesubj; counts_diffsubj ]', 'stacked', 'FaceColor', 'flat' );
b(1).FaceColor = color_darkgray;
b(2).FaceColor = color_lightgray;
xlim( [ 1-0.5 bin_max/bin_width+2+0.5 ] )
xticks( [1.5 2:1/bin_width:bin_max/bin_width+2 ]-0.5 )
xticklabels( [0 arrayfun(@(x) 10^x, bin_min:1:bin_max) ] )
%ylim( [ 0 5*10^5 ] )
set(gca, 'FontSize', fs )
ylabel('# pairs');
xlabel('distance (# SNVs)')
title('pairwise distance between colonies')
% Save
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 4]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 4]);
print([ dir_save '/' 'cluster-hists-pairwise_subj-ALL-split-same-diff.png' ],'-dpng')
% Add legend
legend({'same subject','different subjects'}, 'Location', 'northwest')
% end fig
hold off
% Save again
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 4]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
print([ dir_save '/' 'cluster-hists-pairwise_subj-ALL-split-same-diff-wleg.png' ],'-dpng')


end