%% FIGURE 2 SCRIPT


%% DATA

load('2_snvs/cluster_names')
clusters_all = load('data/cluster_step_variables','clusters_final'); clusters_all = clusters_all.clusters_final;
clusters_all_slst_super = cellfun(@(x) x(1), clusters_all_slst);
SampleNamesLong_all = load('data/cluster_step_variables','SampleNamesLong_final'); SampleNamesLong_all = SampleNamesLong_all.SampleNamesLong_final;
sample_clusters = cellfun(@(x) str2num(x(end-1:end)), SampleNamesLong_all );
subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;


%% LINEAGE/SLST ABUNDANCES BY SUBJECT 

% 2A: Shows relative abundances of lineages colored by SLST for each
% subject

% Implementation notes:
% % Unclustered colonies are gray
% % For now, treating all colonies as independent (even though later we will learn they are not)

% Generate data table: number of colonies in each lineage by subject
my_subjects_list = unique(subjects_all);
my_lineages_list = [ 1:1:numel(clusters_all), 0]; % always includes unclustered and puts it at the end
% Generate table: number of colonies in each zone from each lineage
my_table_subject_lineages = zeros( numel(my_subjects_list), numel(my_lineages_list) );
my_table_subject_lineages_pores = zeros( numel(my_subjects_list), numel(my_lineages_list) );
my_table_subject_lineages_surface = zeros( numel(my_subjects_list), numel(my_lineages_list) );
for s=1:numel(my_subjects_list)
    for c=1:numel(my_lineages_list)
        my_table_subject_lineages(s,c) = sum( subjects_all == my_subjects_list(s) & sample_clusters == my_lineages_list(c) );
    end
end

% Generate relative abundance table
my_table_subject_lineages_relativeabundance = my_table_subject_lineages./sum(my_table_subject_lineages,2);
my_table_subject_num_colonies = sum(my_table_subject_lineages,2);
my_table_subject_num_colonies_pores = sum(my_table_subject_lineages_pores,2);
my_table_subject_num_colonies_surface = sum(my_table_subject_lineages_surface,2);


%%

% Filters
min_col_per_subject = 20;

% Figure prep
% Load SLST colors
load( [ 'data/fig_format_slst_colors.mat' ] )
% Determine lineage colors
my_lineage_list_colors = zeros( numel(my_lineages_list), 3 );
for c=1:numel(my_lineages_list)
    % Get SLST of the lineage
    next_lineage = my_lineages_list(c);
    if next_lineage == 0 % unclustered
        next_lineage_slst = 'XX';
    else % cluster
        next_lineage_slst = clusters_all_slst{ my_lineages_list(c) };
    end
    next_lineage_slst_short = next_lineage_slst(1);
    my_lineage_list_colors(c,:) = slst_list_colors_short( ismember( slst_list_short, next_lineage_slst_short ), : );
end
% Sort subjects by number of colonies
[~, sort_subjects] = sort( my_table_subject_num_colonies, 'descend' );
% Filter subjects
sort_subjects = sort_subjects( my_table_subject_num_colonies(sort_subjects) >= min_col_per_subject );
% Sort lineages by SLST
my_lineage_list_slst = clusters_all_slst; my_lineage_list_slst{end+1} = 'XX';
[~, sort_lineages] = sort( my_lineage_list_slst );


%% STRAIN ABUNDANCE PER SUBJECT

% Implementation notes
% % Same deal as above but collapse lineages of same strain type into one
% bar on the bar chart
% % Does not include strain type of unclustered samples (would need to get 
% SLST via assemblies)

% Figure prep
my_table_subject_strains_relativeabundance = zeros( numel(my_subjects_list), numel(slst_list_short) );
for c=1:numel(my_lineages_list)
    % Get SLST of the lineage
    next_lineage = my_lineages_list(c);
    if next_lineage == 0 % unclustered
        next_lineage_slst = 'XX';
    else % cluster
        next_lineage_slst = clusters_all_slst{ my_lineages_list(c) };
    end
    next_lineage_slst_short = next_lineage_slst(1);
    next_lineage_slst_index = find( ismember( slst_list_short, next_lineage_slst_short) );
    my_table_subject_strains_relativeabundance(:,next_lineage_slst_index) = my_table_subject_strains_relativeabundance(:,next_lineage_slst_index) + my_table_subject_lineages_relativeabundance(:,c);
end


%%

% Implementation notes
% % Same deal as above but collapse lineages of same strain type into one
% bar on the bar chart
% % Does not include strain type of unclustered samples (would need to get 
% SLST via assemblies) in top bar chart
% % Shows proportion of unclustered samples in bottom bar chart

legend_bool = false;

% Bar chart
fig=figure(3);
clf(3)
fs=24;
lw=1.5;
%
% Number of colonies
subplot(6,5,[1:1:10])
hold on 
box on
my_table_subject_num_colonies_unclus = arrayfun(@(x) sum(sample_clusters(x==subjects_all)==0), my_subjects_list');
my_table_subject_num_colonies_clus = my_table_subject_num_colonies - my_table_subject_num_colonies_unclus;
my_table_subject_num_colonies_combined = [my_table_subject_num_colonies_clus(sort_subjects), my_table_subject_num_colonies_unclus(sort_subjects)];
axis_max = 200;
axis_jump = my_table_subject_num_colonies(1)-axis_max;
my_table_subject_num_colonies_combined_plot = my_table_subject_num_colonies_combined; 
my_table_subject_num_colonies_combined_plot(1,1) = my_table_subject_num_colonies_combined_plot(1,1)-axis_jump;
b = bar( my_table_subject_num_colonies_combined_plot, 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', lw );
b(1).CData = 0.25*[ 1 1 1 ];
b(2).CData = 0.75*[ 1 1 1 ];
% Legend
if legend_bool
    legend({'clustered','unclustered'},'Location','northeastoutside')
end
% Axes
ax = gca;
ax.YAxis.FontSize = fs-4;
% X
xticks( 1:1:numel(my_subjects_list) )
xticklabels( {} )
% Y
ylabel({'number','of colonies'},'FontSize',fs)
%ylabel('# colonies')
ylim([0 axis_max])
yticks([0:50:axis_max-100, axis_max-39, axis_max ])
temp_ytick_labels = [0:50:axis_max-100 max(my_table_subject_num_colonies)-39 max(my_table_subject_num_colonies)];
temp_ytick_labels = arrayfun(@(x) {num2str(x)}, temp_ytick_labels );
yticklabels( temp_ytick_labels )
%%%yticklabels([0:50:axis_max-50 max(my_table_subject_num_colonies)])
% Formatting
hold off
%
% Relative abundances of strains colored by SLST
subplot(6,5,[11:1:30])
hold on
box on
my_table_subject_strains_relativeabundance_nounclus = my_table_subject_strains_relativeabundance(:,1:end-1);
my_table_subject_strains_relativeabundance_nounclus = my_table_subject_strains_relativeabundance_nounclus./sum(my_table_subject_strains_relativeabundance_nounclus,2);
b = bar(my_table_subject_strains_relativeabundance_nounclus( sort_subjects, : ), 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', lw );
for t = 1:numel(slst_list_short)-1
    b(t).CData = slst_list_colors_short(t,:);
end
% Legend
if legend_bool
    leg=legend(slst_list_short(1:end-1),'Location','northeastoutside');
    leg.Title.String = 'strain type';
end
% Axes
ax = gca;
ax.YAxis.FontSize = fs-4;
ax.XAxis.FontSize = fs;
% X
xlim([0 numel(sort_subjects)+1])
xticks( 1:1:numel(sort_subjects) )
xticklabels( 1:1:numel(sort_subjects) )
%xticklabels( char(my_subjects_list(sort_subjects)+64) )
xlabel('subject','FontSize',fs)
% Y
yticks([0:0.5:1])
ylabel({'relative abundance', 'of strain type'},'FontSize',fs)
hold off

% Save fig
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 8]);
if legend_bool
    print(fig,[ '3_straintype_biogeo' '/Bar_StrainAbundance_Unclus_leg.png' ],'-dpng')
else
    print(fig,[ '3_straintype_biogeo' '/Bar_StrainAbundance_Unclus.png' ],'-dpng')
end



%% LINEAGE/SLST ABUNDANCES BY SUBJECT 

% 2B: Shows number of lineages per superSLST per subject

% Generate data table: number of lineages in each strain by subject
my_table_subject_strains = zeros( numel(my_subjects_list), numel(slst_list_short)-1 );
for c=1:numel(clusters_all)
    next_subject = find( char(my_subjects_list+64)==clusters_all_subjects{c} );
    next_strain = find( slst_list_short==clusters_all_slst_super(c) );
    my_table_subject_strains(next_subject,next_strain) = my_table_subject_strains(next_subject,next_strain) + 1;
end

%%

% Skip strain-type E beacuse no lineages
slst_list_filter = slst_list_short(1:end-1)~='E';
slst_list_filter_names = slst_list_short~='E'& slst_list_short~='X';

% Figure
fig=figure(10);
clf(10);
heatmap_text = my_table_subject_strains(sort_subjects,slst_list_filter)';
heatmap_vals = heatmap_text;
imagesc( heatmap_vals )
cmap = flipud(gray(5));
colormap( cmap(1:4,:) );
caxis([-0.5 3.5])
cbar = colorbar('XTick', 0:1:3);
cbar.YLabel.String = 'number of lineages';
% x axis
xticks( 1:1:numel(sort_subjects) )
%xticklabels( char(64+my_subjects_list(sort_subjects)) )
xlabel('subject')
% y axis
yticklabels( slst_list_short(slst_list_filter_names) )
ylabel('strain')
% gridlines
lw=1.5;
for i=1:size(heatmap_text,2)+1
    line( [i-0.5 i-0.5], ylim, 'Color', 'k', 'LineWidth', lw )
end
for j=1:size(heatmap_text,1)+1
    line( xlim, [j-0.5 j-0.5], 'Color', 'k', 'LineWidth', lw )
end
% general
set(gca,'FontSize',fs)

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 5]);
print(fig,[ '3_straintype_biogeo' '/Heatmap_Strain_LineageSubject.png' ],'-dpng')

