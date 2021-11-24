%% ANALYZE GENE CLUSTERS FROM CDHIT OUTPUT


%% Summary

% This script identifies between-lineage gene content differences based on
% gene clusters from cd-hit. 


%% Directory setup

% Main directory:
dir_main = char(pwd);
path( dir_main,path );
% Directory for Lieberman Lab scripts:
dir_lab_scripts = '../../lab_scripts';
path(dir_lab_scripts,path);
% Directory for my scripts:
dir_scripts_aro = [dir_main '/miniscripts' ];
path(dir_scripts_aro,path);

% Where to find cd-hit output
dir_cdhit_output = 'input_cdhit_gene_clusters';

% Where to find data from assemblies
dir_data_assemblies = '../../data/data_assemblies/';


%% Load information on lineages

% Lineage names
load( '../../1_snv_analysis/2_snvs/cluster_names.mat' )
clusters_all_subjects_nums = cellfun(@(x) str2double(x), clusters_all_subjects_nums );
num_lineages = numel(cluster_names);

% Get superSLST
clusters_all_slst_super = cellfun(@(x) x(1), clusters_all_slst );

% Get SNP distances between pairs of lineages
temp_table = readtable('input_btwn_lineage_snv_distance_mat/dist_mat.csv');
clusters_dist_mat = cell2mat(table2cell( removevars(temp_table,{'Var1','Var55'}) ));


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE GENE CLUSTER SIMILARITY BETWEEN PAIRS OF LINEAGES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read cd-hit output 

% Initialize
gene_clusters_by_lineage = [];
num_gene_clusters = 0; 

% Note: gene clusters in output indexed starting at 0; gene clusters here
% indexed starting at 1!

% Read file line by line
pid_cutoff = 95;
fid = fopen( [ dir_cdhit_output '/lineages_all_db_' num2str(pid_cutoff) '.clstr' ] );
next_line = fgetl(fid);
while ischar(next_line)
    if contains( next_line, 'Cluster' )
        num_gene_clusters = num_gene_clusters+1;
        gene_clusters_by_lineage = [ gene_clusters_by_lineage; zeros(1,num_lineages) ];
    elseif contains( next_line, 'L_' )
        next_line_split = strsplit( next_line, '_' );
        this_lineage = str2double(next_line_split{2});
        gene_clusters_by_lineage( num_gene_clusters, this_lineage ) = 1;
    end
    next_line = fgetl(fid);
end


%% Compute pairwise similarity between clusters

% Initialize
pairwise_num_genes_same = zeros( num_lineages, num_lineages );
pairwise_num_genes_diff = zeros( num_lineages, num_lineages );
pairwise_fracrel_genes = zeros( num_lineages, num_lineages );
pairwise_fractot_genes = zeros( num_lineages, num_lineages );

for i=1:num_lineages
    for j=i:num_lineages
        % Computes a few metrics for comparing genes between pairs of lineages
        [ pair_num_same, pair_num_diff, pair_fracrel, pair_fractot ] = get_pairwise_gene_content( gene_clusters_by_lineage(:,i), gene_clusters_by_lineage(:,j) );
        pairwise_num_genes_same(i,j) = pair_num_same; pairwise_num_genes_same(j,i) = pair_num_same; % num gene clusters in common
        pairwise_num_genes_diff(i,j) = pair_num_diff; pairwise_num_genes_diff(j,i) = pair_num_diff; % num gene clusters that are unique to one lineage
        pairwise_fracrel_genes(i,j) = pair_fracrel; pairwise_fracrel_genes(j,i) = pair_fracrel; % num gene clusters in common / num total gene clusters
        pairwise_fractot_genes(i,j) = pair_fractot; pairwise_fractot_genes(j,i) = pair_fractot; % num gene clusters in common / num total gene clusters detected across all lineages
    end
end


%% Cladeogram
% https://www.mathworks.com/help/bioinfo/ref/clustergram.html

% Simple way to visualize gene cluster differences between lineages
% but this is slow

clustergram( gene_clusters_by_lineage, ...
    'ColumnLabels', arrayfun(@(x) {x}, clusters_all_slst_super ), 'ColumnLabelsRotate', 0, ...
    'Symmetric', false, 'Colormap', 'parula' )


%%

%%%%%%%%%%%%%%%%
% MAKE FIGURES %
%%%%%%%%%%%%%%%%

% Compares gene diffs to SNV diffs btwn lineages
% Compares gene cluster diffs between lineages on the same vs diff subjects
% Compares gene clusters in a lineage to gene clusters drawn at random


%% Directory for figs

dir_figs = 'output_figures';
if ~exist( dir_figs, 'dir' )
    mkdir( dir_figs )
end


%% Gene differences vs SNP differences: scatter plots

% Filters
my_pairs = [1:1:num_lineages]>[1:1:num_lineages]';
is_same_superslst = clusters_all_slst_super==clusters_all_slst_super';
is_same_subject = clusters_all_subjects_nums == clusters_all_subjects_nums';

% Fit
[ fitobject, gof ] = fit( log10(clusters_dist_mat( my_pairs )), log10(pairwise_num_genes_diff( my_pairs )), 'Poly1' );
coeff_values = coeffvalues(fitobject)
coeff_intervals = confint(fitobject)
rsquared = gof.rsquare % 0.86
% Fitted values
fittedX = linspace(1, 5, 100);
fittedY = polyval(coeff_values, fittedX);  

% Appearance
plot_alpha = 0.25;
plot_markersize = 100;
colors = lines(2);
color_1 = colors(1,:);
color_2 = colors(2,:);
lw = 1;

% Sort by superSLST same/diff
fig=figure(10);
clf(10)
hold on
box on
scatter( clusters_dist_mat( my_pairs & is_same_superslst ), pairwise_num_genes_diff( my_pairs & is_same_superslst ), plot_markersize, 'MarkerEdgeAlpha', plot_alpha, 'LineWidth', lw )
scatter( clusters_dist_mat( my_pairs & ~is_same_superslst ), pairwise_num_genes_diff( my_pairs & ~is_same_superslst ), plot_markersize, 'MarkerEdgeAlpha', plot_alpha, 'LineWidth', lw )
plot(10.^fittedX,10.^fittedY, 'Color', 0*[1 1 1], 'LineWidth', 1, 'HandleVisibility', 'off') 
xlim( [10 100000] )
ylim( [10 1000] )
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('# SNV differences')
ylabel('# gene cluster differences')
[~, icons] = legend({'same strain-type','different strain-type'}, 'Location', 'northwest', 'FontSize', 14 );
icons = findobj(icons,'Type','patch');
icons = findobj(icons,'Marker','none','-xor');
set(icons,'MarkerSize',10); % makes marker size in legend bigger
set(gca,'FontSize',20)
% text(2000,25,[ 'best fit line:' ],'FontSize',14)
% text(2000,20,[ 'slope: 0.50 (0.49-0.51)' ],'FontSize',14)
% text(2000,16,[ 'intercept: 0.65 (0.61-0.69)' ],'FontSize',14)
% text(2000,12.8,[ 'R^2: 0.86' ],'FontSize',14)
text(45,20,[ 'best fit line: log10(#clusters) = 0.50 * log10(#SNVs) + 0.65' ],'FontSize',14)
text(45,16,[ 'R^2: 0.86' ],'FontSize',14)
hold off

print(fig,[ dir_figs '/scatter_pairwise_strain_pid' num2str(pid_cutoff) '.png' ],'-r400','-dpng')


%% Gene differences: histograms

fig=figure(12);
clf(12)
subplot(2,1,1)
histogram( pairwise_num_genes_diff( my_pairs & is_same_subject ), 0:10:1000, 'FaceColor', color_1, 'FaceAlpha', .5 ) 
xlabel('num gene cluster differences')
ylim([0 10])
ylabel('num lineage pairs')
title('same subject')
set(gca,'FontSize',20)
subplot(2,1,2)
histogram( pairwise_num_genes_diff( my_pairs & ~is_same_subject ), 0:10:1000, 'FaceColor', color_1, 'FaceAlpha', .5) 
xlabel('num gene cluster differences')
ylim([0 50])
ylabel('num lineage pairs')
title('different subjects')
set(gca,'FontSize',20)

print(fig,[ dir_figs '/pairwise_genediffs_pid' num2str(pid_cutoff) '.png' ],'-r400','-dpng')

% Distributions different?
[h,p,ks2stat] = kstest2( ...
    pairwise_num_genes_diff( my_pairs & is_same_subject ), ...
    pairwise_num_genes_diff( my_pairs & ~is_same_subject ) ...
    ) % 0.7177


%% Genes found in a subject's lineages vs random draw of lineages

for s=1:10

    % Which subject
    my_subject = s;
    my_subject_lineages = find( clusters_all_subjects_nums==my_subject );
    my_subject_num_lineages = numel( my_subject_lineages );

    % Observed gene clusters
    my_subject_num_clusters_obs = sum( sum( gene_clusters_by_lineage( :,my_subject_lineages ), 2 ) > 0 );

    % Randomly sampled gene clusters
    num_trials = 10000;
    my_subject_num_clusters_rnd = arrayfun( @(x) ...
        sum( sum( gene_clusters_by_lineage( :,randperm( num_lineages, my_subject_num_lineages ) ), 2 ) > 0 ), ...
        1:1:num_trials );

    % Probability observed number of gene clusters is 
    prob_more_clusters_than_random = sum( my_subject_num_clusters_obs<my_subject_num_clusters_rnd )/num_trials;

    % Figure
    fig=figure(13);
    clf(13)
    hold on
    histogram( my_subject_num_clusters_rnd, 2000:50:3500 )
    line( [my_subject_num_clusters_obs my_subject_num_clusters_obs], ylim, 'LineWidth', 4)
    text( 3250, 0.9*max(ylim), num2str(prob_more_clusters_than_random) )
    xlabel('number of gene clusters')
    ylabel('num random resamplings of lineages')
    title(['Subject ' num2str(my_subject) ' (num lineages = ' num2str(my_subject_num_lineages) ')' ])
    set(gca,'FontSize',14)
    hold off

    print(fig,[ dir_figs '/gene-cluster-resampling_subj-' num2str(my_subject) '.png' ],'-r400','-dpng')

end

%% same, but only for subjects with >10 samples

fig=figure(14);
clf(14)
hold on

for s=1:3

    % Which subject
    my_subject = s;
    my_subject_lineages = find( clusters_all_subjects_nums==my_subject );
    my_subject_num_lineages = numel( my_subject_lineages );

    % Observed gene clusters
    my_subject_num_clusters_obs = sum( sum( gene_clusters_by_lineage( :,my_subject_lineages ), 2 ) > 0 );

    % Randomly sampled gene clusters
    num_trials = 100000;
    my_subject_num_clusters_rnd = arrayfun( @(x) ...
        sum( sum( gene_clusters_by_lineage( :,randperm( num_lineages, my_subject_num_lineages ) ), 2 ) > 0 ), ...
        1:1:num_trials );

    % Probability observed number of gene clusters is 
    prob_more_clusters_than_random = sum( my_subject_num_clusters_obs<my_subject_num_clusters_rnd )/num_trials

    % Figure
    subplot(3,1,s)
    bin_0 = 2000;
    bin_w = 50;
    bin_f = 3500;
    next_probs = 100*histcounts( my_subject_num_clusters_rnd, bin_0:bin_w:bin_f )/num_trials;
    bar( next_probs, 1, 'FaceColor', color_1, 'FaceAlpha', 0.5 )
%    histogram( my_subject_num_clusters_rnd, 2000:50:3500 )
    % x axis
    if s==3
        xlabel('number of gene clusters')
    end
    xlim([0 30])
    xticks([0:10:30])
    xticklabels([bin_0:bin_w*10:bin_f])
    % y axis
    if s==2
        ylabel('percent of resamplings')
    end
    ylim([0 30])
    % title
    %title(['subject ' num2str(my_subject) ' (num lineages = ' num2str(my_subject_num_lineages) ')' ])
    title(['Subject ' num2str(my_subject) ])
    set(gca,'FontSize',20)
    % line
    next_line = (my_subject_num_clusters_obs-bin_0)/bin_w;
    line( [next_line next_line], ylim, 'LineWidth', 4)
    % legend
    legend( { 'random resampling (all subjects)', 'observed' }, 'Location', 'northwest', 'FontSize', 14 )

end

hold off

print(fig,[ dir_figs '/gene-cluster-resampling_combined.png' ],'-r400','-dpng')


% prob_more_clusters_than_random =
% 
%     0.8504
% 
% 
% prob_more_clusters_than_random =
% 
%     0.2006
% 
% 
% prob_more_clusters_than_random =
% 
%     0.9847

    
%% same, but only for subjects with >10 samples, and only resample from these three subjects

fig=figure(15);
clf(15)
hold on

lineages_to_resample = find( ismember( clusters_all_subjects_nums, [1 2 3 ]) );

for s=1:3

    % Which subject
    my_subject = s;
    my_subject_lineages = find( clusters_all_subjects_nums==my_subject );
    my_subject_num_lineages = numel( my_subject_lineages );

    % Observed gene clusters
    my_subject_num_clusters_obs = sum( sum( gene_clusters_by_lineage( :,my_subject_lineages ), 2 ) > 0 );

    % Randomly sampled gene clusters
    num_trials = 100000;
    my_subject_num_clusters_rnd = arrayfun( @(x) ...
        sum( sum( gene_clusters_by_lineage( :,lineages_to_resample( randperm( numel(lineages_to_resample), my_subject_num_lineages ) ) ), 2 ) > 0 ), ...
        1:1:num_trials );

    % Probability observed number of gene clusters is 
    prob_more_clusters_than_random = sum( my_subject_num_clusters_obs<my_subject_num_clusters_rnd )/num_trials

    % Figure
    subplot(3,1,s)
    bin_0 = 2000;
    bin_w = 50;
    bin_f = 3500;
    next_probs = 100*histcounts( my_subject_num_clusters_rnd, bin_0:bin_w:bin_f )/num_trials;
    bar( next_probs, 1, 'FaceColor', color_1, 'FaceAlpha', 0.5 )
%    histogram( my_subject_num_clusters_rnd, 2000:50:3500 )
    % x axis
    if s==3
        xlabel('number of gene clusters')
    end
    xlim([0 30])
    xticks([0:10:30])
    xticklabels([bin_0:bin_w*10:bin_f])
    % y axis
    if s==2
        ylabel('percent of resamplings')
    end
    ylim([0 30])
    % title
    %title(['subject ' num2str(my_subject) ' (num lineages = ' num2str(my_subject_num_lineages) ')' ])
    title(['Subject ' num2str(my_subject) ])
    set(gca,'FontSize',20)
    % line
    next_line = (my_subject_num_clusters_obs-bin_0)/bin_w;
    line( [next_line next_line], ylim, 'LineWidth', 4)
    % legend
    legend( { 'random resampling (Subjects 1/2/3)', 'observed' }, 'Location', 'northwest', 'FontSize', 14 )

end

hold off

print(fig,[ dir_figs '/gene-cluster-resampling_combined_123only.png' ],'-r400','-dpng')


% prob_more_clusters_than_random =
% 
%     0.8370
% 
% 
% prob_more_clusters_than_random =
% 
%     0.0918
% 
% 
% prob_more_clusters_than_random =
% 
%     0.9894

%%

%%%%%%%%%%%%%%%%%%%
% MAKE SUPP TABLE %
%%%%%%%%%%%%%%%%%%%

%% Directory

dir_supp_tab = 'output_supp_table';
if ~exist( dir_supp_tab, 'dir' )
    mkdir( dir_supp_tab )
end


%% Get info for supplemental table

% Initialize
homolog_table = cell( num_gene_clusters, num_lineages );
gene_cluster_annotations = cell( num_gene_clusters, 2 );

% Re-read cd-hit output
fid = fopen( [ dir_cdhit_output '/lineages_all_db_' num2str(pid_cutoff) '.clstr' ] );
next_line = fgetl(fid);
gene_cluster_index = 0; 
while ischar(next_line)
    if contains( next_line, 'Cluster' )
        gene_cluster_index = gene_cluster_index+1;
        gene_clusters_by_lineage = [ gene_clusters_by_lineage; zeros(1,num_lineages) ];
    elseif contains( next_line, 'L_' )
        next_line_split = strsplit( next_line, '_' );
        next_line_split_2 = strsplit( next_line_split{4}, '...' );
        this_lineage = str2double(next_line_split{2});
        this_gene = [ next_line_split{3} '_' next_line_split_2{1} ];
        if isempty( homolog_table{ gene_cluster_index, this_lineage } )
            homolog_table{ gene_cluster_index, this_lineage } = this_gene;
        else
            homolog_table{ gene_cluster_index, this_lineage } = [ homolog_table{ gene_cluster_index, this_lineage } ' ' this_gene ];
        end
        if contains( next_line, '*' )
            gene_cluster_annotations{gene_cluster_index,1} = this_lineage;
            if isempty( gene_cluster_annotations{gene_cluster_index,2} )
                gene_cluster_annotations{gene_cluster_index,2} = this_gene;
            end
        end
    end
    next_line = fgetl(fid);
end


%% Get annotations for each gene cluster

gene_cluster_lineages = cell2mat( gene_cluster_annotations(:,1) );
gene_cluster_locustags = gene_cluster_annotations(:,2);

for i=1:num_lineages
   
    % Import fasta file
    my_lineage_fasta = fastaread( [ '../1_contig_gene_filtering/input_fastas_prokka/clade_' num2str(i) '_prokka_out.faa' ] );
    my_lineage_fasta_headers = {my_lineage_fasta.Header};
    
    % Locustags from this lineage
    my_indices = find( gene_cluster_lineages == i );
    for j=1:numel(my_indices)
        my_locustag = gene_cluster_locustags{ my_indices(j) };
        my_header = my_lineage_fasta_headers{ ...
            find( cellfun(@(x) contains( x,my_locustag ), my_lineage_fasta_headers ) ) ...
            };
        my_annotation = strtrim( extractAfter( my_header, my_locustag ) ); 
        % Save annontation
        gene_cluster_annotations{my_indices(j),3} = my_annotation;
    end
    
end


%% Clean up annotationns

for i=1:num_gene_clusters
    my_annotation_temp = gene_cluster_annotations{i,3};
    if contains( my_annotation_temp, ',' )
        gene_cluster_annotations{i,3} = strrep( my_annotation_temp, ',', ';' );
    end
end


%% Save

save([ dir_supp_tab '/supp_table_info.mat' ],'gene_cluster_annotations','homolog_table')


%% Write supplemental table

fid = fopen( [ dir_supp_tab '/supp_table_gene_clusters.csv' ], 'w' );

% Header row
fprintf( fid, 'gene_cluster_index,prokka_annotation,prokka_annotation_source_lineage,prokka_annotation_source_locustag,' );
for j=1:num_lineages
    fprintf( fid, [ cluster_names_new{j} ',' ] );
end
fprintf( fid, '\n' );

% One row per gene cluster
for i=1:num_gene_clusters
    fprintf( fid, [ num2str(i) ',' ...
        gene_cluster_annotations{i,3} ',' ...
        cluster_names_new{gene_cluster_annotations{i,1}} ',' gene_cluster_annotations{i,2} ',' ] );
    for j=1:num_lineages
        fprintf( fid, [ homolog_table{i,j} ',' ] );
    end
    fprintf( fid, '\n' );
end

