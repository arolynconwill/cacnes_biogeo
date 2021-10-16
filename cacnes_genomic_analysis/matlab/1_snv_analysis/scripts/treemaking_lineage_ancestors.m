%% SUMMARY

% This script prepares data for making a tree of lineage ancestors by
% finding ancestral calls for each lineage over all variable positions.


%% Set up directories and environment

% Directory for trees
dir_treemaking = '7_trees_lineageancestors';
if ~exist( dir_treemaking, 'dir' )
    mkdir( dir_treemaking )
end


%% Load lineage names and other basic info

dir_clusters = '2_snvs';
load( [ dir_clusters '/cluster_names.mat' ] )
num_clades = numel(cluster_names);

% Nucleotides: 1=A, 2=T, 3=C, 4=G
NTs='ATCG';


%% Get inferred ancestor for each clade

load([ 'data/cluster_step_variables.mat' ],'p_all');
Calls_clade_ancestors = zeros( numel(clusters_all), numel(p_all) ); % initialize
Clade_slst = cell( numel(clusters_all), 1); % initialize

for k=1:numel(clusters_all)

    % Load data for this cluster
    this_cluster_name = cluster_names{k};
    load([ dir_clusters '/' this_cluster_name '/' 'data_' this_cluster_name '.mat' ],'Calls_lineage_ancestor');
    load([ dir_clusters '/' this_cluster_name '/' 'data_' this_cluster_name '.mat' ],'slst');
    
    % Save inferred lineage ancestor
    Calls_clade_ancestors(k,:) = Calls_lineage_ancestor;
    
    % Get SLST 
    slst_list = unique(slst);
    slst_list_counts = zeros( numel(slst_list), 1 );
    for i=1:numel(slst_list)
       slst_list_counts(i) = length(find(strcmp(slst_list{i},slst)));
    end
    [ ~, this_clade_slst_index ] = max( slst_list_counts );
    this_clade_slst = slst_list{ this_clade_slst_index };
    % Record
    Clade_slst{k} = this_clade_slst;
    
end
Calls_clades_all_original = Calls_clade_ancestors;


%% Add SLST to clade names

% Clade SLST supergroup
Clade_slst_super = cellfun(@(x) x(1), Clade_slst);

% Add SLST to clade name...
names_tree_slst = {};
for i=1:numel(cluster_names)
    names_tree_slst{i} = [ cluster_names{i} '_SLST-' Clade_slst{i} ];
end

% Save
save([dir_treemaking '/' 'clade_slst_info'],'Clade_slst','Clade_slst_super','names_tree_slst')


%% Filter positions a bit

% Params
max_fraction_ambigious_samples = .1; % want positions where most clades have calls

% Find nonvariable positions
pos_nonvar_1 = ( sum( Calls_clade_ancestors == 1 ) ==  sum( Calls_clade_ancestors ~= 0 ) );
pos_nonvar_2 = ( sum( Calls_clade_ancestors == 2 ) ==  sum( Calls_clade_ancestors ~= 0 ) );
pos_nonvar_3 = ( sum( Calls_clade_ancestors == 3 ) ==  sum( Calls_clade_ancestors ~= 0 ) );
pos_nonvar_4 = ( sum( Calls_clade_ancestors == 4 ) ==  sum( Calls_clade_ancestors ~= 0 ) );
pos_nonvar = pos_nonvar_1 + pos_nonvar_2 + pos_nonvar_3 + pos_nonvar_4;
pos_nonvar = ( pos_nonvar > 0 );
num_pos_nonvar = sum( pos_nonvar );

% Find ambiguous positions
pos_ambig = ( sum( Calls_clade_ancestors == 0, 1)/num_clades > max_fraction_ambigious_samples );
num_pos_ambig = sum( pos_ambig );

% Downsize calls matrix
pos_remove = pos_ambig | pos_nonvar;
num_pos_remove = sum( pos_remove ); 

% Downsize
Calls_clades_all_smaller = Calls_clade_ancestors(:,~pos_remove); 
p_all_clades = p_all( ~pos_remove );


%% Distance matrix (for scale bars)

% Generate distance matrix
dist_mat = zeros( numel(clusters_all), numel(clusters_all) );
for i=1:numel(clusters_all)
    for j=1:numel(clusters_all)
        dist_mat(i,j) = sum( Calls_clades_all_smaller(i,:)~=Calls_clades_all_smaller(j,:) & Calls_clades_all_smaller(i,:)>0 & Calls_clades_all_smaller(j,:)>0 );
    end
end

% Write csv file
fid = fopen( [ dir_treemaking '/' 'dist_mat.csv' ], 'w' );
% first line
fprintf(fid, ',');
for i=1:numel(clusters_all)
    fprintf(fid, [ cluster_names{i} '/' cluster_names_new{i}  ',' ] );
end
fprintf(fid,'\n');
% next lines
for i=1:numel(clusters_all)
    fprintf(fid, [ cluster_names{i} '/' cluster_names_new{i} ',' ] );
    for j=1:numel(clusters_all)
        fprintf(fid, [ num2str(dist_mat(i,j)) ',' ] );
    end
    fprintf(fid,'\n');
end
fclose(fid);


%% Make tree

calls_for_tree = zeros(size(Calls_clades_all_smaller));
calls_for_tree(Calls_clades_all_smaller>0) = NTs(Calls_clades_all_smaller(Calls_clades_all_smaller>0));
calls_for_tree(Calls_clades_all_smaller<1) = '?';
names_tree = names_tree_slst;


%% Save to make on cluster

dir_run_on_cluster = [ dir_treemaking '/run_on_cluster' ];
if ~exist( dir_run_on_cluster, 'dir' )
    mkdir( dir_run_on_cluster )
end
if ~exist( [ dir_run_on_cluster '/myscripts' ], 'dir' )
    mkdir( [ dir_run_on_cluster '/myscripts' ] )
end
save([ dir_run_on_cluster '/' 'clade_tree_cluster_data.mat'],'calls_for_tree','names_tree');
! cp scripts/myscripts_treemaking/dnapars 7_trees_lineageancestors/run_on_cluster/dnapars
! cp scripts/myscripts_treemaking/make_tree_clade_cluster.m 7_trees_lineageancestors/run_on_cluster/make_tree_clade_cluster.m
! cp scripts/myscripts_treemaking/myjob.slurm 7_trees_lineageancestors/run_on_cluster/myjob.slurm
! cp scripts/myscripts_treemaking/myscripts/generate_parsimony_tree_old.m 7_trees_lineageancestors/run_on_cluster/myscripts/generate_parsimony_tree_old.m
! cp scripts/myscripts_treemaking/myscripts/generate_phylip_input_old.m 7_trees_lineageancestors/run_on_cluster/myscripts/generate_phylip_input_old.m
! cp scripts/myscripts_treemaking/myscripts/relabeltree_old.m 7_trees_lineageancestors/run_on_cluster/myscripts/relabeltree_old.m
! cp scripts/myscripts_treemaking/myscripts/sort_cluster_order.m 7_trees_lineageancestors/run_on_cluster/myscripts/sort_cluster_order.m

