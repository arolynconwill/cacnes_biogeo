%% SUMMARY

% This script makes a supplemental table with inforamtion on each lineage


%% Directory

dir_tables = '8_tables';
if ~exist( dir_tables, 'dir' )
    mkdir( dir_tables )
end

%% Directory setup

path(path,'scripts/myscripts_denovomuts')


%% Load lineage info

% Load lineage names and SLST info
dir_clusters = '2_snvs';
load( [ dir_clusters '/'  'cluster_names.mat' ]) % lineage names
cluster_sizes = cellfun(@(x) numel(x), clusters_all);
load( [ dir_clusters '/' 'sample_names.mat' ]) 
slst_all = cellfun(@(x) x(end-7:end-6), SampleNamesLong_all, 'UniformOutput', false );
for i=1:numel(slst_all)
    temp = slst_all{i};
    if temp(1)=='-'
        slst_all{i} = temp(2);
    else
        slst_all{i} = temp(1);
    end
end
clusters_all_slst = cellfun(@(x) mode(cell2mat(slst_all(x))), clusters_all );

subjects_list_nums = unique( cellfun(@(x) str2double(x), clusters_all_subjects_nums ) );


%% Write a table

fid = fopen( [ dir_tables '/supp_table_lineageinfo.csv' ], 'w' );
fprintf( fid, [ 'lineage_name, subject, strain_type, num_colonies,\n'] );
for s=1:numel(subjects_list_nums)
    this_subject_clusters = find( cellfun(@(x) str2double(x), clusters_all_subjects_nums) == subjects_list_nums(s) );
    for i=1:numel(this_subject_clusters)
        fprintf( fid, [ cluster_names_new{this_subject_clusters(i)} ',' clusters_all_subjects_nums{this_subject_clusters(i)} ',' clusters_all_slst(this_subject_clusters(i)) ',' num2str(numel(clusters_all{this_subject_clusters(i)})) ',\n'] );
    end
end
fclose(fid);




