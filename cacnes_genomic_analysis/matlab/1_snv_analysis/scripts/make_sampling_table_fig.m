%% Summary

% Makes main figure outlining sampling


%% Directory

dir_tables = '8_tables';
if ~exist( dir_tables, 'dir' )
    mkdir( dir_tables )
end


%% Load data

subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
specimen_number_all = load('data/cluster_step_variables','specimen_number_final'); specimen_number_all = specimen_number_all.specimen_number_final;

load('data/spec_mult.mat')

type_extract = 1;
type_porestrip = 2;
type_scrape = 3;

% Reindex subjects
subjects_list = unique(subjects_all);
subjects_list_num_cols = arrayfun(@(x) sum(x==subjects_all), subjects_list );
[~, sort_subjects] = sort( subjects_list_num_cols, 'descend' );
subjects_list_sorted = subjects_list( sort_subjects );


%% Make a table

fid = fopen([dir_tables '/main_table_sampling.csv'],'w');
fprintf(fid,'subject_letter,subject_number,num_colonies,num_scrape_samples_multcol,num_scrape_samples_singcol,num_pore_singpore_multcol,num_pore_singpore_singcol,num_pore_multpore_multcol,num_pore_multpore_singcol, \n');
for s=1:numel(subjects_list)
    this_subject = subjects_list_sorted(s);
    this_subject_bool = ( subjects_all==this_subject );
    % number of colonies
    num_colonies = sum( this_subject_bool );
    % scrape samples
    scrape_samples = unique( specimen_number_all( ( types_all==type_scrape ) & this_subject_bool ) );
    scrape_num = numel( scrape_samples );
    scrape_samples_cols = arrayfun(@(x) sum(specimen_number_all==x), scrape_samples );
    num_scrape_singcol = sum( scrape_samples_cols==1 );
    num_scrape_multcol = sum( scrape_samples_cols>1 );
    % extracts and strip samples
    pore_samples = unique( specimen_number_all( ( types_all==type_extract | types_all==type_porestrip ) & this_subject_bool ) );
    pore_num = numel( pore_samples );
    pore_samples_cols = arrayfun(@(x) sum(specimen_number_all==x), pore_samples );
    pore_samples_sing = arrayfun(@(x) spec_mult_porenums(spec_mult_specnums==x)==1, pore_samples );
    num_pore_singpore_singcol = sum( pore_samples_cols==1 & pore_samples_sing==1 );
    num_pore_singpore_multcol = sum( pore_samples_cols>1 & pore_samples_sing==1 );
    num_pore_multpore_singcol = sum( pore_samples_cols==1 & pore_samples_sing==0 );
    num_pore_multpore_multcol = sum( pore_samples_cols>1 & pore_samples_sing==0 );
    % add row to csv
    fprintf(fid, [ ...
        char( this_subject+64 ) ', ' ...
        num2str(s) ',' ... 
        num2str(num_colonies) ', ' ...
        num2str(num_scrape_multcol) ', ' ...
        num2str(num_scrape_singcol) ', ' ...
        num2str(num_pore_singpore_multcol) ', ' ...
        num2str(num_pore_singpore_singcol) ', ' ...
        num2str(num_pore_multpore_multcol) ', ' ...
        num2str(num_pore_multpore_singcol) ', ' '\n' ...
        ] );
end
fclose(fid);





