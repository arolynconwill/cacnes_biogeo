%% Summary

% Makes supplemental figure outlining sampling


%% Directory

dir_tables = '8_tables';
if ~exist( dir_tables, 'dir' )
    mkdir( dir_tables )
end


%% Load data

subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
specimen_number_all = load('data/cluster_step_variables','specimen_number_final'); specimen_number_all = specimen_number_all.specimen_number_final;
times_all = load('data/cluster_step_variables','times_final'); times_all = times_all.times_final;

type_extract = 1;
type_porestrip = 2;
type_scrape = 3;

% Reindex subjects
subjects_list = unique(subjects_all);
subjects_list_num_cols = arrayfun(@(x) sum(x==subjects_all), subjects_list );
[~, sort_subjects] = sort( subjects_list_num_cols, 'descend' );
subjects_list_sorted = subjects_list( sort_subjects );


%% Make a table

fid = fopen([dir_tables '/supp_table_sampling.csv'],'w');
fprintf(fid,'subject_letter, subject_number, time_in_months, num_extract_samples, num_extract_colonies, num_porestrip_samples, num_porestrip_colonies, num_scrape_samples, num_scrape_colonies, \n');
for s=1:numel(subjects_list)
    this_subject = subjects_list_sorted(s);
    this_subject_times = unique(times_all( subjects_all==this_subject ));
    this_subject_time_min = min( this_subject_times );
    for t=1:numel(this_subject_times)
        % get info for this subject at this time
        this_time = this_subject_times(t);
        this_time_zeroed = this_time - this_subject_time_min;
        this_subject_this_time_bool = ( subjects_all==this_subject ) & ( times_all==this_time );
        num_extract_samples = numel(unique( specimen_number_all( ( types_all==type_extract ) & this_subject_this_time_bool ) ));
        num_porestrip_samples = numel(unique( specimen_number_all( ( types_all==type_porestrip ) & this_subject_this_time_bool ) ));
        num_scrape_samples = numel(unique( specimen_number_all( ( types_all==type_scrape ) & this_subject_this_time_bool ) ));
        num_extract_colonies = sum( ( types_all==type_extract ) & this_subject_this_time_bool );
        num_porestrip_colonies = sum( ( types_all==type_porestrip ) & this_subject_this_time_bool );
        num_scrape_colonies = sum( ( types_all==type_scrape ) & this_subject_this_time_bool );
        % add row to csv
        fprintf(fid, [ ...
            char( this_subject+64 ) ', ' ...
            num2str(s) ',' ... 
            num2str(this_time_zeroed) ', ' ...
            num2str(num_extract_samples) ', ' ...
            num2str(num_extract_colonies) ', ' ...
            num2str(num_porestrip_samples) ', ' ...
            num2str(num_porestrip_colonies) ', ' ...
            num2str(num_scrape_samples) ', ' ...
            num2str(num_scrape_colonies) ', ' '\n' ...
            ] );
    end
end
fclose(fid);





