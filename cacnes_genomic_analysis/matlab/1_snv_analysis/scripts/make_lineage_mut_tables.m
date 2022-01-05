%% SUMMARY

% This script makes lineage trees with easily interpretable sample names.


%% Directory

dir_save = '8_tables';
if ~exist( dir_save, 'dir' )
    mkdir( dir_save )
end

%% Set up directories and environment

% Add my scripts
dir_my_scripts = [ pwd '/' 'scripts/myscripts_misc'];
path(dir_my_scripts,path);

% Nucleotides: 1=A, 2=T, 3=C, 4=G
NTs = 'ATCG';


%% Load lineage info

% Load lineage names
dir_clusters = '2_snvs';
load( [ dir_clusters '/'  'cluster_names.mat' ],'cluster_names','cluster_names_new','clusters_all','clusters_all_subjects') % lineage names
lineage_names = cluster_names; clear cluster_names; 
lineage_names_new = cluster_names_new; clear cluster_names_new;
cluster_sizes = cellfun(@(x) numel(x), clusters_all);

% List of path to files containing annotation_full for each lineage
paths_to_files = cellfun(@(lineage_name) [ dir_clusters '/' lineage_name '/' 'data_' lineage_name '.mat' ], lineage_names, 'UniformOutput', false );

% Load coverage data
load( 'data_other/chromosomal_coverage.mat' )

% Reindex subjects
subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
subjects_list = unique(subjects_all);
subjects_list_num_cols = arrayfun(@(x) sum(x==subjects_all), subjects_list );
[~, sort_subjects] = sort( subjects_list_num_cols, 'descend' );
subjects_list_sorted = subjects_list( sort_subjects );


%% Other info

% metadata keys
% sample type
type_key={'pore_extract', 'pore_strip', 'scrape'};
% zone name
zone_key={{'forehead','Fo','Tz'},{'chin','Ch','Cn','Ce'},{'right_cheek','Rc'},{'left_cheek','Lc'},{'nose','No','no'},{'near_mouth','Mo','nc'},{'back','Ba'},{'neck','Nk'},{'shoulder','Sh'}};
zone_key = cellfun(@(x) x{1}, zone_key, 'UniformOutput', false);

% plasmid presence/absence
load('data/data_plasmid_presence.mat','has_plasmid')
load('data/data_plasmid_presence.mat','has_transposon')
load('data/data_plasmid_presence.mat','SampleNamesSimple_keep'); SampleNamesSimple_plasmid = SampleNamesSimple_keep; clear SampleNamesSimple_keep;

% read lengths
load('data/read_lengths.mat')


%% Loop through lineages

for c=1:numel(lineage_names)
    
    % load lineage data
    this_lineage_name = lineage_names{c};
    this_lineage_name_new = lineage_names_new{c};
    load([dir_clusters '/' this_lineage_name '/data_' this_lineage_name '.mat']);
    lineage_num_str = num2str(c);
    if numel(lineage_num_str) == 1
        lineage_num_str = [ '0' lineage_num_str ];
    end
    % pre-processing
    ingroup_indices = find( ~outgroup_isolates );
    temp1=cellfun(@(x) strsplit(x,'_'), SampleNamesSimple, 'UniformOutput', false);
    temp2=cellfun(@(x) x{2}, temp1, 'UniformOutput', false);
    temp3=cellfun(@(x) strsplit(x,'-'), temp2, 'UniformOutput', false);
    SampleNamesSimple_multiplicity = cellfun(@(x) x{1}, temp3, 'UniformOutput', false);
    temp4=cellfun(@(x) x{3}, temp1, 'UniformOutput', false);
    temp5=cellfun(@(x) strsplit(x,'-'), temp4, 'UniformOutput', false);
    SampleNamesSimple_colnum = cellfun(@(x) x{2}, temp5, 'UniformOutput', false);
    clear temp1; clear temp2; clear temp3; clear temp4; clear temp5;

    %%

    % get sample metadata
    sample_metadata_list = { 'sample_name_regular', 'sample_name_simple', 'subject', ...
        'sampling_method', 'pore_multiplicity', 'sample_number', 'colony_number', 'sebaceous_skin_region', 'time_in_months', ...
        'median_coverage', 'read_length', 'plasmid_presence', 'transposon_presence' };
    sample_metadata_count = numel( sample_metadata_list );

    sample_metadata = struct; % initialize
    for i=1:numel(ingroup_indices)
        next_index = ingroup_indices(i);
        next_simple_name = SampleNamesSimple{next_index};
        next_name = SampleNames{next_index};
        % populate structure
        sample_metadata(next_index).sample_name_regular = SampleNames{next_index};
        sample_metadata(next_index).sample_name_simple = next_simple_name;
        sample_metadata(next_index).subject = num2str(find(subjects(next_index)==subjects_list_sorted)); % char( subjects(next_index)+64 );
        sample_metadata(next_index).sampling_method = type_key{types(next_index)};
        sample_metadata(next_index).pore_multiplicity = SampleNamesSimple_multiplicity{next_index};
        sample_metadata(next_index).sample_number = num2str(specimen_numbers(next_index));
        sample_metadata(next_index).colony_number = SampleNamesSimple_colnum{next_index};
        sample_metadata(next_index).sebaceous_skin_region = zone_key{zones(next_index)};
        sample_metadata(next_index).time_in_months = num2str(times(next_index));
        sample_metadata(next_index).median_coverage = num2str(coverage_chromosome_median(ismember(SampleNames_all,SampleNames(next_index))));
        % read length
        sample_metadata(next_index).read_length = num2str(read_lengths( ismember( SampleNames_rl, next_name ) ));
        % plasmid presence
        sample_metadata(next_index).plasmid_presence = num2str(has_plasmid( ismember( SampleNamesSimple_plasmid, next_simple_name ) ));
        % transposon presence
        sample_metadata(next_index).transposon_presence = num2str(has_transposon( ismember( SampleNamesSimple_plasmid, next_simple_name ) ));
    end

    if numel(goodpos) > 0
        
        % get SNV metadata
        snv_metadata_list = { 'position_on_ref_genome', 'protein_id', 'locus_tag', 'annontation', 'mutation_type', 'mutation_info', 'ancestral_allele' };
        snv_metadata_count = numel( snv_metadata_list );

        if numel(goodpos) > 1
            snv_metadata = struct;
            for j=1:length(annotation_full)
                snv_metadata(j).position_on_ref_genome = num2str(annotation_full(j).pos);
                snv_metadata(j).protein_id = annotation_full(j).protein_id;
                snv_metadata(j).locus_tag = annotation_full(j).locustag;
                snv_metadata(j).annontation = parse_gene_product( annotation_full(j).annotation );
                if isequal( annotation_full(j).type, 'P' )
                    snv_metadata(j).mutation_type = 'I';
                else
                    snv_metadata(j).mutation_type = annotation_full(j).type;
                end
                temp_muts = annotation_full(j).muts;
                if ~isempty( temp_muts )
                    snv_metadata(j).mutation_info = temp_muts{1};
                else
                    snv_metadata(j).mutation_info = '';
                end
                snv_metadata(j).ancestral_allele = NTs(anc_nti_goodpos(j));
            end
        else % only one SNV
            snv_metadata(j).position_on_ref_genome = num2str(annotation_full.pos);
            if ismember( 'protein_id', fieldnames(annotation_full) )
                snv_metadata(j).protein_id = annotation_full.protein_id;
            else
                snv_metadata(j).protein_id = '';
            end
            if ismember( 'locus_tag', fieldnames(annotation_full) )
                snv_metadata(j).locus_tag = annotation_full.locustag;
            else
                snv_metadata(j).locus_tag = '';
            end
            snv_metadata(j).annontation = parse_gene_product( annotation_full.annotation );
            if isequal( annotation_full.type, 'P' )
                snv_metadata(j).mutation_type = 'I';
            else
                snv_metadata(j).mutation_type = annotation_full.type;
            end
            if ismember( 'muts', fieldnames(annotation_full) )
                snv_metadata(j).mutation_info = temp_muts{1};
            else
                snv_metadata(j).mutation_info = '';
            end
            snv_metadata(j).ancestral_allele = NTs(anc_nti_goodpos);
        end

        % get SNV alleles
        snv_table = Calls_for_analysis(:,~outgroup_isolates);
        snv_table( snv_table > 0 ) = NTs( snv_table( snv_table > 0 ) );
        snv_table( snv_table == 0 ) = 'N';
        
    end
    
    
    %% build CSV: old version with sample metadata and SNV metadata

    fid = fopen( [ dir_save '/' 'temp_supp_table_' lineage_num_str '_everything.csv' ], 'w' );
    % lineage name
    column_buffer_0 = char;
    for i=1:snv_metadata_count+numel(ingroup_indices)
        column_buffer_0(end+1) = ',';
    end
    fprintf( fid, [ this_lineage_name_new ',' column_buffer_0 '\n'] );
    % sample metadata
    column_buffer_1 = char;
    for i=1:snv_metadata_count
        column_buffer_1(end+1) = ',';
    end
    for r=1:sample_metadata_count
        fprintf( fid, column_buffer_1 );
        fprintf( fid, [ sample_metadata_list{r} ',' ] );
        for n=1:numel(ingroup_indices)
            fprintf( fid, [ sample_metadata(n).(sample_metadata_list{r}) ','] );
        end
        fprintf( fid, '\n' );
    end
    % SNV metadata headers
    for m=1:snv_metadata_count
        fprintf( fid, [ snv_metadata_list{m} ',' ] );
    end
    column_buffer_2 = char;
    for i=1:numel(ingroup_indices)+1
        column_buffer_2(end+1) = ',';
    end
    fprintf( fid, [ column_buffer_2 '\n' ] );
    if numel(goodpos) > 0 
        % SNV metadata and calls
        for r=1:numel(goodpos)
            for m=1:snv_metadata_count
                fprintf( fid, [ snv_metadata(r).(snv_metadata_list{m}) ',' ] );
            end
            fprintf( fid, ',' );
            for n=1:numel(ingroup_indices)
                fprintf( fid, [ snv_table(r,n) ','] );
            end    
            fprintf( fid, '\n' );
        end
    else % no SNVs
        % no SNV metadata and calls
        fprintf( fid, column_buffer_1 );
        fprintf( fid, ',' );
        fprintf( fid, 'no SNVs detected in this lineage,' );
        for n=1:numel(ingroup_indices)-1
            fprintf( fid, ',' );
        end    
        fprintf( fid, '\n' );
    end
    fclose(fid);

    % check no issues with csv
    %readtable( [ dir_save '/' 'supp_table_' lineage_names{c} '.csv' ] )

    
    %% build CSV: new version SNV table w/o sample metadata
    
    fid = fopen( [ dir_save '/' 'temp_supp_table_' lineage_num_str '_snvtable.csv' ], 'w' );
    % lineage name
    column_buffer_0 = char;
    for i=1:snv_metadata_count+numel(ingroup_indices)
        column_buffer_0(end+1) = ',';
    end
    fprintf( fid, [ this_lineage_name_new ',' column_buffer_0 '\n'] );
    % sample metadata
    column_buffer_1 = char;
    for i=1:snv_metadata_count
        column_buffer_1(end+1) = ',';
    end
    for r=2%1:sample_metadata_count
        fprintf( fid, column_buffer_1 );
        fprintf( fid, [ sample_metadata_list{r} ',' ] );
        for n=1:numel(ingroup_indices)
            fprintf( fid, [ sample_metadata(n).(sample_metadata_list{r}) ','] );
        end
        fprintf( fid, '\n' );
    end
    % SNV metadata headers
    for m=1:snv_metadata_count
        fprintf( fid, [ snv_metadata_list{m} ',' ] );
    end
    column_buffer_2 = char;
    for i=1:numel(ingroup_indices)+1
        column_buffer_2(end+1) = ',';
    end
    fprintf( fid, [ column_buffer_2 '\n' ] );
    if numel(goodpos) > 0 
        % SNV metadata and calls
        for r=1:numel(goodpos)
            for m=1:snv_metadata_count
                fprintf( fid, [ snv_metadata(r).(snv_metadata_list{m}) ',' ] );
            end
            fprintf( fid, ',' );
            for n=1:numel(ingroup_indices)
                fprintf( fid, [ snv_table(r,n) ','] );
            end    
            fprintf( fid, '\n' );
        end
    else % no SNVs
        % no SNV metadata and calls
        fprintf( fid, column_buffer_1 );
        fprintf( fid, ',' );
        fprintf( fid, 'no SNVs detected in this lineage,' );
        for n=1:numel(ingroup_indices)-1
            fprintf( fid, ',' );
        end    
        fprintf( fid, '\n' );
    end
    fclose(fid);
    
    
    %% build CSV: new version colony metadata only
    
    fid = fopen( [ dir_save '/' 'temp_supp_table_' lineage_num_str '_colonytable.csv' ], 'w' );
    % lineage name
    fprintf( fid, [ this_lineage_name_new ',' column_buffer_0 '\n'] );
    % sample metadata headers
    for r=2:sample_metadata_count
        fprintf( fid, [ sample_metadata_list{r} ',' ] );
    end
    fprintf( fid, '\n' );
    % sample metadata
    for n=1:numel(ingroup_indices)
        for r=2:sample_metadata_count
             fprintf( fid, [ sample_metadata(n).(sample_metadata_list{r}) ','] );
        end
        fprintf( fid, '\n' );
    end
    fclose(fid);

end


%% Also make a table for unclustered samples

% Load data
SampleNames_all = load('2_snvs/sample_names','SampleNames_all'); SampleNames_all = SampleNames_all.SampleNames_all;
SampleNamesSimple_all = load('2_snvs/sample_names','SampleNamesSimple_all'); SampleNamesSimple_all = SampleNamesSimple_all.SampleNamesSimple_all;
unclustered_all = load('data/cluster_step_variables','unclustered_final'); unclustered_all = unclustered_all.unclustered_final;
subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
zones_all = load('data/cluster_step_variables','zones_final'); zones_all = zones_all.zones_final;
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
times_all = load('data/cluster_step_variables','times_final'); times_all = times_all.times_final;
specimen_number_all = load('data/cluster_step_variables','specimen_number_final'); specimen_number_all = specimen_number_all.specimen_number_final;

this_lineage_name = 'unclustered_colonies';
SampleNames = SampleNames_all(unclustered_all);
SampleNamesSimple = SampleNamesSimple_all(unclustered_all);
temp1=cellfun(@(x) strsplit(x,'_'), SampleNamesSimple, 'UniformOutput', false);
temp2=cellfun(@(x) x{2}, temp1, 'UniformOutput', false);
temp3=cellfun(@(x) strsplit(x,'-'), temp2, 'UniformOutput', false);
SampleNamesSimple_multiplicity = cellfun(@(x) x{1}, temp3, 'UniformOutput', false);
temp4=cellfun(@(x) x{3}, temp1, 'UniformOutput', false);
temp5=cellfun(@(x) strsplit(x,'-'), temp4, 'UniformOutput', false);
SampleNamesSimple_colnum = cellfun(@(x) x{2}, temp5, 'UniformOutput', false);
clear temp1; clear temp2; clear temp3; clear temp4; clear temp5;
subjects = subjects_all(unclustered_all);
types = types_all(unclustered_all);
specimen_numbers = specimen_number_all(unclustered_all);
zones = zones_all(unclustered_all);
times = times_all(unclustered_all);

[ ~,reorder ] = sort( SampleNames );

sample_metadata = struct; % initialize
for i=1:numel(unclustered_all)
    next_index = reorder(i);
    next_simple_name = SampleNamesSimple{next_index};
    next_name = SampleNames{next_index};
    % populate structure
    sample_metadata(i).sample_name_regular = SampleNames{next_index};
    sample_metadata(i).sample_name_simple = next_simple_name;
    sample_metadata(i).subject = num2str(find(subjects(next_index)==subjects_list_sorted)); %char( subjects(next_index)+64 );
    sample_metadata(i).sampling_method = type_key{types(next_index)};
    sample_metadata(i).pore_multiplicity = SampleNamesSimple_multiplicity{next_index};
    sample_metadata(i).sample_number = num2str(specimen_numbers(next_index));
    sample_metadata(i).colony_number = SampleNamesSimple_colnum{next_index};
    sample_metadata(i).sebaceous_skin_region = zone_key{zones(next_index)};
    sample_metadata(i).time_in_months = num2str(times(next_index));
    sample_metadata(i).median_coverage = num2str(coverage_chromosome_median(ismember(SampleNames_all,SampleNames(next_index))));
    % read length
    sample_metadata(i).read_length = num2str(read_lengths( ismember( SampleNames_rl, next_name ) ));
    % plasmid presence
    sample_metadata(i).plasmid_presence = num2str(has_plasmid( ismember( SampleNamesSimple_plasmid, next_simple_name ) ));
    % transposon presence
    sample_metadata(next_index).transposon_presence = num2str(has_transposon( ismember( SampleNamesSimple_plasmid, next_simple_name ) ));
end

% build CSV
fid = fopen( [ dir_save '/' 'temp_supp_table_99_everything.csv' ], 'w' );
% lineage name
column_buffer_0 = char;
for i=1:snv_metadata_count+numel(unclustered_all)
    column_buffer_0(end+1) = ',';
end
fprintf( fid, [ 'unclustered_colonies,' column_buffer_0 '\n'] );
% sample metadata
column_buffer_1 = char;
for i=1:snv_metadata_count
    column_buffer_1(end+1) = ',';
end
for r=2:sample_metadata_count
    fprintf( fid, column_buffer_1 );
    fprintf( fid, [ sample_metadata_list{r} ',' ] );
    for n=1:numel(unclustered_all)
        fprintf( fid, [ sample_metadata(n).(sample_metadata_list{r}) ','] );
    end
    fprintf( fid, '\n' );
end

% build CSV: new version colony metadata only

fid = fopen( [ dir_save '/' 'temp_supp_table_99_colonytable.csv' ], 'w' );
% lineage name
fprintf( fid, [ 'unclustered_colonies,' column_buffer_0 '\n'] );
% sample metadata headers
for r=2:sample_metadata_count
    fprintf( fid, [ sample_metadata_list{r} ',' ] );
end
fprintf( fid, '\n' );
% sample metadata
for n=1:numel(unclustered_all)
    for r=2:sample_metadata_count
         fprintf( fid, [ sample_metadata(n).(sample_metadata_list{r}) ','] );
    end
    fprintf( fid, '\n' );
end
fclose(fid);


%% Put everything into one big csv file

cd( dir_save )
! echo "SNV table and metadata for all lineages" > supp_table_all-lineages_everything.csv
! echo "" >> supp_table_all-lineages_everything.csv
! for f in temp_supp_table_*_everything.csv; do cat $f; echo; done >> supp_table_all-lineages_everything.csv
cd( '..' )

cd( dir_save )
! echo "SNV table for all lineages" > supp_table_all-lineages_snvtable.csv
! echo "" >> supp_table_all-lineages_snvtable.csv
! for f in temp_supp_table_*_snvtable.csv; do cat $f; echo; done >> supp_table_all-lineages_snvtable.csv
cd( '..' )

cd( dir_save )
! echo "Colony metadata for all lineages" > supp_table_all-lineages_colonytable.csv
! echo "" >> supp_table_all-lineages_colonytable.csv
! for f in temp_supp_table_*_colonytable.csv; do cat $f; echo; done >> supp_table_all-lineages_colonytable.csv
cd( '..' )