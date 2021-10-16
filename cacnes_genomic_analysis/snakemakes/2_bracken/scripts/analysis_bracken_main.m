%% Bracken analysis 

% 1) Record the number and % of reads assigned to C. acnes
% 2) Identify colonies that are actually C. granulosum


%% Setup

% Make folder in which to save data
if ~exist( [pwd '/3-matlab'], 'dir')
    mkdir( [pwd '/3-matlab'] )
end


%% Get sample names

% Where to fine bracken files (downloaded from cluster)
path_to_sample_data = '2-kraken2';

% Get list of all files, one per sample
list_of_files = dir( path_to_sample_data ); % get table of stuff in directory
list_of_files_names = {list_of_files.name}; % get file names
list_of_files_bracken = cell2mat( cellfun(@(x) contains(x,'.bracken'), list_of_files_names, 'UniformOutput', false) );
list_of_files_names = list_of_files_names(list_of_files_bracken); % remove . and .. and .DS_store etc.

% Get sample names from file names 
SampleNames_bracken = cellfun(@(x) x(1:end-8), list_of_files_names, 'UniformOutput', false);
for i=1:numel(SampleNames_bracken)
    next_name = SampleNames_bracken{i};
    if next_name(2) == '-'
        next_name(2) = ':'; % Fixing because of weird issue where filenames changed when downloaded from c3ddb
    end
    SampleNames_bracken{i} = next_name;
end


%% Get taxa of interest
% Get taxa >1% of reads in at least one sample according to bracken

% Initialize list of taxa to track
list_taxa = {}; % cell
list_taxa_ids = []; % nums

% Minimum frac to pull from file
min_abundance = 0.0001;

for i=1:numel( list_of_files_names )
    fprintf(1,[ 'Reading file ' num2str(i) '/' num2str(numel(list_of_files_names)) '.' '\n' ])
    next_filename = [ path_to_sample_data '/' list_of_files_names{i} ];
    [ next_taxa, next_taxa_ids, ~, ~ ] = ...
        read_bracken_output( next_filename, min_abundance );
    % Identify taxa not already in list
    taxa_to_add = ~ismember( next_taxa_ids, list_taxa_ids ); % taxa not already in list
    if sum(taxa_to_add)>0 % if there are taxa to add
        % Add new IDs
        taxa_ids_to_add = next_taxa_ids( taxa_to_add);
        list_taxa_ids = [ list_taxa_ids; taxa_ids_to_add ];
        % Add new names
        taxa_names_to_add = next_taxa( taxa_to_add,: );
        for n=1:numel(taxa_ids_to_add)
            list_taxa{end+1} = strtrim(taxa_names_to_add(n,:));
        end
    end
end


%% Record %s of taxa of interest
% Also save % and # of C acnes reads

% Initialize matrices of fractions and reads
fracs_mat = zeros( numel(list_of_files_names), numel(list_taxa_ids) );
reads_mat = zeros( numel(list_of_files_names), numel(list_taxa_ids) );

% Read files again and pull data
for i=1:numel( list_of_files_names )
    fprintf(1,[ 'Reading file ' num2str(i) '/' num2str(numel(list_of_files_names)) '.' '\n' ])
    % Get data
    next_filename = [ path_to_sample_data '/' list_of_files_names{i} ];
    [ next_taxa, next_taxa_ids, next_fracs, next_reads ] = ...
        read_bracken_output( next_filename, min_abundance );
    % Find indices
    [ ~, next_indices ] = ismember( next_taxa_ids, list_taxa_ids );
    % Record data
    fracs_mat(i,next_indices) = next_fracs;
    reads_mat(i,next_indices) = next_reads;
end

% Save!
save('3-matlab/data_bracken_all', 'SampleNames_bracken', 'list_taxa_ids', 'fracs_mat', 'reads_mat' )


%% Save data for filtering colonies based on C acnes identity metrics
% C acnes percent
% C acnes number of reads

% Find C. acnes fractions and reads
index_cacnes = find( ismember(list_taxa,'Cutibacterium acnes') );
cacnes_fracs = fracs_mat( :,index_cacnes );
cacnes_reads = reads_mat( :,index_cacnes );

% Save!
save('3-matlab/data_bracken_cacnes', 'SampleNames_bracken', 'cacnes_fracs', 'cacnes_reads' )


%% IDENTIFY C GRANULOSUM COLONIES

% Find samples that are C. granulosum

identity_cutoff = 0.75; % fraction C granulosum must exceed this number
% obviously this is generous; expecting some of these to get filtered later

index_cgranulosum = find( ismember(list_taxa,'Cutibacterium granulosum' ) );
bracken_fracs_cgranulosum = fracs_mat( :,index_cgranulosum );
bracken_reads_cgranulosum = reads_mat( :,index_cgranulosum );
is_cgranulosum = ( bracken_fracs_cgranulosum > identity_cutoff );

bracken_fracs_cgranulosum = bracken_fracs_cgranulosum( is_cgranulosum );
bracken_reads_cgranulosum = bracken_reads_cgranulosum( is_cgranulosum );
bracken_samplenames_cgranulosum = SampleNames_bracken( is_cgranulosum );
numel( bracken_samplenames_cgranulosum )

% Save C granulosum info

save( '3-matlab/data_bracken_cgran', ...
    'bracken_samplenames_cgranulosum', ...
    'bracken_fracs_cgranulosum', 'bracken_reads_cgranulosum' )

% Make samples.csv

csv_table = readtable( 'samples.csv', 'Delimiter', ',' );
csv_table_names = csv_table.Sample;

fid = fopen( '3-matlab/samples_cgran.csv', 'w' );
fprintf(fid, 'Path,Sample,ReferenceGenome,ProviderName,Subject \n' );
for i=1:numel(bracken_samplenames_cgranulosum)
    next_name = bracken_samplenames_cgranulosum{i};
    next_index = find( ismember( csv_table_names, next_name ) );
    fprintf(fid, ...
        [ csv_table.Path{next_index} ',' ...
        next_name ',' ...
        'C_GRANULOSUM_REF_GENOME_NAME' ',' ...
        csv_table.ProviderName{next_index} ',' ...
        csv_table.Subject{next_index} ' \n' ] ...
        );
end
fclose(fid);
