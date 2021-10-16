function [ spec_mult_specnums, spec_mult_porenums ] = ...
    get_pore_specimen_multiplicity( specimen_number_all, SampleNames_all, types_all )

% Get data on how many pores each pore specimen represents

%% Versions

% Arolyn, 2020.10.31: tdfreadunix to tdfread for compatibility with 2019b

%% Execution

%% Strips

% Get specimen multiplicity (# of pores) for Subject A strips
StripCoordinateFileA = 'data_sampleinfo/strip_info/specimen_log_strips_A.csv';
csvAllA=tdfread(StripCoordinateFileA,','); % structure with field names
csvSpecimenNumberA = csvAllA.SpecimenNumber; % specimen number column
csvMultA = csvAllA.Multiple; % strip number column

% Get specimen multiplicity (# of pores) for Subject F strips
StripCoordinateFileF = 'data_sampleinfo/strip_info/specimen_log_strips_F.csv';
csvAllF=tdfread(StripCoordinateFileF,','); % structure with field names
csvSpecimenNumberF = csvAllF.SpecimenNumber; % specimen number column
csvMultF = csvAllF.Multiple; % strip number column


%% Extracts

% Get specimen multiplicity (# of pores) for extracts
extract_spec = specimen_number_all( types_all == 1 );
extract_names = SampleNames_all( types_all == 1 );
% Parse extract tags
extract_tag = cellfun(@(x) strtok( x,'_' ), extract_names, 'UniformOutput', false); % remove time
extract_tag = cellfun(@(x) x(1:end-2), extract_tag, 'UniformOutput', false); % remove colony number
extract_tag = cellfun(@(x) fliplr(strtok(fliplr(x),'-')), extract_tag, 'UniformOutput', false)'; % removes subject and specimen number
extract_single = cellfun(@(x) ~contains(x,'Bs'), extract_tag) & cellfun(@(x) ~contains(x,'Ws'), extract_tag);
extract_single( ismember( extract_spec, [ 12 37 87 406 407 ] ) ) = 1; % see notes from Tami
extract_mult = extract_single + 2.5*~extract_single; % 1 for singles and 2.5 for multiples
% Only keep unique specimen numbers
[ extract_spec_unique, extract_spec_unique_indices ] = unique( extract_spec );
extract_mult_unique = extract_mult( extract_spec_unique_indices );
% % Manual updates with information in sampling notes that isn't reflected
% in sample names % used previously before line 37 added with info from
% Tami
% extract_mult_unique( extract_spec_unique==395 ) = 2.5; % "contaminated"
% extract_mult_unique( extract_spec_unique==406 ) = 2.5; % "maybe mixed"
% extract_mult_unique( extract_spec_unique==407 ) = 2.5; % "maybe mixed"


%% Save

% Save everything
spec_mult_specnums = [ csvSpecimenNumberA; csvSpecimenNumberF; extract_spec_unique ];
spec_mult_porenums = [ csvMultA; csvMultF; extract_mult_unique ];

save('data_sampleinfo/specimen_poremult/spec_mult','spec_mult_specnums','spec_mult_porenums')


%% Write CSV

fid=fopen( 'data_sampleinfo/specimen_poremult/spec_mult.csv', 'w' );
fprintf(fid, 'specimen_type,specimen_number,multiplicity,\n');
for i=1:numel(csvSpecimenNumberA)
    fprintf( fid, ['strip,' num2str(csvSpecimenNumberA(i)) ',' num2str(csvMultA(i)) ',\n'] );
end
for i=1:numel(csvSpecimenNumberF)
    fprintf( fid, ['strip,' num2str(csvSpecimenNumberF(i)) ',' num2str(csvMultF(i)) ',\n'] );
end
for i=1:numel(extract_spec_unique)
    fprintf( fid, ['extract,' num2str(extract_spec_unique(i)) ',' num2str(extract_mult_unique(i)) ',\n'] );
end
fclose(fid);

fid=fopen( 'data_sampleinfo/specimen_poremult/spec_mult_extracts.csv', 'w' );
fprintf(fid, 'specimen_type,specimen_number,multiplicity,\n');
for i=1:numel(extract_spec_unique)
    fprintf( fid, ['extract,' num2str(extract_spec_unique(i)) ',' num2str(extract_mult_unique(i)) ',\n'] );
end
fclose(fid);


end