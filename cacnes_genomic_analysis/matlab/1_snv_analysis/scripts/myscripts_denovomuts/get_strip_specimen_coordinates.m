function [ spec_coor_specnums, spec_coor_stripnum, spec_coor_x, spec_coor_y ] = ...
    get_strip_specimen_coordinates

% Get data on how many pores each pore specimen represents

%% Versions

% Arolyn, 2020.10.31: tdfreadunix to tdfread for compatibility with 2019b

%% Execution

%% Get coordinates

fprintf(1,'Warning! Double-check that both A and F strips are in millimeters!\n')

% Get specimen multiplicity (# of pores) for Subject A strips
StripCoordinateFileA = 'data_sampleinfo/strip_info/specimen_log_strips_A.csv';
csvAllA=tdfread(StripCoordinateFileA,','); % structure with field names
csvSpecimenNumberA = csvAllA.SpecimenNumber; % specimen number column
csvStripSegmentA = csvAllA.StripNumber; % strip number column
csvStripCoorXA = csvAllA.PosX; % strip x pos column
csvStripCoorYA = csvAllA.PosY; % strip y pos column

% Get specimen multiplicity (# of pores) for Subject F strips
StripCoordinateFileF = 'data_sampleinfo/strip_info/specimen_log_strips_F.csv';
csvAllF=tdfread(StripCoordinateFileF,','); % structure with field names
csvSpecimenNumberF = csvAllF.SpecimenNumber; % specimen number column
csvStripSegmentF = csvAllF.StripNumber; % strip number column
csvStripCoorXF = csvAllF.PosX; % strip x pos column
csvStripCoorYF = csvAllF.PosY; % strip y pos column


%% Save

% Save everything
spec_coor_specnums = [ csvSpecimenNumberA; csvSpecimenNumberF ];
spec_coor_stripnum = [ csvStripSegmentA; csvStripSegmentF ];
spec_coor_x = [ csvStripCoorXA; csvStripCoorXF ];
spec_coor_y = [ csvStripCoorYA; csvStripCoorYF ];

save('data_sampleinfo/strip_info/strip_coordinates',...
    'spec_coor_specnums','spec_coor_stripnum','spec_coor_x','spec_coor_y')

end