function [ SampleTimesUnitMonths, SampleNamesWithTimes ] = get_sample_times( SampleNamesWithoutTimes, SpecimenFile )

% This function finds the dates on which each sample was obtained.

% Input:
%SampleNamesWithoutTimes = SampleNames_all;
%SpecimenFile = 'Specimen_log_simple_tucker.csv';

% Output:
% % SampleTimes: Time each sample was taken, based on information in the
% specimen spreadsheet
% % SampleNamesWithTimes: Original sample names with times appended onto
% the end


%% Version History
% % 2018.11.06, Arolyn: Original version


%% Extract specimen numbers from sample names

% Reference: Set of characters that are numbers
char_digits = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

sampleSpecimenNumbers = zeros(numel(SampleNamesWithoutTimes),1); % initialize

for s=1:numel(SampleNamesWithoutTimes)
    nextName = SampleNamesWithoutTimes{s};
    % Assumes format is [Subject_Letter]:[SpecimenNumber3Digits]...
    nextSpecimen = nextName(3:5);
    % Checks if this was the format
    if ismember(nextSpecimen(1),char_digits) && ismember(nextSpecimen(2),char_digits) && ismember(nextSpecimen(3),char_digits)
    % if all three characters are digits
        nextSpecimenNum = str2num(nextSpecimen);
        sampleSpecimenNumbers(s) = nextSpecimenNum;
    else
    % if not 
        fprintf(1, ['Warning! No specimen number in this sample name: ' nextName] )
    end
end


%% Match specimen numbers in spreadsheet; extract strings with times

% Import CSV file
csvAll=tdfreadunix(SpecimenFile,','); % structure with field names

% Separate out specimen numbers and sampling times
csvSpecimenNumbers = csvAll.SpecimenNumber; % double array
csvSamplingTimes = csvAll.SamplingTime; % character array


%% Convert strings with times into actual dates

% Note: Indexed by specimen!!!

% Get years
SpecimenTimesYear = csvSamplingTimes(:,1:4); % get characters for years
SpecimenTimesYear = str2num(SpecimenTimesYear); % turn them into numbers

% Get months
SpecimenTimesMonth = csvSamplingTimes(:,6:7); % get characters for months
SpecimenTimesMonth = str2num(SpecimenTimesMonth); % turn them into numbers

% Define t=0 as the first sample ever
SampleTimesYearStart = min(SpecimenTimesYear); % 2015
SampleTimesMonthStart = min(SpecimenTimesMonth(SampleTimesYearStart==SpecimenTimesYear)); % 8

% Calculate times from t=0 in months
SpecimenTimesUnitMonths = 12*(SpecimenTimesYear-SampleTimesYearStart)+(SpecimenTimesMonth-SampleTimesMonthStart);


%% Append strings with times onto sample names and store sampling times for all samples

% Sample names with times appended
SampleNamesWithTimes = {}; % initialize
% Time of sample (months since first sample)
SampleTimesUnitMonths = -ones( numel(SampleNamesWithoutTimes),1 ); % initalize
% initialize as -1 so that it will be obvious if there is a problem...

for s=1:numel(SampleNamesWithoutTimes) % loop through all sample names
    % Find which specimen this sample belongs to
    nextSpecimenNumber = sampleSpecimenNumbers(s);
    nextSpecimenIndex = find( csvSpecimenNumbers == nextSpecimenNumber ); % index should actually be equal to the specimen number
    nextSpecimenTime = SpecimenTimesUnitMonths(nextSpecimenIndex);
    if nextSpecimenTime < 10 && nextSpecimenTime >= 0 
        nextSpecimenTimeString = [ '0' num2str(nextSpecimenTime) ];
    else
        nextSpecimenTimeString = num2str(nextSpecimenTime);
    end
    % Store sampling time
    SampleTimesUnitMonths(s) = nextSpecimenTime;
    % Append string to sample name
    SampleNamesWithTimes{end+1} = [ SampleNamesWithoutTimes{s} '_t' nextSpecimenTimeString ];
end
