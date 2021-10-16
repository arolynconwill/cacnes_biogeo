function [ mySamplesetImportantPositionsIndices, mySamplesetImportantPositionsGenome ] ...
    = find_sampleset_important_positions( mySamplesetCoverage, mySamplesetCalls, p_all )

%% Summary

% This function takes calls and coverage data for a set of samples and
% returns a set of positions that are likely to be important in
% distinguishing these samples from each other. This is basically rough SNP
% calling.

% Inputs: 
% mySubjectCallsCoverage = coverage table
% mySamplesetCalls = calls table
% p_all = positions on genome

% Outputs: 
% mySamplesetImportantPositionsIndices = indices of important positions
% mySamplesetImportantPositionsGenome = positions of important positions on genome


%% Basic info on sample set

mySamplesetSize = length(mySamplesetCalls(1,:));
mySamplesetNumPositions = numel(p_all);


%% Find positions that are important for this sample set

% 1) Eliminate positions that are trivially non-diverse
% i.e. all calls (excluding N's) are the same

% First: Ignore trivially non-diverse positions, i.e. calls over all 
% samples are the same whenever they aren't N
myPositionsNonvariable = (...
    sum( mySamplesetCalls==1, 2 ) == sum( mySamplesetCalls~=0, 2) | ...
    sum( mySamplesetCalls==2, 2 ) == sum( mySamplesetCalls~=0, 2) | ...
    sum( mySamplesetCalls==3, 2 ) == sum( mySamplesetCalls~=0, 2) | ...
    sum( mySamplesetCalls==4, 2 ) == sum( mySamplesetCalls~=0, 2) ...
    );

% Second: Ignore broadly ambiguous positions, i.e. calls are N over almost
% all of the samples
max_fraction_Ns_across_samples = 0.67;
myPositionsAmbiguous = (...
    sum( mySamplesetCalls==0, 2 ) >= mySamplesetSize*max_fraction_Ns_across_samples ...
    );

% Third: Ignore positions that generally have low coverage across samples
min_median_coverage_across_samples = 10;
myPositionsLowCoverage = (...
    median( transpose(mySamplesetCoverage), 2) < min_median_coverage_across_samples ...
    );

% Combine all criteria to yield a list of positions that are good enough
% quality to consider here
myPositionsToConsider = (...
    ~myPositionsNonvariable & ...
    ~myPositionsAmbiguous & ...
    ~myPositionsLowCoverage ...
    );

% Now look at positions where there is variability in base calls
% "MAF" across samples
mySamplesetCallsMAFs = -ones( length(mySamplesetCalls),1 ); % initialize
for thisPos=1:mySamplesetNumPositions
    % Major allele at this position
    thisPositionMajorAllele = mode(mySamplesetCalls(thisPos,:)); % most common call, could be N
    % Major allele frequency at this position
    if thisPositionMajorAllele == 0 
        mySamplesetCallsMAFs(thisPos) = 2; % to indicate major allele over all samples was N
    else
        mySamplesetCallsMAFs(thisPos) = sum( mySamplesetCalls(thisPos,:) == thisPositionMajorAllele ) ...
        / sum( mySamplesetCalls(thisPos,:) ~= 0 ); % save MAF taken over all samples that weren't N
    end
end
max_maf_across_sample_set = 0.999; % currently allows positions with rare SNPs, but can change to focus on SNPs that are prevalent within the sample set
myPositionsPossibleMutations = mySamplesetCallsMAFs < max_maf_across_sample_set;

% Important positions = positions to consider that also have possible
% mutations
mySamplesetImportantPositionsBooleans = ( ...
    myPositionsToConsider & ...
    myPositionsPossibleMutations ...
    );
% NOTE: Further optimization could happen here...
fprintf(1,['...considering ' num2str(sum(mySamplesetImportantPositionsBooleans)) ' positions...\n'])

% Find the indices of these important positions
mySamplesetImportantPositionsIndices = find( mySamplesetImportantPositionsBooleans ); % indices for *_all
mySamplesetImportantPositionsGenome = p_all( mySamplesetImportantPositionsIndices ); % positions on genome


end
