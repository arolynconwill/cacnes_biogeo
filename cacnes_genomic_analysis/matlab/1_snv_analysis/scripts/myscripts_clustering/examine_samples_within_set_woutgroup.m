function mySetPositionsImportantNum = ...
    examine_samples_within_set_woutgroup( this_cluster, clusters_all, outgroup, this_subject, ...
    SampleNames_all, Calls_all, coverage_all, maf_all, distance_matrix, p_all, ...
    figures_boolean, min_pos_to_examine_cluster, ...
    dirFigName )

% This function examines samples within one cluster, and filters out
% samples that do not meet filtering criteria. There is also an option to
% make figures showing calls/coverage/similarity within this set of
% samples.
% Note that this could also be used to look at an arbitrary set of samples,
% but you would need to slightly adjust the inputs and labeling.

% Inputs!
% this_cluster; % which cluster to look at
% clusters_all; % clusters % passed through
% outgroups; % outgroups % passed through
% SampleNames_all; % sample names % passed through
% Calls_all; % calls % passed through
% maf_all; % major allele frequencies % passed through
% coverage_all; % coverage % passed through
% distance_matrix; % distance matrix % passed through
% p_all; % positions % passed through
% figures_boolean = true; % whether or not to make figure files % DEFINE!
% min_pos_to_examine_cluster; % only remove samples if sufficient positions were identified in analysis
% dirFigName; % name of directory to store figures

% Note that you could modify the inputs to be more generic, i.e. instead of
% clusters_all, just take a list of indices corresponding to some set of
% interest according to the indexing in *_all.

% Outputs! 
% samples_to_remove = indices of samples not to include in clustering
% mySetPositionsImportantNum = number of important positions considered


%% Version history:
% 2018.10.04 Arolyn: Initial version
% 2018.11.14 Arolyn: New input added for the minimum number of important
% positions that must be found to examine the cluster
% 2018.11.14 Arolyn: New input added for the name of the directory for
% storing figures
% 2018.11.15 Arolyn: Now we can look at outgroup samples
% 2019.11.22 Arolyn: Save positions to examine later
% 2019.12.05 Arolyn: Remove "filtering"
% ...


%% Collect data corresponding to this subject only

% Name for this set of samples from this subject
mySamplesetName = [ 'Cluster-' num2str(this_cluster) '-' this_subject ]; % name for labeling plots
fprintf(1,['(Filtering) Currently examining: ' mySamplesetName '\n']); % print current set to console

% Names and indices for this set
mySetIndicesNum = union( clusters_all{ this_cluster }, outgroup{ this_cluster } ); % indices of samples in *_all
mySetIndicesBinary = ismember( linspace(1,numel(SampleNames_all),numel(SampleNames_all)), mySetIndicesNum );
% Also without the outgroup
mySetIndicesNumNoOutgroup = clusters_all{ this_cluster }; % indices of samples in *_all
mySetIndicesBinaryNoOutgroup = ismember( linspace(1,numel(SampleNames_all),numel(SampleNames_all)), mySetIndicesNumNoOutgroup );
%mySetIndicesBinary = zeros( numel(SampleNames_all),1 ); % initialize
%mySetIndicesBinary( mySetIndicesNum ) = 1; % 1 or 0 for in / not in cluster over *_all indices
% Add OUTGROUP tag to outgroup samples
SampleNames_all_OutgroupTag = SampleNames_all;
indicesOutgroupSamples = outgroup{ this_cluster };
numOutgroupSamples = length( indicesOutgroupSamples );
for s=1:numOutgroupSamples
    nextIndex = indicesOutgroupSamples(s);
    SampleNames_all_OutgroupTag{nextIndex} = [ 'OUTGROUP_' SampleNames_all{nextIndex} ];
end
mySetSampleNames = SampleNames_all_OutgroupTag( mySetIndicesBinary ); % names of samples

% Distance matrix and similarity matrix
% Grab subset of big distance matrix corresponding to this sample set
myMatrixDistance = distance_matrix( mySetIndicesNum,mySetIndicesNum );
% Calculate similarity matrix
myMatrixSimilarity = max(myMatrixDistance(:)) - myMatrixDistance;

% Order samples by average similarity to other samples
% Average similarity to other samples for each sample
mySampleSimilarityAverages = mean( myMatrixSimilarity ); % average over rows of matrix
% Get order based on sample average similarity
[ ~, mySampleSortOrder ] = sort( mySampleSimilarityAverages ); % by index
% Presumably this will still work with an outgroup, so these should be the
% least similar samples to everything else (one would hope...)

% Calls for this sample set
% Take calls used in distance matrix for this sample set only
mySetCallsDist = Calls_all( :, mySetIndicesBinary );
mySetCallsDist( (p_all>=1793873 & p_all<=1803257) | (p_all>=1411528 & p_all<=1416658), : ) = 0;
% Also without the outgroup
mySetCallsDistNoOutgroup = Calls_all( :, mySetIndicesBinaryNoOutgroup );
mySetCallsDistNoOutgroup( (p_all>=1793873 & p_all<=1803257) | (p_all>=1411528 & p_all<=1416658), : ) = 0;
% This step removes the same positions that were removed in
% identify_clusters; the rationale behind not passing through Calls_dist is
% that Calls_dist is only created in identify_clusters when there is not
% already a distance matrix.

% Coverage for this sample set
mySetCallsCoverage = transpose( coverage_all( :,mySetIndicesBinary ) );
% Also without the outgroup
mySetCallsCoverageNoOutgroup = transpose( coverage_all( :,mySetIndicesBinaryNoOutgroup ) );

% Major allele frequencies for this sample set
mySetCallsMAFs = maf_all( :, mySetIndicesBinary);


%% Identify positions of interest for this set of samples

% Find a set of positions that are likely to be important in distinguishing
% samples from this set from each other
% DO NOT USE OUTGROUP SAMPLES HERE!!!
[ mySetPositionsImportant, ~ ] = find_sampleset_important_positions( mySetCallsCoverageNoOutgroup, mySetCallsDistNoOutgroup, p_all );
% Note: There are parameters inside find_sampleset_important_positions that
% you may need to change!!!

% Number of interesting positions found
mySetPositionsImportantNum = numel(mySetPositionsImportant);

% Return a warning if there are no important positions found
if mySetPositionsImportantNum < 1
    fprintf(1,'Warning! No interesting positions found!!\n')
end


%% Make plots of sample similarity, coverage, and calls

% Make a special directory for all the figures
cluster_filtering_directory_name = [ dirFigName ];
if ~exist(cluster_filtering_directory_name,'dir')
    mkdir(cluster_filtering_directory_name)
end

if figures_boolean && mySetPositionsImportantNum>0 % only make plots if requested and if there are actually interesting positions to plot

    % Save positions
    save([cluster_filtering_directory_name '/' 'Data_' mySamplesetName '_Positions'],'mySetPositionsImportant')
        
    % Plot 1: Similarity with Sample Set
    plot_sampleset_similarities( myMatrixSimilarity, mySetSampleNames, mySampleSortOrder, mySamplesetName, cluster_filtering_directory_name )
    % Heatmap of histograms of similarity of each sample to all other
    % samples within the set
    
    % Plot 2: Coverage within Sample Set
    plot_sampleset_coverage_morepos( mySetCallsCoverage, mySetPositionsImportant, p_all, mySetSampleNames, mySampleSortOrder, mySamplesetName, cluster_filtering_directory_name )
    % Heatmap of coverage of a subset of positions on the genome that 
    % might be important for inferring relationships within this set of 
    % samples; coverage currently normalized over each position (but this
    % could change)

    % Plot 3: MAFs within Sample Set
    plot_sampleset_mafs_morepos( mySetCallsMAFs, mySetPositionsImportant, p_all, mySetSampleNames, mySampleSortOrder, mySamplesetName, cluster_filtering_directory_name )
    % Heatmap of MAFs of a subset of positions on the genome
    
    % Plot 4: Calls within Sample Set
    plot_sampleset_calls_morepos( mySetCallsDist, mySetPositionsImportant, p_all, mySetSampleNames, mySampleSortOrder, mySamplesetName, cluster_filtering_directory_name )
    % Heatmap of calls at a subset of positions on the genome that might be
    % important for inferring relationships within this set of samples;
    % four colors = four bases; gray = N
    
end


