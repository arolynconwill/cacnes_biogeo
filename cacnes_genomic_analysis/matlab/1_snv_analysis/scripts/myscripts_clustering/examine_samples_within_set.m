function [ mySetSamplesToRemoveIndexedAll, mySetPositionsImportantNum ] = ...
    examine_samples_within_set( this_cluster, clusters_all, this_subject, ...
    SampleNames_all, Calls_all, coverage_all, maf_all, Quals_all, distance_matrix, p_all, ...
    figures_boolean, min_pos_to_examine_cluster, filtering_mean_maf_per_sample, filtering_min_cov_to_consider_pos, ...
    dirFigName )

%% Summary

% This function examines samples within one cluster, and filters out
% samples suspected of being contaminated (i.e. have mixed alleles at
% positions that differentiate samples within the cluster).

% There is also an option (figures_boolean) to save figures summarizing
% which positions on the genome were identified as candidate SNPs and which
% samples were identified as potentially contaminated (and why). 

% Inputs:
% this_cluster; % which cluster to look at
% clusters_all; % clusters % passed through
% SampleNames_all; % sample names % passed through
% Calls_all; % calls % passed through
% maf_all; % major allele frequencies % passed through
% coverage_all; % coverage % passed through
% distance_matrix; % distance matrix % passed through
% p_all; % positions % passed through
% figures_boolean; % whether or not to make figure files 
% min_pos_to_examine_cluster; % only perform filtering samples if sufficient candidate SNPs were identified % passed through
% filtering_mean_maf_per_sample; % minimum allowed mean major allele frequency over candidate SNPs for a sample to pass filtering  % passed through
% filtering_min_cov_to_consider_pos; % minimum coverage for a position to be included in sample filtering  % passed through
% dirFigName; % name of directory to store figures

% Outputs:
% samples_to_remove = indices of samples not to include in clustering
% mySetPositionsImportantNum = number of important positions considered


%% Collect data corresponding to this cluster only

% Name for this set of samples from this subject
mySamplesetName = [ 'Cluster-' num2str(this_cluster) '-' this_subject ]; % name for labeling plots
fprintf(1,['(Filtering) Currently examining: ' mySamplesetName '\n']); % print current set to console

% Names and indices for this set
mySetIndicesNum = clusters_all{ this_cluster }; % indices of samples in *_all
mySetIndicesBinary = ismember( linspace(1,numel(SampleNames_all),numel(SampleNames_all)), mySetIndicesNum );
mySetSampleNames = SampleNames_all( mySetIndicesBinary ); % names of samples

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

% Calls for this sample set
% Take calls used in distance matrix for this sample set only
mySetCallsDist = Calls_all( :, mySetIndicesBinary );

% Coverage for this sample set
mySetCallsCoverage = transpose( coverage_all( :,mySetIndicesBinary ) );

% Major allele frequencies for this sample set
mySetCallsMAFs = maf_all( :, mySetIndicesBinary);

% Quals for this sample set
mySetQuals = Quals_all( :, mySetIndicesBinary);


%% Identify positions of interest for this set of samples

% Find a set of positions that are likely to be important in distinguishing
% samples from this set from each other
[ mySetPositionsImportant, ~ ] = find_sampleset_important_positions( mySetCallsCoverage, mySetCallsDist, p_all );
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
cluster_filtering_directory_name = dirFigName;
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
    plot_sampleset_coverage( mySetCallsCoverage, mySetPositionsImportant, p_all, mySetSampleNames, mySampleSortOrder, mySamplesetName, cluster_filtering_directory_name )
    % Heatmap of coverage of a subset of positions on the genome that 
    % might be important for inferring relationships within this set of 
    % samples; coverage currently normalized over each position (but this
    % could change)

    % Plot 3: MAFs within Sample Set
    plot_sampleset_mafs( mySetCallsMAFs, mySetPositionsImportant, p_all, mySetSampleNames, mySampleSortOrder, mySamplesetName, cluster_filtering_directory_name )
    % Heatmap of MAFs of a subset of positions on the genome
    
    % Plot 4: MAFs within Sample Set
    plot_sampleset_quals( mySetQuals, mySetPositionsImportant, p_all, mySetSampleNames, mySampleSortOrder, mySamplesetName, cluster_filtering_directory_name )
    % Heatmap of MAFs of a subset of positions on the genome

    % Plot 5: Calls within Sample Set
    plot_sampleset_calls( mySetCallsDist, mySetPositionsImportant, p_all, mySetSampleNames, mySampleSortOrder, mySamplesetName, cluster_filtering_directory_name )
    % Heatmap of calls at a subset of positions on the genome that might be
    % important for inferring relationships within this set of samples;
    % four colors = four bases; gray = N
    
end


%% Filter samples within this set

% Filter out samples
if mySetPositionsImportantNum>0
    mySetSamplesToRemove = filter_samples_within_set( ...
        mySetCallsCoverage, mySetCallsMAFs, mySetQuals, mySetCallsDist, ...
        mySetPositionsImportant, mySamplesetName, mySetSampleNames, mySampleSortOrder, ...
        figures_boolean, cluster_filtering_directory_name, ...
        filtering_mean_maf_per_sample, filtering_min_cov_to_consider_pos, ...
        min_pos_to_examine_cluster );
else
    mySetSamplesToRemove = [];
    fprintf(1,['Skipping filtering!\n'])
end
% Note: There are parameters inside filter_samples_within_set that you
% may want to modify, depending on your sample set!!!

% Reindex samples to remove so that it corresponds to _all indexing in
% identify_clusters: 
mySetSamplesToRemoveIndexedAll = mySetIndicesNum(mySetSamplesToRemove);

% Print out list of names of samples to be removed
% But only if sufficient positions were considered
if mySetPositionsImportantNum >= min_pos_to_examine_cluster
    fprintf(1,[ 'Removing the following ' num2str(numel(mySetSamplesToRemoveIndexedAll)) ' samples:\n' ])
    for thisSample=1:length(mySetSamplesToRemoveIndexedAll)
        fprintf(1,SampleNames_all{mySetSamplesToRemoveIndexedAll(thisSample)})
        fprintf(1,'\n')
    end
end


end
