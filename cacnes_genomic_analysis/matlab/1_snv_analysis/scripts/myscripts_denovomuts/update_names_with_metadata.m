function SampleNamesLong_all = ...
    update_names_with_metadata( SampleNames_all, slst_all, clusters_all, unclustered_all )

% Update SampleNames_Long
SampleNamesFormattedSLST = {};
for i=1:numel(SampleNames_all)
    SampleNamesFormattedSLST{end+1} = [ SampleNames_all{i} '_SLST-' slst_all{i} ];
end
SampleNamesLong_all = SampleNamesFormattedSLST;

% Figure out which cluster each sample is in
SampleClusters_final = -ones(numel(SampleNamesLong_all),1);
% New tree names that include cluster assignment
SampleNamesLong_all_withClade = {};
% Loop through all samples and check if they are in a cluster
for nextNameIndex=1:numel(SampleNamesLong_all)
    nextName = SampleNamesLong_all{nextNameIndex}; % name of the next sample
    nextNameClusterFind = cell2mat(cellfun(@(x) sum(ismember(x,nextNameIndex)), clusters_all, 'UniformOutput', false)); % boolean membership of each cluster
    if sum(nextNameClusterFind)>0 % if sample is in a cluster
        nextClusterNum = find(nextNameClusterFind); % cluster number
        SampleClusters_final(nextNameIndex) = nextClusterNum; % save cluster number
        % Update name to indicate cluster number
        if nextClusterNum < 10
            SampleNamesLong_all_withClade{end+1} = [ nextName '_C-0' num2str(nextClusterNum) ];
        else
            SampleNamesLong_all_withClade{end+1} = [ nextName '_C-' num2str(nextClusterNum) ];
        end
    elseif ismember(nextNameIndex,unclustered_all)
        % Update name to indicate this sample is not in a cluster
        SampleNamesLong_all_withClade{end+1} = [ nextName '_C-00' ];
        SampleClusters_final(nextNameIndex) = 0;
    else
        % Update name to indicate this sample was filtered
        SampleNamesLong_all_withClade{end+1} = [ nextName '_filtered' ];
    end
end

% Check naming
for i=1:numel(SampleNamesLong_all_withClade); fprintf(1,[SampleNamesLong_all_withClade{i} '\n']); end

SampleNamesLong_all = SampleNamesLong_all_withClade;