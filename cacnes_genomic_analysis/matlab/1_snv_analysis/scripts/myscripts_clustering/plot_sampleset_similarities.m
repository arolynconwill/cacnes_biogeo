function plot_sampleset_similarities_v2( myMatrixSimilarity, plotSampleNames, plotSampleOrder, mySamplesetName, cluster_filtering_directory_name )

% This function makes a heatmap of sample similarities relative to other
% samples within a set of samples

% Inputs: 
% myMatrixSimilarity = similarity matrix
% plotSampleNames = names of samples in order of similarity matrix
% plotSampleOrder = order to plot samples 
% mySamplesetName = name of sample set for labeling the plot

% Outputs: 
% None! (just writes a file)


%% Version history:
% 2018.10.04 Arolyn: Initial version
% 2018.11.14 Arolyn: Longer figure saved if over 100 samples
% ...


%% Get histogram of sample similarities relative to all other samples in this set

% Number of samples
mySubjectNumSamples = numel(plotSampleNames);

% Define bin edges for similarity histograms
myHeatmapDataSimilarityEdgeMin = 0;
if max(max(myMatrixSimilarity)) > 1000
    myHeatmapDataSimilarityEdgeMax = round( max(max(myMatrixSimilarity))/1000 + 0.5 )*1000; % round up to nearest 1000
else
    myHeatmapDataSimilarityEdgeMax = round( max(max(myMatrixSimilarity))/100 + 0.5 )*100; % round up to nearest 100
end
if myHeatmapDataSimilarityEdgeMax > 2000
	myHeatmapDataSimilarityEdgeJump = round( myHeatmapDataSimilarityEdgeMax/50/100 )*100; % round ~50 jumps to nearest 100
else
	myHeatmapDataSimilarityEdgeJump = round( myHeatmapDataSimilarityEdgeMax/50/1 )*1; % round ~50 jumps to nearest 1
end
if 	myHeatmapDataSimilarityEdgeJump == 0
    myHeatmapDataSimilarityEdgeJump = 1; % set to one in the case that it rounds down to zero
end
% myHeatmapDataSimilarityEdgeJump = max( 1, max( round( myHeatmapDataSimilarityEdgeMax/50/5 )*5, round( myHeatmapDataSimilarityEdgeMax/50/100 )*100 ) ); % bin width for ~50 bins
myHeatmapDataSimilarityEdges = myHeatmapDataSimilarityEdgeMin:myHeatmapDataSimilarityEdgeJump:myHeatmapDataSimilarityEdgeMax;

% Histogram similarity to other samples for each sample
myHeatmapDataSimilarity = zeros( mySubjectNumSamples, numel(myHeatmapDataSimilarityEdges)-1 );
for i=1:length(myMatrixSimilarity) % loop through samples
    myHeatmapDataSimilarity(i,:) = histcounts( myMatrixSimilarity(plotSampleOrder(i),:), myHeatmapDataSimilarityEdges );
end 


%% Make figure: Similarity histograms heatmap

% Labels for horizontal axis
xvalues_start = myHeatmapDataSimilarityEdges(2)/2+myHeatmapDataSimilarityEdges(1)/2;
xvalues_jump = myHeatmapDataSimilarityEdges(2)-myHeatmapDataSimilarityEdges(1);
xvalues_end = myHeatmapDataSimilarityEdges(end)/2+myHeatmapDataSimilarityEdges(end-1)/2;
xvalues = xvalues_start:xvalues_jump:xvalues_end;
% Labels for vertical axis
yvalues = plotSampleNames(plotSampleOrder);

% Make figure
figure(1) 
clf(1,'reset')
figSim=heatmap( xvalues, yvalues, myHeatmapDataSimilarity );
figSim.Title = [ 'Similarity Histograms: ' mySamplesetName ];
figSim.XLabel = 'Similarity to Other Samples in Set';
figSim.YLabel = 'Sample Name';
set(gca, 'FontSize', 12, 'FontName', 'Verdana')

% Save figure
if length(plotSampleNames) > 100 % bigger plot if many samples
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [16 24]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 16 24]);
else
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [16 12]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 16 12]);
end
print([pwd '/' cluster_filtering_directory_name '/' 'Figure_' mySamplesetName '_Similarity.png'],'-dpng')

