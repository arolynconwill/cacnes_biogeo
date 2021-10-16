function plot_sampleset_coverage_v2_morepos( mySamplesetCoverage, p_important, p_all, plotSampleNames, plotSampleOrder, mySamplesetName, cluster_filtering_directory_name )

% This function makes a heatmap of sample coverage across interesting
% positions for a set of samples; normalized on each position (but this
% could change...); note that this only plots a subset of important 
% positions (or else the figure would be too wide)

% Inputs: 
% mySamplesetCoverage = coverage table
% p_important = subset of important positions
% p_all = positions on genome
% plotSampleNames = names of samples in order of coverage table
% plotSampleOrder = order in which to plot samples
% mySamplesetName = name for plot label

% Outputs: 
% None! (just writes a file)


%% Version history:
% 2018.10.04 Arolyn: Initial version
% 2018.11.14 Arolyn: Longer figure saved if over 100 samples
% 2018.11.14 Arolyn: Change to absolute coverage, not normalized coverage
% ...


%% Get coverage data for all samples at all positions

% Number of samples
mySubjectNumSamples = numel(plotSampleNames);

% % Coverage normalized by average coverage at that position
% mySamplesetCoverageNormalized = -1*ones( size(mySamplesetCoverage) ); % initialize
% % Mean coverage at each position
% mySamplesetCoverageMeansPosition = mean(mySamplesetCoverage,1);
% % Normalized coverage at each position for each sample
% for thisSample=1:mySubjectNumSamples
%     mySamplesetCoverageNormalized(thisSample,:) = mySamplesetCoverage(thisSample,:)./mySamplesetCoverageMeansPosition;
% end
% 
% % Normalized coverage at important positions only
% mySamplesetCoverageNormalizedImportantPositions = mySamplesetCoverageNormalized(:,p_important);


%% CHANGE: Try plotting coverage without normalization, but rail at 100

mySamplesetCoverageImportantPositions = mySamplesetCoverage(:,p_important);
coverage_rail = 50; % max coverage to show on plot
% Set anything higher than coverage_rail to coverage_rail
mySamplesetCoverageImportantPositions( mySamplesetCoverageImportantPositions > coverage_rail ) = coverage_rail;


%% Make a figure

% Subset of positions to plot, sampled across all important positions
pos_min = 1;
pos_max = numel(p_important);
if numel(p_important) > 500 % High so that everything gets plotted
    pos_skip = round(pos_max/100);
else
    pos_skip = 1;
end
% Labels for horizontal and vertical axes
xvalues = p_all( p_important( pos_min:pos_skip:pos_max ) ); % positions on genome
yvalues = plotSampleNames( plotSampleOrder ); % sample names

% Make figure
figure(1)
clf(1,'reset')
%figCov=heatmap(xvalues,
%yvalues,mySamplesetCoverageNormalizedImportantPositions(plotSampleOrder,pos_min:pos_skip:pos_max)); % normalized
figCov=heatmap(xvalues, yvalues,mySamplesetCoverageImportantPositions(plotSampleOrder,pos_min:pos_skip:pos_max)); % not normalized
% Color scheme
%myColormap = flipud(hot(20)); % old
myColormap = hot(30);
myColormap = myColormap(1:25,:); % avoid white
figCov.Colormap = myColormap;
caxis([ 0 coverage_rail ]); % force colormap over whole range of values
% Figure title and axis labels
figCov.Title = [ 'Coverage at Important Positions: ' mySamplesetName ];
figCov.XLabel = ['Positions on Genome'];
figCov.YLabel = 'Sample Names';
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
print([pwd '/' cluster_filtering_directory_name '/' 'Figure_' mySamplesetName '_Coverage.png'],'-dpng')

