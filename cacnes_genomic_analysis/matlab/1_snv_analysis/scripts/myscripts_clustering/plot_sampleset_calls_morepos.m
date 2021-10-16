function plot_sampleset_calls_v2_morepos( mySamplesetCalls, p_important, p_all, plotSampleNames, plotSampleOrder, mySamplesetName, cluster_filtering_directory_name )

% This function makes a heatmap of sample calls across interesting
% positions for a set of samples; note that this only plots a subset of
% important positions (or else the figure would be too wide)

% Inputs: 
% mySamplesetCalls = calls table
% p_important = subset of important positions
% p_all = positions on genome
% plotSampleNames = names of samples in order of calls table
% plotSampleOrder = order in which to plot samples 
% mySamplesetName = name for plot labeling

% Outputs: 
% None! (just writes a file)


%% Version history:
% 2018.10.04 Arolyn: Initial version
% 2018.11.14 Arolyn: Longer figure saved if over 100 samples
% 2018.11.14 Arolyn: Define colormap max/min so that colors still work even
% if there are no N's
% ...


%% Get calls data for all samples at important positions

% Number of samples
mySubjectNumSamples = numel(plotSampleNames);

% Calls at important positions
mySamplesetCallsImportantPositions = transpose(mySamplesetCalls(p_important,:));


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
figCalls=heatmap(xvalues,yvalues,mySamplesetCallsImportantPositions(plotSampleOrder,pos_min:pos_skip:pos_max));
% Color scheme
myColormap = parula(20); % initial colormap
myColormap = myColormap([ 2 7 13 19 ],: ); % grab colors
myColormap = [ [.925 .925 .925]; myColormap ]; % add gray
figCalls.Colormap = myColormap;
caxis([ 0 4 ]); % force colormap over whole range of values
colorbar('off')
% Figure title and axis labels
figCalls.Title = [ 'Calls at Important Positions: ' mySamplesetName ];
figCov.XLabel = ['Positions on Genome '];
figCalls.YLabel = 'Sample Names';
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
print([pwd '/' cluster_filtering_directory_name '/' 'Figure_' mySamplesetName '_Calls.png'],'-dpng')
