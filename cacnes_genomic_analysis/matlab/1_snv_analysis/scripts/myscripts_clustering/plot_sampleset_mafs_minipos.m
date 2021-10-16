function plot_sampleset_mafs_v2( mySamplesetMAFs, p_important, p_all, plotSampleNames, plotSampleOrder, mySamplesetName, cluster_filtering_directory_name )

% This function makes a heatmap of sample MAFs across interesting
% positions for a set of samples.

% Inputs: 
% mySamplesetMAFs = major allele frequency (MAF) table
% p_important = subset of important positions
% p_all = positions on genome
% plotSampleOrder = order in which to plot samples
% plotSampleNames = names of samples in order of MAF table
% mySamplesetName = name for plot label

% Outputs: 
% None! (just writes a file)


%% Version history:
% 2018.10.04 Arolyn: Initial version
% 2018.11.14 Arolyn: Longer figure saved if over 100 samples
% ...


%% Get MAF data for all samples at all positions

% Number of samples
mySubjectNumSamples = numel(plotSampleNames);

% MAFs at important positions only
mySamplesetMAFsImportantPositions = transpose(mySamplesetMAFs(p_important,:));


%% Make a figure

% Subset of positions to plot, sampled across all important positions
pos_min = 1;
pos_max = numel(p_important);
if numel(p_important) > 100
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
figCov=heatmap(xvalues, yvalues, mySamplesetMAFsImportantPositions(plotSampleOrder,pos_min:pos_skip:pos_max));
% Color scheme
%myColormap = flipud(summer(20)); % not a huge fan of this colormap... could change it but I'm lazy...
myColormap = flipud(cool(20)); % changed colormap
myColormap = flipud(myColormap(10:20,:));
figCov.Colormap = myColormap;
caxis([ 0.5 1 ]);
% Figure title and axis labels
figCov.Title = [ 'MAFs at Important Positions: ' mySamplesetName ];
figCov.XLabel = 'Positions on Genome';
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
print([pwd '/' cluster_filtering_directory_name '/' 'Figure_' mySamplesetName '_MAFs.png'],'-dpng')

