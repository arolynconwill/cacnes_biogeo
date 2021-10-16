function mySubjectSamplesToRemove = filter_samples_within_set( ...
    mySubjectCallsCoverage, mySubjectCallsMAFs, mySubjectQuals, mySubjectCallsDist, ...
    p_important, mySamplesetName, plotSampleNames, plotSampleOrder, ...
    figures_boolean, cluster_filtering_directory_name, ...
    filtering_mean_maf_per_sample, filtering_min_cov_to_consider_pos, ...
    min_pos_to_examine_cluster )

%% Summary

% This function examines a set of samples (belonging to one subject) and
% determines if any of the samples should be filtered out. 

% Inputs: 
% mySubjectCallsCoverage = coverage table
% mySubjectCallsDist = calls table
% p_important = important positions to consider
% mySamplesetName = name of sample set for plot labeling
% plotSampleNames = names of samples in order of tables
% plotSampleOrder = order in which to plot samples
% filtering_mean_maf_per_sample = no more than this fraction of sketchy positions in a sample
% filtering_min_cov_to_consider_pos = minimum coverage for a position to be included in sample filtering 
% min_pos_to_examine_cluster = minimum number of positions to allow for actual filtering

% Outputs: 
% mySubjectSamplesToRemove = indices of samples to remove


%% Look at data for positions deemed important only

% Number of total samples for this subject
mySubjectNumSamples = length(mySubjectCallsCoverage(:,1));

% Coverage and calls: indexed by samples x positions
mySubjectCoveragePositionsImportant = mySubjectCallsCoverage(:,p_important); % coverage table
mySubjectMAFsPositionsImportant = transpose(mySubjectCallsMAFs(p_important,:)); % maf table
mySubjectQualsPositionsImportant = transpose(mySubjectQuals(p_important,:)); % quals table
mySubjectCallsDistPositionsImportant = transpose(mySubjectCallsDist(p_important,:)); % calls table


%% Filtering: Looks at the portion of important positions that have a base call of N despite having decent coverage

% Filtering
if numel(p_important)>=min_pos_to_examine_cluster % Not useful if too few important positions are found

    % To keep track of the fraction of sketchy positions for each sample
    mySamplesFracSketchy = -ones( mySubjectNumSamples,1 ); % initialize
    mySamplesMeanMAF = -ones( mySubjectNumSamples,1 ); % initialize 
    mySamplesMedianMAF = -ones( mySubjectNumSamples,1 ); % initialize 
    mySamples25pMAF = -ones( mySubjectNumSamples,1 ); % initialize 

    % Compute fraction of sketchy positions in each sample
    for thisSample=1:mySubjectNumSamples % loop through samples
        mySamplePositionsCalls = mySubjectCallsDistPositionsImportant(thisSample,:); % get calls
        mySamplePositionCallsN = ( mySamplePositionsCalls == 0 ); % N's
        mySamplePositionCoverage = mySubjectCoveragePositionsImportant(thisSample,:); % get coverage
        mySamplePositionCoverageOk = ( mySamplePositionCoverage > filtering_min_cov_to_consider_pos ); % enough coverage
        % Old version: Looks at fraction of N's despite good coverage
        mySamplePositionsSketchy = ( mySamplePositionCallsN & mySamplePositionCoverageOk ); % N despite enough coverage
        mySamplesFracSketchy(thisSample) = sum(mySamplePositionsSketchy) / sum(mySamplePositionCoverageOk); % fraction N despite enough coverage
        % New version: Looks at low MAF despite good coverage
        mySamplePositionsMAFs = sort( mySubjectMAFsPositionsImportant( thisSample,mySamplePositionCoverageOk ) );
        mySamplesMeanMAF(thisSample) = mean( mySamplePositionsMAFs );
        mySamplesMedianMAF(thisSample) = median( mySamplePositionsMAFs );
        index_25p = floor(numel(mySamplePositionsMAFs)/4);
        if index_25p == 0
            index_25p = 1;
        end
        mySamples25pMAF(thisSample) = mySamplePositionsMAFs( index_25p );
    end

    % New new version: Looks at quals (regardless of coverage)
    mySamplesMeanQual = mean(mySubjectQualsPositionsImportant,2);
    mySubjectSamplesToRemove = find( mySamplesMeanMAF < filtering_mean_maf_per_sample ); % Changed to MAF on 2019.12.02; fixed 2019.12.29 (previously had mySamplesMedianMAF by mistake)

else
    
    fprintf(1,['Too few positions (' num2str(numel(p_important)) ') to examine cluster...\n'])
    mySubjectSamplesToRemove = [];
    
    return % avoid plots if insufficient important positions
    
end


%% Make a plot for filtering
% Does not make a plot by default

if figures_boolean
    metric_bad = (mySamplesMeanMAF < filtering_mean_maf_per_sample);
    metric_good = (mySamplesMeanMAF >= filtering_mean_maf_per_sample);
    % Save subsets so that good/bad samples can be different colors
    mySamplesMetricGood = mySamplesMeanMAF;
    mySamplesMetricGood(metric_bad) = 0;
    mySamplesMetricBad = mySamplesMeanMAF;
    mySamplesMetricBad(metric_good) = 0;
    % Make figure
    figSketchy=figure(2);
    clf(2,'reset')
    hold on
    if numel(p_important) >= min_pos_to_examine_cluster
        % Bar colors different if filtering is taking place
        barh(mySamplesMetricGood(fliplr(plotSampleOrder)),'FaceColor',rgb('LightSkyBlue'),'EdgeColor',rgb('LightSkyBlue'));
        barh(mySamplesMetricBad(fliplr(plotSampleOrder)),'FaceColor',rgb('Violet'),'EdgeColor',rgb('Violet'));
    else
        % Bar colors green if filtering is not taking place (ie not enough
        % positions)
        barh(mySamplesMetricGood(fliplr(plotSampleOrder)),'FaceColor',rgb('LightGreen'),'EdgeColor',rgb('LightSkyBlue'));
        barh(mySamplesMetricBad(fliplr(plotSampleOrder)),'FaceColor',rgb('LightGreen'),'EdgeColor',rgb('Violet'));
    end
    xlim([.5 1.05])
    %set(gca, 'XTick',1:1:numel(plotSampleNames),'XTickLabel',plotSampleNames(plotSampleOrder),'XTickLabelRotation',90)
    set(gca, 'YTick',1:1:numel(plotSampleNames),'YTickLabel',plotSampleNames(fliplr(plotSampleOrder)),'TickLabelInterpreter','none')
    %xlabel('Sample Name', 'FontSize', 20) % x-axis label
    %ylabel('Fraction N Despite Good Coverage', 'FontSize', 20) % y-axis label
    ylabel('Sample Name', 'FontSize', 20) % x-axis label
    xlabel([ 'Mean MAF when Coverage > ' num2str(filtering_min_cov_to_consider_pos) ], 'FontSize', 20) % y-axis label
    title(['Filtering: ' mySamplesetName],'Interpreter','none') % title % Interpreter -> none to avoid intepreting underscore as subscript
    set(gca, 'FontSize', 14, 'FontName', 'Verdana')
    line([filtering_mean_maf_per_sample filtering_mean_maf_per_sample],ylim,'Color','k','LineStyle','--','LineWidth',1)
    %line([thresholdSketchy thresholdSketchy],ylim,'Color','k','LineStyle','--','LineWidth',1)
    % Record number of important positions on plot (top right)
    txt_position_label = [ '#pos = ' num2str(numel(p_important)) ];
    text(0.55,0.95*numel(plotSampleNames),txt_position_label,'FontSize',12,'HorizontalAlignment','left')
    txt_filter_label = [ '#failed = ' num2str(numel(mySubjectSamplesToRemove)) '/' num2str(numel(metric_bad)) ];
    text(0.55,0.9*numel(plotSampleNames),txt_filter_label,'FontSize',12,'HorizontalAlignment','left')
    hold off
    
    %figure(3)
    %hist(mySamplesFracSketchy)

    % Save figure
    if length(plotSampleNames) > 100 % bigger plot if many samples
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [6 32]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 6 32]);
    else
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [6 18]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 6 18]);
    end
    print([pwd '/' cluster_filtering_directory_name '/' 'Figure_' mySamplesetName '_Filtering.png'],'-dpng')
end

%%

% Monster figure for looking at different filter options...

if figures_boolean

    metric_bad = (mySamplesMeanMAF < filtering_mean_maf_per_sample);
    metric_good = (mySamplesMeanMAF >= filtering_mean_maf_per_sample);
    
    figure(3) 
    clf(3)
    % Block: Current filter - mean MAF
    subplot( 5, 7, [ 1 2 3 8 9 10 15 16 17 22 23 24 29 30 31 ] )
        % Save subsets so that good/bad samples can be different colors
        mySamplesMetricGood = mySamplesMeanMAF;
        mySamplesMetricGood(metric_bad) = 0;
        mySamplesMetricBad = mySamplesMeanMAF;
        mySamplesMetricBad(metric_good) = 0;
        hold on
        box on
        if numel(p_important) >= min_pos_to_examine_cluster
            % Bar colors different if filtering is taking place
            barh(mySamplesMetricGood(fliplr(plotSampleOrder)),'FaceColor',rgb('LightSkyBlue'),'EdgeColor',rgb('LightSkyBlue'));
            barh(mySamplesMetricBad(fliplr(plotSampleOrder)),'FaceColor',rgb('Violet'),'EdgeColor',rgb('Violet'));
        else
            % Bar colors green if filtering is not taking place (ie not enough
            % positions)
            barh(mySamplesMetricGood(fliplr(plotSampleOrder)),'FaceColor',rgb('LightGreen'),'EdgeColor',rgb('LightSkyBlue'));
            barh(mySamplesMetricBad(fliplr(plotSampleOrder)),'FaceColor',rgb('LightGreen'),'EdgeColor',rgb('Violet'));
        end
        xlim([.5 1.05])
        set(gca, 'YTick',1:1:numel(plotSampleNames),'YTickLabel',plotSampleNames(fliplr(plotSampleOrder)),'TickLabelInterpreter','none')
        ylabel('Sample Name', 'FontSize', 14) % x-axis label
        xlabel([ 'Mean MAF when cov > ' num2str(filtering_min_cov_to_consider_pos) ], 'FontSize', 20) % y-axis label
        title(['Filtering: Mean MAF when Good Cov'],'Interpreter','none') % title % Interpreter -> none to avoid intepreting underscore as subscript
        set(gca, 'FontSize', 8, 'FontName', 'Verdana')
        line([filtering_mean_maf_per_sample filtering_mean_maf_per_sample],ylim,'Color','k','LineStyle','--','LineWidth',1)
        % Record number of important positions on plot (top right)
        txt_position_label = [ '#pos = ' num2str(numel(p_important)) ];
        text(0.55,0.95*numel(plotSampleNames),txt_position_label,'FontSize',12,'HorizontalAlignment','left')
        txt_filter_label = [ '#failed = ' num2str(sum(metric_bad)) '/' num2str(numel(metric_bad)) ];
        text(0.55,0.9*numel(plotSampleNames),txt_filter_label,'FontSize',12,'HorizontalAlignment','left')
        hold off
    % Block: Alternative filter -  Mean Qual
    subplot( 5, 7, [ 4 11 18 25 32 ] )
        % Save subsets so that good/bad samples can be different colors
        thresholdAltQual = 75;
        mySamplesMetricGood = mySamplesMeanQual;
        mySamplesMetricGood(mySamplesMeanQual<=thresholdAltQual) = 0;
        mySamplesMetricBad = mySamplesMeanQual;
        mySamplesMetricBad(mySamplesMeanQual>thresholdAltQual) = 0;
        hold on
        box on
        if numel(p_important) >= min_pos_to_examine_cluster
            % Bar colors different if filtering is taking place
            barh(mySamplesMetricGood(fliplr(plotSampleOrder)),'FaceColor',rgb('LightSkyBlue'),'EdgeColor',rgb('LightSkyBlue'));
            barh(mySamplesMetricBad(fliplr(plotSampleOrder)),'FaceColor',rgb('Violet'),'EdgeColor',rgb('Violet'));
        else
            % Bar colors green if filtering is not taking place (ie not enough
            % positions)
            barh(mySamplesMetricGood(fliplr(plotSampleOrder)),'FaceColor',rgb('LightGreen'),'EdgeColor',rgb('LightSkyBlue'));
            barh(mySamplesMetricBad(fliplr(plotSampleOrder)),'FaceColor',rgb('LightGreen'),'EdgeColor',rgb('Violet'));
        end
        xlabel('Mean Qual', 'FontSize', 14) % y-axis label
        yticks([])
        title(['(Alt): Mean Qual'],'Interpreter','none') % title % Interpreter -> none to avoid intepreting underscore as subscript
        set(gca, 'FontSize', 8, 'FontName', 'Verdana')
        line([thresholdAltQual thresholdAltQual],ylim,'Color',rgb('DarkGray'),'LineStyle','--','LineWidth',1)
        txt_filter_label = [ '#failed = ' num2str(sum(mySamplesMeanQual<=75)) '/' num2str(numel(mySamplesFracSketchy)) ];
        text(0.1*max(xlim),0.925*numel(plotSampleNames),txt_filter_label,'FontSize',8,'Color','k','HorizontalAlignment','left')
        hold off
    % Block: Old filter - Frac N
    subplot( 5, 7, [ 5 12 19 26 33 ] )
        % Save subsets so that good/bad samples can be different colors
        thresholdAltFracN = 0.3;
        mySamplesMetricGood = mySamplesFracSketchy;
        mySamplesMetricGood(mySamplesFracSketchy > thresholdAltFracN) = 0;
        mySamplesMetricBad = mySamplesFracSketchy;
        mySamplesMetricBad(mySamplesFracSketchy <= thresholdAltFracN) = 0;
        hold on
        box on
        if numel(p_important) >= min_pos_to_examine_cluster
            % Bar colors different if filtering is taking place
            barh(mySamplesMetricGood(fliplr(plotSampleOrder)),'FaceColor',rgb('LightSkyBlue'),'EdgeColor',rgb('LightSkyBlue'));
            barh(mySamplesMetricBad(fliplr(plotSampleOrder)),'FaceColor',rgb('Violet'),'EdgeColor',rgb('Violet'));
        else
            % Bar colors green if filtering is not taking place (ie not enough
            % positions)
            barh(mySamplesMetricGood(fliplr(plotSampleOrder)),'FaceColor',rgb('LightGreen'),'EdgeColor',rgb('LightSkyBlue'));
            barh(mySamplesMetricBad(fliplr(plotSampleOrder)),'FaceColor',rgb('LightGreen'),'EdgeColor',rgb('Violet'));
        end
        xticks( 0.5:0.5:1 )
        xlim([0 1.05])
        xlabel('Frac N', 'FontSize', 14) % y-axis label
        yticks([])
        title(['(Old): Frac N'],'Interpreter','none') % title % Interpreter -> none to avoid intepreting underscore as subscript
        set(gca, 'FontSize', 8, 'FontName', 'Verdana')
        line([thresholdAltFracN thresholdAltFracN],ylim,'Color',rgb('DarkGray'),'LineStyle','--','LineWidth',1)
        txt_filter_label = [ '#failed = ' num2str(sum(mySamplesFracSketchy > 0.3)) '/' num2str(numel(mySamplesFracSketchy)) ];
        text(0.1*max(xlim),0.925*numel(plotSampleNames),txt_filter_label,'FontSize',8,'Color','k','HorizontalAlignment','left')
        hold off
    % Block: Histogram - Mean MAF 
    subplot( 5, 7, [ 6 7 ] )
        hold on
        box on
        histogram( mySamplesMeanMAF, 0.5:0.025:1, 'FaceColor', rgb('LightGray') )
        title(['Hist: Mean MAF (y:zoom)'],'Interpreter','none') % title % Interpreter -> none to avoid intepreting underscore as subscript
        xlim([0.5-.025/2 1+.025/2])
        xlabel('Mean MAF', 'FontSize', 14) % y-axis label
        ylim([0 sum(mySamplesMeanMAF<=filtering_mean_maf_per_sample)+1])
        ylabel('Num Samples', 'FontSize', 14)
        line([filtering_mean_maf_per_sample filtering_mean_maf_per_sample],ylim,'Color','k','LineStyle','--','LineWidth',1)
        set(gca, 'FontSize', 8, 'FontName', 'Verdana')
        hold off
    % Block: Scatter - Frac sketchy & Mean MAF
    subplot( 5, 7, [ 13 14 20 21 ] )
        hold on
        box on
        scatter( mySamplesFracSketchy(metric_good), mySamplesMeanMAF(metric_good), 50, rgb('Gray') )
        scatter( mySamplesFracSketchy(metric_bad), mySamplesMeanMAF(metric_bad), 50, rgb('LightGray') )
        title(['Scatter: Mean MAF vs Frac N'],'Interpreter','none') % title % Interpreter -> none to avoid intepreting underscore as subscript
        xlim([0 1])
        xlabel('Frac N')
        ylim([.5 1])
        ylabel('Mean MAF')
        line(xlim,[filtering_mean_maf_per_sample filtering_mean_maf_per_sample],'Color','k','LineStyle','--','LineWidth',1)
        line([.3 .3],ylim,'Color',rgb('LightGray'),'LineStyle','--','LineWidth',1)
        hold off
    % Block: Scatter: - Frac sketchy and Mean FQ
    subplot( 5, 7, [ 27 28 34 35 ] )
        hold on
        box on
        scatter( mySamplesFracSketchy(metric_good), mySamplesMeanQual(metric_good), 50, rgb('Gray') )
        scatter( mySamplesFracSketchy(metric_bad), mySamplesMeanQual(metric_bad), 50, rgb('LightGray') )
        title(['Scatter: Mean Qual vs Frac N'],'Interpreter','none') % title % Interpreter -> none to avoid intepreting underscore as subscript
        xlim([0 1])
        xlabel('Frac N')
        ylim([0 max(mySamplesMeanQual)])
        ylabel('Mean Qual')
        line(xlim,[75 75],'Color','k','LineStyle','--','LineWidth',1)
        line([.3 .3],ylim,'Color',rgb('LightGray'),'LineStyle','--','LineWidth',1)
        hold off
    % Main
    suptitle(['Filtering: ' mySamplesetName])

    % Save figure
    if length(plotSampleNames) > 100 % bigger plot if many samples
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [15 10]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 15 10]);
    else
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [15 10]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 15 10]);
    end
    print([pwd '/' cluster_filtering_directory_name '/' 'Figure_' mySamplesetName '_Filtering_Types.png'],'-dpng')

end



end