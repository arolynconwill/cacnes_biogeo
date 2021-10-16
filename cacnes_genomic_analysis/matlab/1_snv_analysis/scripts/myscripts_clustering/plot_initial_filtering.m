function plot_initial_filtering( dir_initialfilt, ...
    cacnes_frac_unfiltered, min_frac_cacnes_bracken, goodsamples_bracken, ...
    coverage_unfiltered_median, min_median_coverage_to_include_sample, goodsamples_coveragefilter, ...
    maf_unfiltered_fraclowmaf, max_low_maf_positions_to_include_sample, min_maf_for_purity, min_cov_for_purity, goodsamples_purityfilter, ...
    goodsamples_manual, goodsamples_all )


%% Bracken filter

figure(10)
clf(10)
hold on
box on
subplot( 2, 1, 1 )
histogram(cacnes_frac_unfiltered,0:.01:1)
ylim( [0 1200] )
line([ min_frac_cacnes_bracken, min_frac_cacnes_bracken],ylim,...
    'Color','r')
text(0.8*max(xlim),max(ylim)*.9,['min_frac_cacnes_bracken = ' num2str(min_frac_cacnes_bracken)],...
    'Interpreter','none','Color','r','HorizontalAlignment','right')
xlabel('Fraction reads assigned to C. acnes by bracken')
ylabel('Number of samples')
title('Initial Filter: Bracken')
set(gca,'FontSize',14)
subplot( 2, 1, 2 )
histogram(cacnes_frac_unfiltered,0:.01:1)
ylim([0 50])
line([ min_frac_cacnes_bracken, min_frac_cacnes_bracken],ylim,...
    'Color','r')
text(0.8*max(xlim),max(ylim)*.9,['min_frac_cacnes_bracken = ' num2str(min_frac_cacnes_bracken)],...
    'Interpreter','none','Color','r','HorizontalAlignment','right')
xlabel('Fraction reads assigned to C. acnes by bracken')
ylabel('Number of samples')
title('(y:zoom)')
set(gca,'FontSize',14)
hold off

% Save
print([dir_initialfilt '/' 'Filter_Hist_Bracken'],'-dpng')


%% Coverage filter

figure(11)
clf(11)
hold on
box on
subplot( 2, 1, 1 )
histogram(coverage_unfiltered_median(goodsamples_bracken),0:2.5:max(coverage_unfiltered_median))
line([ min_median_coverage_to_include_sample, min_median_coverage_to_include_sample],ylim,...
    'Color','r')
text(0.05*max(xlim),0.9*max(ylim),[ 'min_median_coverage_to_include_sample = ' num2str(min_median_coverage_to_include_sample)],... 
    'Interpreter','none','Color','r','HorizontalAlignment','left')
xlabel('Mean coverage of C. acnes reference genome')
ylabel('Number of samples')
title('Initial Filter: Coverage+Bracken')
set(gca,'FontSize',14)
subplot( 2, 1, 2 )
histogram(coverage_unfiltered_median(goodsamples_bracken),0:2.5:max(coverage_unfiltered_median))
ylim([0 40])
line([ min_median_coverage_to_include_sample, min_median_coverage_to_include_sample],ylim,...
    'Color','r')
text(0.05*max(xlim),0.9*max(ylim),[ 'min_median_coverage_to_include_sample = ' num2str(min_median_coverage_to_include_sample)],... 
    'Interpreter','none','Color','r','HorizontalAlignment','left')
xlabel('Mean coverage of C. acnes reference genome')
ylabel('Number of samples')
title('(y:zoom)')
set(gca,'FontSize',14)
hold off

% Save
print([dir_initialfilt '/' 'Filter_Hist_Coverage+Bracken'],'-dpng')


%% Purity filter
figure(12)
clf(12)
hold on
box on
subplot( 2, 1, 1 )
histogram(maf_unfiltered_fraclowmaf(goodsamples_bracken & goodsamples_coveragefilter),...
    0:max_low_maf_positions_to_include_sample/5:max(maf_unfiltered_fraclowmaf))
line([ max_low_maf_positions_to_include_sample, max_low_maf_positions_to_include_sample],ylim,...
    'Color','r')
text(0.15*max(xlim),0.9*max(ylim),[ 'max_low_maf_positions_to_include_sample = ' num2str(max_low_maf_positions_to_include_sample)],... 
    'Interpreter','none','Color','r','HorizontalAlignment','left')
text(0.15*max(xlim),0.85*max(ylim),[ '(min_maf_for_purity = ' num2str(min_maf_for_purity) ...
    ' with min_cov_for_purity = ' num2str(min_cov_for_purity) ')'],... 
    'Interpreter','none','Color','r','HorizontalAlignment','left')
xlabel('Fraction of positions with a low maf')
ylabel('Number of samples')
title('Initial Filter: Purity+Coverage+Bracken')
set(gca,'FontSize',14)
subplot( 2, 1, 2 )
histogram(maf_unfiltered_fraclowmaf(goodsamples_bracken & goodsamples_coveragefilter),...
    0:max_low_maf_positions_to_include_sample/5:max(maf_unfiltered_fraclowmaf))
ylim([0 60])
line([ max_low_maf_positions_to_include_sample, max_low_maf_positions_to_include_sample],ylim,...
    'Color','r')
text(0.15*max(xlim),0.9*max(ylim),[ 'max_low_maf_positions_to_include_sample = ' num2str(max_low_maf_positions_to_include_sample)],... 
    'Interpreter','none','Color','r','HorizontalAlignment','left')
text(0.15*max(xlim),0.85*max(ylim),[ '(min_maf_for_purity = ' num2str(min_maf_for_purity) ...
    ' with min_cov_for_purity = ' num2str(min_cov_for_purity) ')'],... 
    'Interpreter','none','Color','r','HorizontalAlignment','left')
xlabel('Fraction of positions with a low maf')
ylabel('y:zoom')
title('Initial Filter: Purity+Coverage+Bracken')
set(gca,'FontSize',14)
hold off

% Save
print([dir_initialfilt '/' 'Filter_Hist_Purity+Coverage+Bracken'],'-dpng')


%% Close figures

close(10)
close(11)
close(12)


%% Report on filtering

fprintf(1,[ 'Samples that passed bracken filter: ' ...
    num2str(sum(goodsamples_bracken)) '/' num2str(numel(goodsamples_bracken)) ...
    ' (' num2str(100*sum(goodsamples_bracken)/numel(goodsamples_bracken)) ' percent).\n' ] )
fprintf(1,[ 'Samples that passed coverage filter: ' ...
    num2str(sum(goodsamples_coveragefilter)) '/' num2str(numel(goodsamples_coveragefilter)) ...
    ' (' num2str(100*sum(goodsamples_coveragefilter)/numel(goodsamples_coveragefilter)) ' percent).\n' ] )
fprintf(1,[ 'Samples that passed BOTH bracken and coverage filters: ' ...
    num2str(sum(goodsamples_coveragefilter & goodsamples_bracken)) '/' num2str(numel(goodsamples_coveragefilter)) ...
    ' (' num2str(100*sum(goodsamples_coveragefilter & goodsamples_bracken)/numel(goodsamples_coveragefilter)) ' percent).\n' ] )
fprintf(1,[ 'Samples that passed purity filter: ' ...
    num2str(sum(goodsamples_purityfilter)) '/' num2str(numel(goodsamples_purityfilter)) ...
    ' (' num2str(100*sum(goodsamples_purityfilter)/numel(goodsamples_purityfilter)) ' percent).\n' ] )
fprintf(1,[ 'Samples removed with contaminated samples filter: ' ...
    num2str(sum(~goodsamples_manual)) '/' num2str(numel(goodsamples_all)) ...
    ' (' num2str(100*sum(~goodsamples_manual)/numel(goodsamples_all)) ' percent).\n' ] )
fprintf(1,[ 'Total samples that passed ALL filtering so far: ' ...
    num2str(sum(goodsamples_all)) '/' num2str(numel(goodsamples_all)) ...
    ' (' num2str(100*sum(goodsamples_all)/numel(goodsamples_all)) ' percent).\n' ] )


end